"""Launch the isolated all-P=1 full thermal calculation, retaining provenance."""
from pathlib import Path
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
info = json.loads((HERE / 'build_manifest.json').read_text())
RUN = ROOT / info['run_directory']
EXE = HERE / 'maindaocl_all_P1'
assert not RUN.exists(), 'Refuse to replace an existing calculation'
assert not (HERE / 'launch.json').exists(), 'Experiment already launched'
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
assert all(digest(ROOT / p) == h for p, h in info['production_sha256'].items())
baseline_hashes = {str(p.relative_to(ROOT)): digest(p) for p in (ROOT/'results/5b7fe321').rglob('*') if p.is_file()}
env = os.environ.copy(); env.update(info['environment'])
command = [str(EXE), '-corona','nthcomp','-Gamma','1.8','-kT_e','60','-kT_bb','0.01',
           '-nh','18','-zeta','3','-frac','-1','-incidence','0.7071','-angsca','0',
           '-reflionx_norm','1','-log_depth','1','-O','1','-Fe','1',
           '-label','All lines P=1: self-consistent atmosphere, lognH=18 Gamma=1.8 logxi=3']
record = dict(info, command=command, cwd=str(ROOT), experiment_executable_sha256=digest(EXE),
              started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
              baseline_input_sha256=baseline_hashes, state='starting', wrapper_pid=os.getpid())

def save():
    (HERE/'launch.json').write_text(json.dumps(record,indent=2)+'\n')
    if RUN.exists():
        (RUN/'launch_provenance.json').write_text(json.dumps(record,indent=2)+'\n')

start=time.monotonic()
try:
    save()
    test=subprocess.run([str(EXE),'--verify-no-survival'],cwd=ROOT,env=env,
                         stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,check=True)
    (HERE/'verification.log').write_text(test.stdout)
    print(test.stdout,flush=True)
    with (HERE/'console.log').open('w') as log:
        proc=subprocess.Popen(command,cwd=ROOT,env=env,stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT,text=True,bufsize=1)
        record.update(pid=proc.pid,state='running');save()
        for line in proc.stdout:
            log.write(line);log.flush();print(line,end='',flush=True)
            if RUN.exists() and not (RUN/'launch_provenance.json').exists(): save()
        code=proc.wait()
    record['return_code']=code
    status=json.loads((RUN/'thermal_status.json').read_text())
    record['thermal_status']=status
    assert code==0 and status['state']=='converged', f'Thermal solve incomplete: {status}'
    record['state']='converged'
except Exception as exc:
    record.update(state='failed',error=str(exc))
    raise
finally:
    record['wall_seconds']=time.monotonic()-start
    record['production_files_unchanged']=all(digest(ROOT/p)==h for p,h in info['production_sha256'].items())
    record['baseline_files_unchanged']=all(digest(ROOT/p)==h for p,h in baseline_hashes.items())
    save()
    assert record['production_files_unchanged'] and record['baseline_files_unchanged']
    print(json.dumps({k:record[k] for k in ('state','run_directory','wall_seconds')},indent=2),flush=True)
