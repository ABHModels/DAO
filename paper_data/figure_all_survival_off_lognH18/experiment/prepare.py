"""Build an isolated, reproducible all-line P=1 thermal experiment.

Only copied sources are edited. Both line-source entry points use the same
modified probability function. The production solver and tolerances are kept.
"""
from pathlib import Path
import difflib
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SNAP = HERE / 'source'
SNAP.mkdir(exist_ok=True)
accepted = json.loads((ROOT / 'results/5b7fe321/launch_provenance.json').read_text())
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
assert all(digest(ROOT / p) == h for p, h in accepted['source_sha256'].items())

patches = []


def save_patch(name, edited):
    original = (ROOT / name).read_text()
    (SNAP / Path(name).name).write_text(edited)
    patches.extend(difflib.unified_diff(original.splitlines(True), edited.splitlines(True),
                   fromfile='original/' + name, tofile='experiment/' + name))


name = 'source/cloudy_interface_v2.cpp'
text = (ROOT / name).read_text()
start = text.index('LineProbabilities line_probabilities(')
body = text.index('{', start)
end = text.index('\n}\n}', body)
text = text[:body] + '''{
    // Experiment: every bound-bound line is injected with P=1 in both
    // the local temperature search and the full-column transfer solve.
    // No trapping-destruction power is generated in this artificial limit.
    // Diagnostic beta/Pelec/Pdest below are imposed values, not opacity fits.
    return {1.0, 0.0, 0.0, 1.0 + r.y, 1.0};
''' + text[end:]
save_patch(name, text)

name = 'source/params.cpp'
text = (ROOT / name).read_text()
text = text.replace('#include <sys/stat.h>', '#include <sys/stat.h>\n#include <filesystem>\n#include <stdexcept>')
assert text.count('unsigned int h = 2166136261u;') == 1
text = text.replace('unsigned int h = 2166136261u;',
                    'hash_input += "|line_survival=all_P1_experiment_v1";\n\t\tunsigned int h = 2166136261u;')
needle = '\t\tmkdir(p.run_dir, 0755);'
assert text.count(needle) == 1
text = text.replace(needle, '''        if (std::filesystem::exists(p.run_dir))
            throw std::runtime_error("Experiment refuses to overwrite an existing result directory");
''' + needle)
needle = '\t\t\tfprintf(fp, "{\\n");'
assert text.count(needle) == 1
text = text.replace(needle, needle + '\n' + r'''            fprintf(fp, "  \"line_survival\": \"all_P1_experiment_v1\",\n");
            fprintf(fp, "  \"survival_factor_enabled\": false,\n");
''')
save_patch(name, text)

name = 'maindaocl.cpp'
text = (ROOT / name).read_text()
text = text.replace('#include <chrono>', '#include <chrono>\n#include "verify_no_survival.h"')
text = text.replace('\ttry {', '''\ttry {
        if (argc == 2 && std::string(argv[1]) == "--verify-no-survival") {
            verify_no_survival();
            return 0;
        }''', 1)
needle = '\t\tstatic const char* EGRID_FILE = "cloudy_energy_grid.dat";'
assert text.count(needle) == 1
text = text.replace(needle, '''        const std::string energy_path = std::string(par.run_dir) + "/cloudy_energy_grid.dat";
        const char* EGRID_FILE = energy_path.c_str();''')
save_patch(name, text)
(HERE / 'experiment.patch').write_text(''.join(patches))

# The existing physics hash can be extended by continuing its FNV-1a state.
h = int('5b7fe321', 16)
for byte in b'|line_survival=all_P1_experiment_v1':
    h = ((h ^ byte) * 16777619) & 0xffffffff
run_hash = f'{h:08x}'
assert not (ROOT / 'results' / run_hash).exists()
manifest = dict(baseline='5b7fe321', run_hash=run_hash,
                run_directory=f'results/{run_hash}',
                production_sha256=accepted['source_sha256'],
                experimental_sha256={str(p.relative_to(ROOT)): digest(p) for p in SNAP.glob('*')},
                environment=accepted['environment'])
(HERE / 'build_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print('Experiment run hash:', run_hash)
