"""Standalone metadata regression; no Cloudy or HEASoft runtime required."""
import json
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest
import os

ROOT = Path(__file__).resolve().parents[1]

class RunMetadataTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix='dao-run-metadata-')
        cls.work = Path(cls.tmp.name)
        main = cls.work / 'main.cpp'
        main.write_text('#include "params.h"\nint main(int argc, char** argv) { read_params(argc, argv); }\n')
        cls.exe = cls.work / 'params'
        subprocess.run(shlex.split(os.environ.get('CXX', 'c++')) + [
            '-std=c++17', '-Wall', '-I' + str(ROOT / 'source'), str(main),
            str(ROOT / 'source/params.cpp'), '-o', str(cls.exe)], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_model(self, args):
        result = subprocess.run([str(self.exe), *args], cwd=self.work,
                                capture_output=True, text=True, check=True)
        run_hash = next(line.split()[2] for line in result.stdout.splitlines()
                        if line.startswith('Run hash:'))
        path = self.work / 'results' / run_hash
        return json.loads((path / 'params.json').read_text()), (path / 'RUN.txt').read_bytes().decode('utf-8')

    def test_label_round_trip_and_same_physics_hash(self):
        args = ['-corona', 'cutoffpl', '-Gamma', '2', '-Ecut', '300']
        baseline, summary = self.run_model(args)
        self.assertEqual(baseline['label'], '')
        self.assertIn('Label:   (none)', summary)
        labels = ['Iron test', 'He said "yes" \\ path', "O'Brien $(touch nope) `date` <script>",
                  'line1\nline2\t\r\b\f\x01', '铁线 α café', 'x' * 4096]
        for label in labels:
            with self.subTest(label=label[:40]):
                meta, summary = self.run_model(args + ['-label', label])
                self.assertEqual(meta['label'], label)
                self.assertEqual(meta['hash'], baseline['hash'])
                self.assertEqual({k:v for k,v in meta.items() if k not in ('label','time')},
                                 {k:v for k,v in baseline.items() if k not in ('label','time')})
                self.assertIn('Label:   ' + label + '\n', summary)
                self.assertIn('Time:    ' + meta['time'], summary)
                self.assertIn('Hash:    ' + meta['hash'], summary)

    def test_model_summaries(self):
        cases = [('powerlaw',['-Gamma','2']), ('nthcomp',['-Gamma','2','-kT_e','60','-kT_bb','.1']),
                 ('comptt',['-kT_e','50','-kT_bb','.05','-taup','1']), ('blackbody',['-kT_bb','.1'])]
        for model, args in cases:
            meta, summary = self.run_model(['-corona',model,*args])
            self.assertIn(model, summary)
            self.assertEqual(meta['corona'],model)
        _, summary = self.run_model(['-corona','blackbody','-kT_bb','.1','-kT_e','60','-test_rt','compps'])
        self.assertIn('compps test', summary)

    def test_missing_label_value_rejected(self):
        result = subprocess.run([str(self.exe),'-corona','powerlaw','-Gamma','2','-label'],
                                cwd=self.work,capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('-label requires text',result.stderr)

if __name__ == '__main__':
    unittest.main()
