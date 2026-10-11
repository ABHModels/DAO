"""Standalone metadata regression; no Cloudy or HEASoft runtime required."""
import json
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest
import os
import re

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

    def run_model(self, args, exe=None):
        result = subprocess.run([str(exe or self.exe), *args], cwd=self.work,
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

    def test_ion_output_switches_do_not_change_physics_hash(self):
        args = ['-corona', 'powerlaw', '-Gamma', '2']
        baseline, _ = self.run_model(args)
        self.assertFalse(baseline['save_oxygen'])
        self.assertFalse(baseline['save_iron'])
        for flags, oxygen, iron in [(['-O'], True, False), (['-Fe'], False, True),
                                    (['-O', '1', '-Fe', 'yes'], True, True),
                                    (['-O', 'true', '-Fe', '0'], True, False),
                                    (['-O', 'no', '-Fe', 'false'], False, False)]:
            meta, _ = self.run_model(args + flags)
            self.assertEqual(meta['hash'], baseline['hash'])
            self.assertEqual(meta['save_oxygen'], oxygen)
            self.assertEqual(meta['save_iron'], iron)
        for flag in ('-O', '-Fe'):
            bad = subprocess.run([str(self.exe), *args, flag, 'bad'],
                                 cwd=self.work, capture_output=True, text=True)
            self.assertNotEqual(bad.returncode, 0)
            self.assertIn('accepts an optional', bad.stderr)
            synthetic = subprocess.run([str(self.exe), *args, '-test_rt', 'compps',
                                        '-kT_e', '60', flag],
                                       cwd=self.work, capture_output=True, text=True)
            self.assertNotEqual(synthetic.returncode, 0)
            self.assertIn('require production mode', synthetic.stderr)

    def test_fixed_reflionx_normalization(self):
        args = ['-corona', 'cutoffpl', '-Gamma', '2', '-Ecut', '300']
        meta, summary = self.run_model(args)
        self.assertEqual(meta['illumination_normalization'], 'reflionx')
        self.assertNotIn('xillver_norm', meta)
        self.assertNotIn('reflionx_norm', meta)
        self.assertIn('Normalization: reflionx flux', summary)
        self.assertIn('incident_corona_flux=7.9577e+16', summary)
        isotropic, _ = self.run_model(args + ['-incidence', '-2'])
        self.assertEqual(isotropic['incidence'], -2)
        self.assertEqual(isotropic['illumination_normalization'], 'reflionx')
        self.assertNotEqual(meta['hash'], isotropic['hash'])

    def test_removed_switches_rejected_before_creating_results(self):
        args = ['-corona', 'cutoffpl', '-Gamma', '2', '-Ecut', '300']
        for flags in (['-xillver_norm'], ['-xillver_norm', '0'],
                      ['-xillver_norm', '1'], ['-reflionx_norm'],
                      ['-reflionx_norm', '0'], ['-reflionx_norm', '1'],
                      ['-sc'], ['-sc', 'bezier3'], ['-rt_thermal_balance'],
                      ['-rt_thermal_balance', '0'], ['-rt_thermal_balance', '1'],
                      ['-rt_thermal_balance', 'true'], ['-rt_thermal_balance', 'yes']):
            with self.subTest(flags=flags), tempfile.TemporaryDirectory() as work:
                result = subprocess.run([str(self.exe), *args, *flags], cwd=work,
                                        capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('unknown or incomplete option', result.stderr)
                self.assertIn(flags[0], result.stderr)
                self.assertFalse((Path(work) / 'results').exists())

    def test_depth_grid_flag_and_hash(self):
        args = ['-corona', 'cutoffpl', '-Gamma', '2', '-Ecut', '300']
        baseline, summary = self.run_model(args)
        self.assertTrue(baseline['log_depth'])
        self.assertEqual(baseline['depth_cells'], 48)
        self.assertIn('48 cells', summary)
        self.assertIn('Depth grid: log', summary)
        for on in ('1', 'true', 'yes'):
            meta, _ = self.run_model(args + ['-log_depth', on])
            self.assertTrue(meta['log_depth'])
            self.assertEqual(meta['hash'], baseline['hash'])
        tanh_hash = None
        for off in ('0', 'false', 'no'):
            meta, summary = self.run_model(args + ['-log_depth', off])
            self.assertFalse(meta['log_depth'])
            self.assertIn('Depth grid: tanh (k=3)', summary)
            self.assertNotEqual(meta['hash'], baseline['hash'])
            if tanh_hash is not None:
                self.assertEqual(meta['hash'], tanh_hash)
            tanh_hash = meta['hash']

    def test_depth_resolution_separates_result_directories(self):
        # Build the same parser with a different grid default. Identical
        # physical parameters must not overwrite a different-resolution run.
        variant = self.work / 'grid80'
        variant.mkdir()
        (variant / 'params.cpp').write_text((ROOT / 'source/params.cpp').read_text())
        header = (ROOT / 'source/rt_grids.h').read_text()
        header, count = re.subn(r'ND_EDGE_DEFAULT\s*=\s*\d+',
                               'ND_EDGE_DEFAULT = 81', header)
        self.assertEqual(count, 1)
        (variant / 'rt_grids.h').write_text(header)
        exe = variant / 'params'
        subprocess.run(shlex.split(os.environ.get('CXX', 'c++')) + [
            '-std=c++17', '-Wall', '-I' + str(ROOT / 'source'),
            str(self.work / 'main.cpp'), str(variant / 'params.cpp'),
            '-o', str(exe)], check=True)
        args = ['-corona', 'cutoffpl', '-Gamma', '2', '-Ecut', '300']
        coarse, _ = self.run_model(args)
        fine, summary = self.run_model(args, exe=exe)
        self.assertEqual(fine['depth_cells'], 80)
        self.assertIn('80 cells', summary)
        self.assertNotEqual(coarse['hash'], fine['hash'])

    def test_free_electron_solver_separates_previous_density_hashes(self):
        args = ['-corona', 'nthcomp', '-Gamma', '2', '-kT_e', '60', '-kT_bb', '.01', '-angsca', '0']
        new, summary = self.run_model(args)
        self.assertTrue(new['rt_thermal_balance'])
        self.assertEqual(new['solver'], 'rt_energy_balance_v1')
        self.assertEqual(new['scattering_density'], 'cloudy_free_electrons')
        self.assertEqual(new['temperature_iteration'], 'local_cell_response_v3')
        self.assertEqual(new['column_energy_tolerance'], 0.01)
        self.assertEqual(new['initialization'], 'boundary_intensity_v1')
        self.assertNotIn('experimental', new)
        self.assertEqual(new['intensity_location'], 'cell-volume-average')
        self.assertEqual(new['E_rt_lo'], 1.0)
        self.assertEqual(new['E_rt_hi'], 1e6)
        self.assertIn('RT thermal balance', summary)
        self.assertNotEqual(new['hash'], '8ae4f853')  # former fixed 1.21*nH run
        # A new initialization/root scheme must never overwrite the stopped
        # Cloudy-ne cold-start benchmark with identical physical parameters.
        cold, _ = self.run_model(args + ['-zeta','1','-frac','-1','-incidence','.7071'])
        self.assertNotEqual(cold['hash'], '6467d7ec')

    def test_test_mode_uses_shared_transfer_and_new_hash(self):
        args = ['-corona','blackbody','-kT_bb','.1','-kT_e','60',
                '-test_rt','compps','-tau','.5']
        current, _ = self.run_model(args)
        self.assertEqual(current['transfer_solver'], 'constant_cell_v1')
        self.assertEqual(current['intensity_location'], 'cell-volume-average')
        self.assertEqual(current['angular_points'], 8)
        self.assertFalse(current['rt_thermal_balance'])  # prescribed T, shared RT
        # Reconstruct the previous test hash: the new run must be separate.
        variant = self.work / 'previous_test_params.cpp'
        source = (ROOT / 'source/params.cpp').read_text()
        source, count = re.subn(r'^.*if \(p.test_rt\) hash_input \+= "\|transfer=constant_cell_v1.*\n',
                                '', source, flags=re.M)
        self.assertEqual(count, 1)
        variant.write_text(source)
        old_exe = self.work / 'previous_test_params'
        subprocess.run(shlex.split(os.environ.get('CXX', 'c++')) + [
            '-std=c++17','-I'+str(ROOT/'source'),str(self.work/'main.cpp'),
            str(variant),'-o',str(old_exe)], check=True)
        previous, _ = self.run_model(args, exe=old_exe)
        self.assertNotEqual(current['hash'], previous['hash'])

    def test_angular_resolution_separates_result_directories(self):
        exe = self.work / 'params10'
        subprocess.run(shlex.split(os.environ.get('CXX', 'c++')) + [
            '-std=c++17','-DDAO_RT_ANGLES=10','-I'+str(ROOT/'source'),
            str(self.work/'main.cpp'),str(ROOT/'source/params.cpp'),
            '-o',str(exe)], check=True)
        base = ['-corona','blackbody','-kT_bb','.1','-kT_e','60']
        for args in (base, base+['-test_rt','compps']):
            eight, _ = self.run_model(args)
            ten, _ = self.run_model(args, exe=exe)
            self.assertEqual(eight['angular_points'], 8)
            self.assertEqual(ten['angular_points'], 10)
            self.assertNotEqual(eight['hash'], ten['hash'])

    def test_both_scattering_modes_use_thermal_solver(self):
        args = ['-corona', 'nthcomp', '-Gamma', '2', '-kT_e', '60', '-kT_bb', '.01']
        default, _ = self.run_model(args)
        directional, _ = self.run_model(args + ['-angsca','1'])
        mean, _ = self.run_model(args + ['-angsca','0'])
        self.assertTrue(default['angsca'])
        self.assertEqual(default['hash'], directional['hash'])
        self.assertNotEqual(mean['hash'], directional['hash'])
        for meta in (default,directional,mean):
            self.assertEqual(meta['solver'], 'rt_energy_balance_v1')
        test, _ = self.run_model(args + ['-test_rt','compps'])
        self.assertFalse(test['rt_thermal_balance'])
        self.assertEqual(test['E_rt_lo'],10)
        result = subprocess.run([str(self.exe), *args, '-test_rt','compps',
                                 '-rt_thermal_balance','1'], cwd=self.work,
                                capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn("unknown or incomplete option '-rt_thermal_balance'",result.stderr)

    def test_invalid_scattering_flag(self):
        args = ['-corona','powerlaw','-Gamma','2']
        for flag in ('-angsca',):
            for values in ([],['maybe'],['2'],['-1']):
                result = subprocess.run([str(self.exe), *args, flag, *values],
                                        cwd=self.work,capture_output=True,text=True)
                self.assertNotEqual(result.returncode,0)
                self.assertIn(flag+' requires',result.stderr)

    def test_invalid_depth_grid_flag_rejected(self):
        args = ['-corona', 'cutoffpl', '-Gamma', '2', '-Ecut', '300', '-log_depth']
        for values in ([], ['maybe'], ['2'], ['-1']):
            result = subprocess.run([str(self.exe), *args, *values], cwd=self.work,
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('-log_depth requires', result.stderr)

if __name__ == '__main__':
    unittest.main()
