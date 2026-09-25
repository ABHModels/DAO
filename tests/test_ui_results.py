"""Results API regression tests using temporary, independent run fixtures."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('dao_ui', ROOT / 'ui.py')
ui = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ui)


class ResultsTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.previous = ui.WORK_DIR
        ui.WORK_DIR = self.tmp.name
        self.addCleanup(setattr, ui, 'WORK_DIR', self.previous)
        self.client = ui.app.test_client()
        self.results = Path(self.tmp.name) / 'results'

    def run_fixture(self, run_id='abc', label='My run'):
        path = self.results / run_id
        path.mkdir(parents=True)
        meta = dict(hash='abc', label=label, corona='cutoffpl', Gamma=2,
                    E_cut=300, nh=15, zeta=3, frac=-1, angsca=True)
        (path / 'params.json').write_text(json.dumps(meta))
        for it, temperatures in [(1, [100., 200.]), (2, [100.01, 200.02])]:
            (path / f'profile_iter{it:03d}.dat').write_text(
                '# depth tau T ne x y xi\n' + ''.join(
                    f'{i} {0.1 * (i+1)} {t} 1e15 0 0 3\n'
                    for i, t in enumerate(temperatures)))
            (path / f'emergent_iter{it:03d}.dat').write_text(
                '# I(mu=-0.5000)\n# I(mu=0.5000)\n'
                '100 1 2 3 4\n1000 10 20 30 40\n')
            (path / f'moments_iter{it:03d}.dat').write_text(
                '# E depth tau T J0 J2 J3\n100 0 .1 100 7 0 0\n1000 0 .1 100 8 0 0\n')
        (path / 'line_labels_iter001.dat').write_text('"Fe 26" 0 6400 0 0 4 F\n')
        return path

    def test_pages_preserve_theme_and_add_navigation(self):
        for route in ['/', '/docs', '/plots', '/compare', '/convergence']:
            response = self.client.get(route)
            self.assertEqual(response.status_code, 200)
            self.assertIn('--bg:#f5f2ed', response.text)
            self.assertIn('href="/compare"', response.text)
        self.assertIn('runLabel', self.client.get('/').text)

    def test_campaigns_same_hash_and_run_cards(self):
        self.run_fixture()
        self.run_fixture('campaign #1/abc', '<custom text>')
        runs = self.client.get('/api/runs').json['runs']
        self.assertEqual([r['id'] for r in runs], ['abc', 'campaign #1/abc'])
        run = runs[1]
        self.assertEqual(run['group'], 'campaign #1')
        self.assertIn('<custom text>', run['label'])
        self.assertIn('campaign #1/abc', run['title'])
        self.assertEqual(run['iters'], [1, 2])
        self.assertEqual(run['convergence']['status'], 'converged')
        self.assertAlmostEqual(run['convergence']['dT'], .0001)
        self.assertTrue(run['groups'])
        rid = quote(run['id'], safe='')
        data = self.client.get(f'/api/data/{rid}/2').json
        self.assertEqual(data['emergent']['angles']['0.5000'], [4, 40])
        self.assertEqual(data['profile']['T_K'], [100.01, 200.02])
        self.assertEqual(data['mean_intensity']['J0'], [7, 8])
        self.assertEqual(set(self.client.get('/api/profiles/'+rid).json['profiles']), {'1','2'})
        self.assertEqual(set(self.client.get('/api/mean_intensity/'+rid).json['moments']), {'1','2'})
        self.assertEqual(self.client.get('/api/line_labels/'+rid).json['lines'][0]['label'], 'Fe 26')

    def test_live_discovery_and_partial_convergence(self):
        self.assertEqual(self.client.get('/api/runs').json['runs'], [])
        path = self.run_fixture()
        self.assertEqual(len(self.client.get('/api/runs').json['runs']), 1)
        (path / 'profile_iter003.dat').write_text('0 .1 100.01 1e15 0 0 3\n')
        run = self.client.get('/api/runs').json['runs'][0]
        self.assertEqual(run['iters'], [1,2,3])
        self.assertIsNone(run['convergence']['dT'])
        self.assertNotEqual(run['convergence']['status'], 'converged')

    def test_convergence_orders_iterations_numerically(self):
        path = self.run_fixture()
        (path / 'profile_iter999.dat').write_text('0 .1 100 1e15 0 0 3\n')
        (path / 'profile_iter1000.dat').write_text('0 .1 110 1e15 0 0 3\n')
        self.assertAlmostEqual(ui._run_convergence(path)['dT'], .1)

    def test_exact_fe_iteration_and_path_containment(self):
        path = self.run_fixture()
        (path / 'fe_lines.dat').write_text(
            '# Fe lines iter=10\nWrong 0 6.5 1 0\n'
            '# Fe lines iter=1\nRight 0 6.4 1 0\n')
        lines = self.client.get('/api/data/abc/1').json['fe_lines']
        self.assertEqual([line['label'] for line in lines], ['Right'])
        (self.results / 'outside').symlink_to(Path(self.tmp.name))
        for endpoint in ['data/outside/1', 'profiles/outside',
                         'mean_intensity/../outside', 'line_labels/outside']:
            self.assertEqual(self.client.get('/api/'+endpoint).status_code, 404)


if __name__ == '__main__':
    unittest.main()
