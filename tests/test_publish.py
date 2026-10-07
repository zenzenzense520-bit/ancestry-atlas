import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from check_publish import allowed
class PublishGuard(unittest.TestCase):
    def test_rejects_personal_data_and_derived_artifacts(self):
        for path in ['input/sample.txt','results/global_pca.json','results/genotypes.npz','qa/report.png','references/panel.psam','archive.zip']:
            self.assertFalse(allowed(path),path)
    def test_rejects_environment_files(self):
        self.assertFalse(allowed('.env'));self.assertFalse(allowed('src/.env.local'))
    def test_accepts_code_docs_and_licensed_fonts(self):
        for path in ['src/qc.py','src/report.body.html','src/fonts/atlas.woff2','src/fonts/OFL.txt','docs/METHODS.md']:
            self.assertTrue(allowed(path),path)
