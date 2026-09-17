"""Regression checks for R1-M1 using temporary copies of the real counts fixture.

Altered flags/count metadata and the one cloned fifth guide are test inputs only,
not additional biological data or scientific validation results.
"""
import json
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np
import pandas as pd
from scipy import sparse

from layered_eval.bundle import validate_bundle
from layered_eval.pipeline import prepare_counts, evaluate_measured


PACKAGE = Path(__file__).resolve().parents[1]
FIXTURE = PACKAGE / 'examples/public_microglia_fixture'


class SourceEligibilityRegression(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='layered-source-regression-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.input = self.root / 'input'
        shutil.copytree(FIXTURE, self.input)
        self.meta = pd.read_csv(self.input / 'cells.csv')
        self.control = sorted(self.meta.loc[self.meta.is_control, 'guide_id'].unique())[0]
        self.target = sorted(self.meta.loc[~self.meta.is_control, 'target_id'].unique())[0]
        self.treated = sorted(self.meta.loc[self.meta.target_id.eq(self.target), 'guide_id'].unique())[0]

    def write_meta(self, metadata=None, directory=None):
        (self.meta if metadata is None else metadata).to_csv(
            (self.input if directory is None else directory) / 'cells.csv', index=False)

    def prepare(self, name='bundle'):
        out = self.root / name
        prepare_counts(self.input / 'config.json', out)
        return out

    def test_excluded_control_cannot_fill_four_guide_requirement(self):
        self.meta.loc[self.meta.guide_id.eq(self.control), 'source_excluded_guide'] = True
        self.write_meta()
        with self.assertRaisesRegex(ValueError, 'four qualified control guides required after'):
            self.prepare()
        self.assertFalse((self.root / 'bundle').exists())

    def test_extra_excluded_control_never_enters_fit_or_reference(self):
        baseline = self.prepare('baseline')
        original = self.meta.copy()
        duplicate = self.meta[self.meta.guide_id.eq(self.control)].copy()
        duplicate.cell_id = duplicate.cell_id.astype(str) + '_test_only'
        duplicate.guide_id = 'excluded_control_test_only'
        duplicate.source_excluded_guide = True
        counts = sparse.load_npz(self.input / 'counts.npz')
        sparse.save_npz(self.input / 'counts.npz', sparse.vstack([
            counts, counts[original.guide_id.eq(self.control).to_numpy()]], format='csr'))
        self.write_meta(pd.concat([original, duplicate], ignore_index=True))
        repaired = self.prepare('with_excluded_control')
        a = json.loads((baseline / 'design.json').read_text())
        b = json.loads((repaired / 'design.json').read_text())
        for key in ('targets', 'target_groups', 'control_groups', 'pca_fit_indices'):
            self.assertEqual(a[key], b[key])
        self.assertNotIn('excluded_control_test_only', b['control_groups'])
        np.testing.assert_allclose(np.load(baseline / 'representations/PCA.npy'),
                                   np.load(repaired / 'representations/PCA.npy')[:len(original)], atol=1e-12)
        for path, label in ((baseline, 'base_scores'), (repaired, 'repaired_scores')):
            evaluate_measured(path, self.root / label, iterations=2, n=16)
        a = pd.read_csv(self.root / 'base_scores/summary.csv')
        b = pd.read_csv(self.root / 'repaired_scores/summary.csv')
        pd.testing.assert_frame_equal(a, b, check_exact=False, atol=1e-12, rtol=0)

    def test_excluded_treated_guide_removes_target_qualification(self):
        self.meta.loc[self.meta.guide_id.eq(self.treated), 'source_excluded_guide'] = True
        self.write_meta()
        out = self.prepare()
        d = json.loads((out / 'design.json').read_text())
        self.assertNotIn(self.target, d['targets'])
        self.assertEqual(len(d['targets']), 3)
        self.assertTrue(validate_bundle(out)['source_exclusion_checked'])

    def test_source_minimum_applies_to_controls_and_treated(self):
        config_path = self.input / 'config.json'
        config = json.loads(config_path.read_text())
        config['min_source_cells_per_guide'] = 60
        config_path.write_text(json.dumps(config))
        self.meta.source_guide_n = 100
        self.meta.loc[self.meta.guide_id.eq(self.control), 'source_guide_n'] = 52
        self.write_meta()
        with self.assertRaisesRegex(ValueError, 'min_source_cells_per_guide'):
            self.prepare('low_control')
        self.meta.source_guide_n = 100
        self.meta.loc[self.meta.guide_id.eq(self.treated), 'source_guide_n'] = 32
        self.write_meta()
        out = self.prepare('low_treated')
        design = json.loads((out / 'design.json').read_text())
        self.assertEqual(design['min_source_cells_per_guide'], 60)
        self.assertNotIn(self.target, design['targets'])
        self.assertEqual(validate_bundle(out)['source_min_cells_checked'], 60)

    def test_boolean_case_and_outer_whitespace_are_normalized(self):
        baseline = self.prepare('baseline')
        self.meta.is_control = self.meta.is_control.map({True: ' TRUE ', False: ' false '})
        self.meta.source_excluded_guide = ' FaLsE '
        self.write_meta()
        out = self.prepare('normalized')
        for name in ('control_groups', 'target_groups', 'pca_fit_indices'):
            self.assertEqual(json.loads((baseline / 'design.json').read_text())[name],
                             json.loads((out / 'design.json').read_text())[name])
        self.assertTrue(validate_bundle(out)['source_exclusion_checked'])

    def test_invalid_boolean_values_rejected_by_both_entrypoints(self):
        baseline = self.prepare('baseline')
        clean = self.meta.copy()
        for column in ('is_control', 'source_excluded_guide'):
            for bad in ('', 'yes', 'unknown', 0, 1, np.nan):
                with self.subTest(column=column, bad=bad):
                    modified = clean.copy()
                    modified[column] = bad
                    self.write_meta(modified)
                    with self.assertRaisesRegex(ValueError, column + ' requires explicit true/false'):
                        self.prepare('must_not_exist')
                    self.write_meta(modified, baseline)
                    with self.assertRaisesRegex(ValueError, column + ' requires explicit true/false'):
                        validate_bundle(baseline)
        self.assertFalse((self.root / 'must_not_exist').exists())

    def test_missing_flag_columns_are_actionable_errors(self):
        baseline = self.prepare('baseline')
        for column in ('is_control', 'source_excluded_guide'):
            with self.subTest(column=column):
                modified = self.meta.drop(columns=[column])
                self.write_meta(modified)
                self.write_meta(modified, baseline)
                for action in (lambda: self.prepare('bad'), lambda: validate_bundle(baseline)):
                    with self.assertRaisesRegex(ValueError, 'missing required metadata column.*' + column):
                        action()

    def test_per_guide_metadata_must_be_consistent(self):
        baseline = self.prepare('baseline')
        clean = self.meta.copy()
        row = self.meta.index[self.meta.guide_id.eq(self.control)][0]
        for column, value in [('is_control', False), ('source_excluded_guide', True),
                              ('source_guide_n', int(self.meta.loc[row, 'source_guide_n']) + 1)]:
            with self.subTest(column=column):
                changed = clean.copy()
                changed.loc[row, column] = value
                self.write_meta(changed)
                self.write_meta(changed, baseline)
                for action in (lambda: self.prepare('bad'), lambda: validate_bundle(baseline)):
                    with self.assertRaisesRegex(ValueError, column + ' is inconsistent within guide_id'):
                        action()

    def test_tampered_bundle_excluded_guide_rejected_in_each_used_role(self):
        baseline = self.prepare('baseline')
        design = json.loads((baseline / 'design.json').read_text())
        clean = self.meta.copy()
        for guide in (self.control, self.treated):
            with self.subTest(guide=guide):
                changed = clean.copy()
                changed.loc[changed.guide_id.eq(guide), 'source_excluded_guide'] = True
                self.write_meta(changed, baseline)
                with self.assertRaisesRegex(ValueError, 'includes source_excluded_guide=True'):
                    validate_bundle(baseline)
        # The same control appears in both roles. Test a reference-only and a
        # fit-only design explicitly so neither role can hide behind the other.
        changed = clean.copy()
        changed.loc[changed.guide_id.eq(self.control), 'source_excluded_guide'] = True
        self.write_meta(changed, baseline)
        fit_only = json.loads(json.dumps(design))
        del fit_only['control_groups'][self.control]
        (baseline / 'design.json').write_text(json.dumps(fit_only))
        with self.assertRaisesRegex(ValueError, 'PCA fit includes source_excluded_guide=True'):
            validate_bundle(baseline)
        reference_only = json.loads(json.dumps(design))
        reference_only['pca_fit_indices'] = [i for i in design['pca_fit_indices']
                                           if clean.guide_id.iloc[i] != self.control]
        (baseline / 'design.json').write_text(json.dumps(reference_only))
        with self.assertRaisesRegex(ValueError, 'reference includes source_excluded_guide=True'):
            validate_bundle(baseline)

    def test_tampered_bundle_below_declared_source_minimum_rejected(self):
        baseline = self.prepare('baseline')
        design = json.loads((baseline / 'design.json').read_text())
        design['min_source_cells_per_guide'] = 60
        (baseline / 'design.json').write_text(json.dumps(design))
        changed = self.meta.copy()
        changed.source_guide_n = 100
        changed.loc[changed.guide_id.eq(self.control), 'source_guide_n'] = 52
        self.write_meta(changed, baseline)
        with self.assertRaisesRegex(ValueError, 'below min_source_cells_per_guide=60'):
            validate_bundle(baseline)

    def test_legacy_bundle_default_remains_explicit(self):
        baseline = self.prepare('baseline')
        design = json.loads((baseline / 'design.json').read_text())
        del design['min_source_cells_per_guide']
        (baseline / 'design.json').write_text(json.dumps(design))
        result = validate_bundle(baseline)
        self.assertEqual(result['source_min_cells_checked'], 20)
        self.assertEqual(result['source_min_cells_origin'], 'legacy_default_20')


if __name__ == '__main__':
    unittest.main()
