import gzip
import os
from pathlib import Path
import runpy
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
from scipy import stats

ROOT=Path(__file__).resolve().parents[1]

class CurrentManuscriptTests(unittest.TestCase):
    def test_categorical_severity_residual_correlation(self):
        with tempfile.TemporaryDirectory() as d,patch.dict(os.environ,{'PROJECT_DIR':d}):
            s=runpy.run_path(str(ROOT/'pipelines/05_covid_hca_single_cell_validation/analysis/08_continuous_severity_sensitivity.py'))
            rng=np.random.default_rng(42)
            severity=np.repeat(s['SEVERITY'],12)
            frame=pd.DataFrame({'severity':severity,'candidate_mean_expr':rng.normal(size=60)+np.repeat(np.arange(5)*5,12),'outcome':rng.normal(size=60)+np.repeat(np.arange(5)*5,12)})
            rho,p,df=s['partial_spearman'](frame,'outcome','categorical')
            ranked=frame[['candidate_mean_expr','outcome']].rank()
            residual=ranked-ranked.groupby(frame.severity).transform('mean')
            expected=residual.corr().iloc[0,1]
            self.assertAlmostEqual(rho,expected,12)
            self.assertEqual(df,54)
            self.assertLess(abs(rho),abs(stats.spearmanr(frame.candidate_mean_expr,frame.outcome).statistic))
            self.assertTrue(0<=p<=1)

    def test_s8_distinguishes_full_fit_and_display(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);run=root/'run'
            for x in ('results','prepared','munged'):(run/x).mkdir(parents=True)
            (run/'results/vte_meta_standard_ldsc.log').write_text('After merging with regression SNP LD, 3 SNPs remain.\nMean Chi^2: 1.2736\nLambda GC: 1.1747\nIntercept: 1.0308 (0.0089)\nRatio: 0.1127 (0.0326)\n')
            (run/'prepared/preparation_qc.tsv').write_text('metric\tvalue\ngwas_rows\t10\nmatched_rows\t4\n')
            with gzip.open(run/'munged/vte_meta_standard.sumstats.gz','wt') as h:h.write('SNP\tZ\na\t1\nb\t2\nc\t10\nd\tnan\n')
            with patch.dict(os.environ,{'PROJECT_DIR':d,'STANDARD_LDSC_RUN_DIR':str(run),'TABLE_S8_OUT_DIR':str(root/'out')}):
                s=runpy.run_path(str(ROOT/'pipelines/03_meta_gwas_qc_pipeline/tables/02_make_table_s8_standard_ldsc.py'));s['main']()
            table=pd.read_csv(root/'out/Table_S8_ldsc_calibration_summary.tsv',sep='\t')
            self.assertEqual(table.SNPs.tolist(),[10,4,3,2])
            self.assertAlmostEqual(table.loc[2,'Attenuation ratio SE'],.0326)
            self.assertTrue(pd.isna(table.loc[3,'LDSC intercept']))

if __name__=='__main__':unittest.main()
