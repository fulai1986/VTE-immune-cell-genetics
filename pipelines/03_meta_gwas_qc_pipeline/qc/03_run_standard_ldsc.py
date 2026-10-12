#!/usr/bin/env python3
"""Prepare alleles from the primary GWAS and run standard munging/LDSC.

LDSC_DIR and LDSC_PYTHON point to the documented external LDSC installation.
Existing run directories are refused to protect completed analyses.
"""
import gzip
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess

PROJECT=Path(os.environ['PROJECT_DIR']).expanduser().resolve()
RUN=Path(os.environ.get('STANDARD_LDSC_RUN_DIR',PROJECT/'gwas/ldsc_standard')).resolve()

def main():
    ldsc=Path(os.environ['LDSC_DIR']).expanduser().resolve()
    python=os.environ['LDSC_PYTHON']
    for script in ('munge_sumstats.py','ldsc.py'):
        if not (ldsc/script).is_file():raise FileNotFoundError(ldsc/script)
    if RUN.exists() and any(RUN.iterdir()):
        raise FileExistsError(f'Completed/partial run exists: {RUN}. Reuse it with meta-tables or select a new STANDARD_LDSC_RUN_DIR.')
    ref=PROJECT/'tools/ldsc_reference/eur_w_ld_chr'
    spec=importlib.util.spec_from_file_location('prepare',Path(__file__).with_name('01_prepare_ldsc_sumstats_ukb_approx_logor.py'))
    prepare=importlib.util.module_from_spec(spec);spec.loader.exec_module(prepare)
    prepare.OUT_DIR=RUN/'prepared'
    prepare.SUMSTATS=prepare.OUT_DIR/'vte_meta_hm3_from_raw.tsv.gz'
    prepare.MATCHED=prepare.OUT_DIR/'vte_meta_hm3_from_raw.ldscore_matched.tsv.gz'
    prepare.QC=prepare.OUT_DIR/'preparation_qc.tsv'
    prepare.main()
    plain=prepare.OUT_DIR/'vte_meta_hm3_from_raw.tsv'
    with gzip.open(prepare.SUMSTATS,'rb') as source,plain.open('wb') as target:shutil.copyfileobj(source,target)
    for folder in ('munged','results'):(RUN/folder).mkdir(parents=True,exist_ok=True)
    munged=RUN/'munged/vte_meta_standard'
    subprocess.run([python,str(ldsc/'munge_sumstats.py'),'--sumstats',str(plain),'--merge-alleles',str(ref/'w_hm3.snplist'),'--signed-sumstats','Z,0','--out',str(munged)],check=True,cwd=RUN)
    subprocess.run([python,str(ldsc/'ldsc.py'),'--h2',str(munged)+'.sumstats.gz','--ref-ld-chr',str(ref)+'/', '--w-ld-chr',str(ref)+'/', '--out',str(RUN/'results/vte_meta_standard_ldsc')],check=True,cwd=RUN)

if __name__=='__main__':main()
