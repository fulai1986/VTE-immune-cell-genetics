# Shared reference preparation

Set `PROJECT_DIR` to the analysis workspace. The root workflow runner places
these steps at the correct points between expression mapping and downstream
association analysis.

1. After the genotype sample map is available, run
   `16_make_onek1k_ld_reference.py --allele-specific-ids` to generate the OneK1K
   PLINK reference.
2. After cis-eQTL mapping, run `18_prepare_smr_qfile_besd.py --make-besd` to export
   allele-matched qfiles and BESD files. These qfiles also support eGene correction.
3. Run `09_liftover_finngen_r13_to_grch37.py` and then
   `10_harmonize_finngen_r13_to_onek1k_ld.py` before FinnGen SMR.

See `docs/DATA_INPUTS.md` and `docs/ENVIRONMENT.md` at the repository root for data,
credential and executable configuration. The preparation code uses the same
variant identities and allele definitions as the manuscript analysis.
