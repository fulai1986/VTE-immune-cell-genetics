# Environment

Use Python 3.12 and R. `requirements.txt` pins the scientific Python packages to
the versions present in the checked analysis environment. `python-docx` and
`pyliftover` are required for Word table export and FinnGen coordinate conversion;
their exact original versions were not recorded, so they remain unpinned.
The Python file is a direct-dependency specification, not a complete transitive
lock file. R package versions available during release validation are recorded
in `r-package-versions.tsv`.

Install R packages using `Rscript tools/install_r_packages.R`. Base R `grid` and
`stats` are also used. Fonts use Arial where available; machines without Arial
may render with a fallback font.

## External tools

| Tool | Environment override | Default workspace location |
| --- | --- | --- |
| PLINK 2 | `PLINK2_BIN` | `tools/plink2/plink2` |
| SMR | `SMR_BIN` | `tools/smr/source/SMR/build/RelWithDebInfo/smr` |
| METAL | `METAL_BIN` | `tools/metal/METAL/build/metal/metal` |

Overrides must point to executable files. Bash, `unzip`, `gzip` and `tar` must be
available on PATH. Encrypted archive readers use the system unzip utility.
Record the versions of locally installed executables with each analysis run.
Install or compile these tools for the target operating system; binary files are
not bundled.

The Table S12 and S13 workbook formatters uses `openpyxl`; no Node.js or proprietary
spreadsheet runtime is required.

Standard LDSC uses a separate pinned Python 3.9 environment; see the Result 3
README. Set `LDSC_DIR` and `LDSC_PYTHON` before `meta-qc`.

Figure 4 PDF export uses Quartz on macOS and Cairo on other platforms to preserve
Greek letters and Unicode symbols. The corresponding R graphics device must be
available.
