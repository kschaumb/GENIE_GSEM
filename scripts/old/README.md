# Original (submission) scripts — archived

These are the analysis scripts used for the original JCPP submission. Each has been superseded by a
`*.R1.*` version in `scripts/` (run via `run_R1_pipeline.sh`), except `3.GSEM_sumstats.R` and `clean_labels.R`,
which are not used in the revised analyses. They are kept so the original results can be traced; they use
the original absolute paths and settings (`gsem/old/gsem_settings/`, `gsem/mvLD/`) and are not maintained.

`1.Munge.R` munged all GWAS for the original analyses; the munged files it produced are reused by the R1 pipeline
for every trait except anxiety (re-munged with corrected N by `scripts/1.Munge.R1.R`).

`gsem_output/` holds the original supplemental table and figure scripts (superseded by
`gsem/gsem_output/Supplemental_Tables.R1.R` and `Supplemental_Figures.R1.qmd`).
