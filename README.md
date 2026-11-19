# GENIE Genomic SEM

Code and output for **"Mapping Shared and Unique Genetic Liability to Eating and Internalizing Psychopathology"**
(Schaumberg et al., under review at *Journal of Child Psychology and Psychiatry*). Preregistration: <https://osf.io/jpzrk/overview>.

The study uses Genomic Structural Equation Modeling ([GenomicSEM](https://github.com/GenomicSEM/GenomicSEM), v0.0.5) on publicly
available GWAS summary statistics to model the shared genetic structure of eating disorder, internalizing, anthropometric,
and pubertal development traits at the disorder and symptom level (Aims 1–2), and conditional genetic associations with
eating disorder phenotypes after accounting for BMI (Aim 3).

## Quick start (revised analyses, R1)

```bash
bash run_R1_pipeline.sh
```

This runs every step behind the revised manuscript, from GWAS settings through the supplemental tables and figures,
and writes logs to `gsem/gsem_output/log_files/`. Rendering the supplemental figures requires [Quarto](https://quarto.org).

## Repository layout

| Path | Contents |
|---|---|
| `run_R1_pipeline.sh` | Single entry point for all revised (R1) analyses |
| `scripts/` | Analysis scripts used in the revision (`*.R1.*`), numbered in run order, plus `0.Packages.R` |
| `scripts/cleaning_scripts/` | Harmonization of each raw GWAS into a common format before munging |
| `scripts/R1_dev_steps/` | Archived one-off steps from the revision (superseded by `run_R1_pipeline.sh`) |
| `data/gwas_sumstats/GWAS_sumstats_settings*.xlsx` | Metadata for every GWAS (source, N, prevalences, file names). `_R1` applies the documented corrections in `scripts/0.Settings_corrections.R1.R` |
| `gsem/mvLD_R1/` | LDSC genetic covariance (S) and sampling covariance (V) matrices used by all R1 models |
| `gsem/gsem_output/` | Model output: `RGs_R1/` (Figure 1), `*_efa_R1/`, `*_cfa_R1/` (Figures 2–3), `aim3_conditional_R1/` (Figure 4), `Supplemental_Tables_R1.xlsx`, `Supplemental_Figures.R1.pdf` |

## R1 pipeline steps

| Step | Script | Output |
|---|---|---|
| Settings corrections | `scripts/0.Settings_corrections.R1.R` | `GWAS_sumstats_settings_R1.xlsx` |
| Settings | `scripts/0.Settings.R1.R` | `gsem/gsem_settings_R1/` |
| Munging (anxiety, corrected N) | `scripts/1.Munge.R1.R` | `data/gwas_sumstats/munged/hm3_munged_R1/` |
| LDSC | `scripts/2.LDSC.R1.R` | `gsem/mvLD_R1/` |
| Genetic correlation matrices | `scripts/4.RGs.R1.R` | `gsem/gsem_output/RGs_R1/` |
| Exploratory factor analysis (ML, geomin) | `scripts/5.EFA.R1.r` | `gsem/gsem_output/{disorder,symptoms}_efa_R1/` |
| Confirmatory factor analysis (DWLS) | `scripts/6.CFA_disorders.R1.R`, `scripts/6.CFA_symtpoms.R1.R` | `gsem/gsem_output/{disorder,symptoms}_cfa_R1/` |
| Aim 3 trait trios and genomic regression | `scripts/8.Conditional_Overlap.R1.r`, `scripts/9.GeneticMR.R1.r` | `gsem/gsem_output/aim3_conditional_R1/` |
| Path diagrams | `scripts/7.CFA_Diagrams.R1.r` | `gsem/gsem_output/*_cfa_R1/diagrams/` |
| Supplemental tables S1–S9 | `gsem/gsem_output/Supplemental_Tables.R1.R` | `Supplemental_Tables_R1.xlsx` |
| Supplemental figures S1–S4 | `gsem/gsem_output/Supplemental_Figures.R1.qmd` | `Supplemental_Figures.R1.pdf` |

## Data

GWAS summary statistics are not redistributed here (~3 GB). Sources, sample sizes, and download references for each trait
are listed in Table S1 and in `data/gwas_sumstats/GWAS_sumstats_settings_R1.xlsx`. To reproduce from scratch, place:

- raw summary statistics in `data/gwas_sumstats/<trait>/` and run the matching script in `scripts/cleaning_scripts/`;
- the munged files (`<code>.sumstats.gz`, HapMap 3 SNPs, INFO > 0.8, MAF > 0.01) in `data/gwas_sumstats/munged/hm3_munged_info_0_8_maf_0_1/` (produced by the original munging script, `scripts/old/1.Munge.R`);
- European LD scores (`eur_w_ld_chr`) in `data/ld_scores/` and `w_hm3.snplist` in `data/variant_lists/`
  (both available from the [GenomicSEM resources](https://github.com/GenomicSEM/GenomicSEM/wiki)).

The LDSC matrices in `gsem/mvLD_R1/` are included, so all factor analyses and Aim 3 models can be re-run without the raw data
(start `run_R1_pipeline.sh` from step 05 by commenting out steps 01–04).

## Software

R with GenomicSEM (v0.0.5), lavaan, psych, GPArotation, dplyr, readr, data.table, openxlsx, ggplot2, gt
(see `scripts/0.Packages.R`), and Quarto for the supplemental figures.

## Key conventions

- Factor analyses use the genetic correlation matrix; EFA uses maximum likelihood with oblique geomin rotation (preregistered).
- In final CFAs, the residual variance of adult BMI is fixed at .01 to resolve a Heywood case (preregistered procedure).
- The Child Physical Development factor is reported so that higher values reflect higher childhood BMI and earlier menarche.
- Effective N is used for all binary traits (sample prevalence 0.5 for liability-scale conversion).
