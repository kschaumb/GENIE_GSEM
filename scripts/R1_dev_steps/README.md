# R1 development steps (archived)

One-off shell scripts used while preparing the JCPP revision (R1), kept for transparency.
They are **superseded by `run_R1_pipeline.sh`** in the project root, which reproduces every
R1 result in a single run. They still use the original working-folder path, so run them only as a record.

| Script | What it did |
|---|---|
| `R1_01_export_cfa_params.sh` | Exported parameter tables from the original (submission) CFA fits |
| `R1_02_run_cfa.sh` | First R1 CFA refit with the Internalizing–Child Physical Development correlation freed |
| `R1_03_run_cfa_symptoms_floor.sh` | Applied the BMI residual floor to the symptom-level CFA |
| `R1_04_tables_and_diagrams.sh` | Rebuilt diagrams and supplemental tables |
| `R1_05_cfa_disorders_bmi_floor.sh` | Restricted the disorder-level residual floor to BMI |
| `R1_06_efa_geomin.sh` | Re-ran EFAs with the preregistered geomin rotation |
| `R1_07_full_rerun_ANX_N.sh` | Full rerun with the corrected anxiety effective N (became `run_R1_pipeline.sh`) |
| `R1_07b_finish_rerun.sh`, `R1_08_render_supp_figures.sh` | Finishing steps, now part of `run_R1_pipeline.sh` |
| `6.CFA_export_params.R1.R` | Helper used by `R1_01` to export parameters from the original CFA fits |
| `R1_compare_rerun.R` | Compared LDSC/CFA results before vs. after the anxiety-N correction |
| `8.Sensitivity.R1.R` | Exploratory BMI ~ AN + OCD model; not reported (fails on a GenomicSEM column name) |

Related outputs are in `gsem/gsem_output/R1_dev_steps/`. The anxiety-N comparison showed unchanged
genetic correlations and standardized estimates changing by ≤ .004.
