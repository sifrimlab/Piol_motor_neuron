# Nanostring_GeoMx

This section contains code to reproduce the count matrix, results and downstream GSEA for Nanostring GeoMx, as well as select pre-rendered reports displaying code and results. The raw data used for this analysis can be found on GEO [GSE269707](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269707)

## Running the Notebook

To run the notebook, follow these steps:

1. **Install and activate the conda environment:**

   ```bash
   conda env create -f R_env.yml
   conda activate R_env
   ```
2. Run Scripts 

The scripts for this workflow should be run in the following order:
1. `NanoString_Exploratory_Analysis_Final.R` convert `.dcc` and `.pcc` files into long count matrix (optional)
2. `NS_make_raw_counts.R` to re-generate metadata and processed count matrix files. These files are also found in `./data` (optional)

Note: you may need to change the input file path, based on the file architecture of your local machine:

Then run the DESeq2 DEG analyses workbooks:
3. `NS_sciatic_nerve_Chatpos_ctrl_v_Chatneg_ctrl_condition_model.Rmd`
4. `NS_sciatic_nerve_Chatpos_ctrl_v_Chatneg_ctrl_segment_model.Rmd`
5. `NS_sciatic_nerve_Chatpos_ctrl_v_Chatpos_FUS.Rmd`
6. `NS_spinal_cord_Chatpos_ctrl_v_Chatneg_ctrl.Rmd`

To perform GSEA on the DEG lists derived from the DESEq2 analyses listed above, run:
7. `NS_FGSEA_on_all_DE_res.R`

Assuming that you have already installed the R environment dependencies for these scripts,
they can be run in the command line as follows, or ran in RStudio:
- Rscripts, e.g.: `Rscript NS_make_raw_counts.R`
- Rmarkdown workbooks, e.g.: `Rscript -e "rmarkdown::render('your_file.Rmd')"`
