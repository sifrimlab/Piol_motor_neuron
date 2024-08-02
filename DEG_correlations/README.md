# DEG_correlations

This section contains R scripts to perform correlations between publicly available motor neuron datasets from GEO and data from the Piol, et al. manuscript. The raw data for these datasets can be downloaded from GEO from the following links:

1. [Gautier, et al. 2023 - GSE228778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228778)
2. [Yadav, et al. 2023 - GSE190442](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190442)
3. [Alkaslasi, et al. 2021 - GSE167597](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167597)
4. [Blum, et al. 2021 - GSE161621](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161621)

Note:
- author annotated versions of the datasets listed above that were used in the Piol, et al. study can also be downloaded from [Spinal Cord Atlas](http://spinalcordatlas.org/), or by contacting the respective authors of these studies.
- you may need to change the input file path, based on the file architecture of your local machine
- Seurat objects of these datasets can be made available upon request

## Running the Notebook

To run the notebook, follow these steps:

1. **Install and activate the conda environment:**

   ```bash
   conda env create -f R_env.yml
   conda activate R_env
   ```
2. Run Scripts 

To reproduce the data from this section, download the raw data outlined above and then run the scripts found in the `./DEG_correlations` folder, in the following order:
1. `GSE161621_Blum_et_al_DE.R` performs Seurat DE between Chat+ vs Chat - nuclei in the data of Blum, et al. 2021
2. `GSE167597_Alkaslasi_et_al_DE.R` performs Seurat DE between Chat+ vs Chat - nuclei in the data of Alkaslasi, et al. 2021
3. `GSE190442_Yadav_et_al_DE.R` performs Seurat DE between Chat+ vs Chat - nuclei in the data of Yadav, et al. 2021
4. `GSE228778_Gautier_et_al_DE.R` performs Seurat DE between Chat+ vs Chat - nuclei in the data of Gautier, et al. 2023
5. `NS_vs_GEO_DEG_list_corr.Rmd` this script takes the DEG lists from the four scripts listed above and correlates the LFC values between matched DEGs derived from the Nanostring comparisons

Assuming that you have already installed the R environment dependencies for these scripts,
they can be run in the command line as follows, or ran in RStudio:
- Rscripts, e.g.: `GSE161621_Blum_et_al_DE.R`
- Rmarkdown workbooks, e.g.: `Rscript -e "rmarkdown::render('NS_vs_GEO_DEG_list_corr.Rmd')"`
