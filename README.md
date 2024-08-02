# Piol_motor_neuron
This repository contains a collection of notebooks, scripts and config files used to reproduce the analyses, results and figures in the *Piol, et al.* manuscript. 

## Contents

- [Config](#Config)
- [10X_Visium](#10X_Visium)
- [Nanostring](Nanostring_GeoMx)
- [DEG_correlations](#DEG_correlations)
- [Data](#Data)
- [Plots](#Plots)
- [LICENSE](./LICENSE)
- [Contact](#Contact)

# Config

To set up the necessary environment and install all required packages to 
reproduce the analyses found in this repository, we recommend using Conda, and
installing the packages with their specified versions from the `.yaml` files
found in `./config/yamls/`, outlined in the steps below: 

1. Ensure you have Conda installed. You can download it from the [official Conda website](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

2. Clone the repository to your local machine:

    ```sh
    git clone https://github.com/sifrimlab/Piol_motor_neuron/
    cd Piol_motor_neuron
    ```

3. Install the necessary packages using the provided Conda environment file,
(example R environment shown here):

    ```sh
    conda env install -f ./configs/yamls/R_env.yaml
    ```

4. Activate the newly created environment:

    ```sh
    conda activate R_env
    ```

This will set up the environment with all the dependencies required for the project. The time to install the Conda dependencies can vary, depending on your machine. If it takes too long, we suggest using [Mamba](https://anaconda.org/conda-forge/mamba) to install these dependencies.

## Yamls

There are three `.yaml` files detailing the packages for different three
environments used for the Piol, et al. analysis. We recommend creating three
separate Conda environments (in the manner as outlined above) using the following
`.yaml` files:
1. [10X Visium software Python environment](./configs/yamls/OtF.yaml) `.yaml` file describing dependencies to reproduce the 10X Visium portion of this study
2. [Nanostring and Seurat software R environment](./configs/yamls/R_env.yaml) `.yaml` file describing dependencies to reproduce the Nanostring preprocessing, DESeq2 DEG analysis, GSEA, GEO data Seurat DEG analysis and correlations portions of this study
3. [Image processing Python environment](./configs/yamls/img_processing.yaml) `.yaml` file describing dependencies to reproduce the imaging processing protion of this study

Note: if you encounter issues installing all of the R dependencies via Conda for
the R environment, the necessary R packages can also be installed in an R 
environment using the `install_packages.R` script found in `./config/R_env/`.

# 10X_Visium

This section contains code to reproduce the results for 10X Visium untargeted spatial transcriptomics from the the Piol, et al. manuscript. The raw data used for this analysis can be found on GEO [GSE269377](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269377).

All code to reproduce this portion of the study is found in the `10XVisium_SC_marker_gene_enrich.ipynb` Jupyter notebook

# Nanostring_GeoMx

This section contains code to reproduce the count matrix, results and downstream GSEA for Nanostring GeoMx, as well as select pre-rendered reports displaying code and results. The raw data used for this analysis can be found on GEO [GSE269707](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269707)

The scripts for this workflow should be run in the following order:
1. `NanoString_Exploratory_Analysis_Final.R` convert `.dcc` and `.pcc` files into long count matrix (optional)
2. `NS_make_raw_counts.R` to re-generate metadata and processed count matrix files. These files are also found in `./data` (optional)

Note: you may need to change the input file path, based on the file architecture of your local machine:

Then run the DESeq2 DEG analyses workbooks:
1. `NS_sciatic_nerve_Chatpos_ctrl_v_Chatneg_ctrl_condition_model.Rmd`
2. `NS_sciatic_nerve_Chatpos_ctrl_v_Chatneg_ctrl_segment_model.Rmd`
3. `NS_sciatic_nerve_Chatpos_ctrl_v_Chatpos_FUS.Rmd`
4. `NS_spinal_cord_Chatpos_ctrl_v_Chatneg_ctrl.Rmd`

To perform GSEA on the DEG lists derived from the DESEq2 analyses listed above, run:
5. `NS_FGSEA_on_all_DE_res.R`

Assuming that you have already installed the R environment dependencies for these scripts,
they can be run in the command line as follows, or ran in RStudio:
- Rscripts, e.g.: `Rscript NS_make_raw_counts.R`
- Rmarkdown workbooks, e.g.: `Rscript -e "rmarkdown::render('NS_spinal_cord_Chatpos_ctrl_v_Chatneg_ctrl.Rmd')"`

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

To reproduce the data from this section, download the raw data outlined above and then run the scripts found in the `./DEG_correlations` folder, in the following order:
1. `GSE161621_Blum_et_al_DE.R` performs Seurat DE between Chat+ vs Chat - nuclei in the data of Blum, et al. 2021
2. `GSE167597_Alkaslasi_et_al_DE.R` performs Seurat DE between Chat+ vs Chat - nuclei in the data of Alkaslasi, et al. 2021
3. `GSE190442_Yadav_et_al_DE.R` performs Seurat DE between Chat+ vs Chat - nuclei in the data of Yadav, et al. 2021
4. `GSE228778_Gautier_et_al_DE.R` performs Seurat DE between Chat+ vs Chat - nuclei in the data of Gautier, et al. 2023
5. `NS_vs_GEO_DEG_list_corr.Rmd` this script takes the DEG lists from the four scripts listed above and correlates the LFC values between matched DEGs derived from the Nanostring comparisons

# Data

This section contains raw data files for Nanostring, as well processed DEG gene lists derived from Seurat analysis of the publicly available motor neuron datasets from GEO 

# Plots

This section contains select visualizations that appear in the Piol, et al. manuscript

# Contact

For any questions or additional information, please contact:

- **Name:** Theo Killian
- **Email:** [theo.killian@kuleuven.be](mailto:theo.killian@kuleuven.be)
- **Affiliation:** Da Cruz/Sifrim Lab Bioinformatics
