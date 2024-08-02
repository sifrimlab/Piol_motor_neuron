# 10X Visium Data Analysis

This folder contains a notebook that executes and reproduces the figures for the motor and non-motor neuron signature enrichment analysis in 10X Visium spinal cord data.

## Running the Notebook

To run the notebook, follow these steps:

1. **Install and activate the conda environment:**

   ```bash
   conda env create -f OtF.yml
   conda activate OtF
   ```

2. **Download the processed 10X Visium `.h5ad` file from GEO:**

   - Dataset: [GSE269377_adata_10X_Visium_mouse_SC_all_Fus_multisample.h5ad](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269377)
