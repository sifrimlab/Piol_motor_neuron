import numpy as np
import pandas as pd
from skimage import io
from pathlib import Path
import matplotlib.pyplot as plt
from pyxelperfect.measure import measureLabeledImage

def createGeneMatrix(labeled_image: np.array, decoded_df: pd.DataFrame, rowname: str = "row", colname: str ="col") -> pd.DataFrame:

    # Wrapper around skimage.measure.regionprops
    measured_df = measureLabeledImage(labeled_image)

    def assignGenesToSpots(labeled_image, decoded_df, cell_properties):
        gene_dict = {} # keys = label integers, values = dict{gene: count}
        n_labels = np.unique(labeled_image)
        for label in n_labels:
            gene_dict[label] = {}
            
        for row in decoded_df.itertuples():
            label = labeled_image[int(row.row), int(row.col)]
            if label != 0:
                gene_dict[label][row.gene] = gene_dict[label].get(row.gene, 0) + 1 

        # then rename the label to me a bit more clear
        for label in n_labels:
            if label == 0:
                continue
            gene_dict[f"{cell_properties[cell_properties['image_label'] == label].iloc[0]['cell_label']}"] = gene_dict.pop(label)
        return gene_dict


    gene_dict = assignGenesToSpots(labeled_image, decoded_df, measured_df)

    # Just to make sure that the background isn't accidentally in there
    try:
        gene_dict.pop(0)
    except KeyError:
        print("No label 0 found.")

    sorted_gene_keys = sorted(list(gene_dict.keys()), key=lambda x: int(x.split('_')[-1]))

    gene_dict_list = [gene_dict[i] for i in sorted(list(gene_dict.keys()), key=lambda x: int(x.split('_')[-1]))]

    keys = set().union(*gene_dict_list)
    final = {k: [d.get(k, 0) for d in gene_dict_list] for k in keys}

    new_df = pd.DataFrame(final)
    new_df["cell_label"] = sorted_gene_keys
    new_df.set_index("cell_label", inplace=True)

    return new_df
