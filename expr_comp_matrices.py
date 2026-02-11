import os
import sys
import scanpy as sc
import pandas as pd
import decoupler as dc
from sccoda.util import cell_composition_data as dat

input_files = [f for f in sys.argv[2:] if f[-5:] == ".h5ad"]
cluster_names = [f for f in sys.argv[2:] if f[-5:] != ".h5ad"]

workdir = sys.argv[1]

# Process each file
for file, cluster in zip(input_files, cluster_names):
    print(f"Processing {file}...")
    
    # Load the h5ad file
    adata = sc.read_h5ad(file)
    adata.obs_names_make_unique()

    data_scanpy_1 = dat.from_scanpy(
        adata,
        cell_type_identifier="cell_class",
        sample_identifier="Sample"
    )

    # Pseudobulk aggregation
    pdata = dc.pp.pseudobulk(
        adata,
        sample_col="Sample",  # Use "Sample" column for aggregation
        groups_col="Sample",
        skip_checks=True,
        mode="sum",
        empty=True
    )
    del adata
    dc.pp.filter_samples(pdata, min_cells=5, min_counts=100)
    
    # Normalize total counts and log-transform
    sc.pp.normalize_total(pdata, target_sum=1e6)
    sc.pp.log1p(pdata)
    
    # Filter genes with mean log-transformed expression >= 3
    gene_filter = pdata.X.mean(axis=0) >= 3
    pdata = pdata[:, gene_filter]
    
    # Scale the data to a maximum value of 10
    sc.pp.scale(pdata.copy(), max_value=10)
    
    # Save the processed data as a CSV file
    output_file1 = os.path.join(f"{workdir}/processed_matrices/{os.path.splitext(os.path.basename(file))[0]}_expression_matrix_ds.csv")
    output_file2 = os.path.join(f"{workdir}/processed_matrices/{os.path.splitext(os.path.basename(file))[0]}_composition_matrix_ds.csv")

    # Filter both datasets to only include common samples
    data = pd.DataFrame(pdata.X, index=pdata.obs_names, columns=pdata.var_names)
    common_samples = data.index.intersection(data_scanpy_1.obs.index)
    del pdata

    data_aligned = data.loc[common_samples]
    data_aligned.to_csv(output_file1)
    print(f"Saved {output_file1}")
    del data, data_aligned
    
    combined_df_corrected =pd.concat([data_scanpy_1.obs, pd.DataFrame(data_scanpy_1.X, index=data_scanpy_1.obs.index, columns=data_scanpy_1.var.index)],axis=1)
    if "leiden" in cluster:
        new_column_names = ['cluster' + str(col) for col in data_scanpy_1.var.index]
    else:
        new_column_names = list(data_scanpy_1.var.index)
    # Update the relevant columns in combined_df_corrected with the new names
    combined_df_corrected.columns = list(data_scanpy_1.obs.columns) + new_column_names

    combined_df_corrected_aligned = combined_df_corrected.loc[common_samples]
    combined_df_corrected_aligned.to_csv(output_file2)
    print(f"Saved {output_file2}")
print("Processing completed.")
