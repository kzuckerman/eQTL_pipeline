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
        cell_type_identifier=cluster,   # use the cell-type/cluster column passed from the WDL
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
    dc.pp.filter_samples(pdata, min_cells=2, min_counts=10)
    
    # Normalize total counts and log-transform
    sc.pp.normalize_total(pdata, target_sum=1e6)
    sc.pp.log1p(pdata)
    
    # Filter genes with mean log-transformed expression >= 3
    gene_filter = pdata.X.mean(axis=0) >= 3
    pdata = pdata[:, gene_filter].copy()  # .copy() so scaling below persists (slice is a view)

    # Scale the data to a maximum value of 10.
    # NOTE: scale in place on `pdata` itself. Passing pdata.copy() scaled a
    # throwaway object and left `pdata` unscaled.
    sc.pp.scale(pdata, max_value=10)
    
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
    
    # Composition matrix = cell-type counts ONLY.
    # data_scanpy_1.obs carries per-sample metadata (condition, batch, etc.).
    # Previously it was concatenated in and every column was later swept up as an
    # eQTL covariate downstream. We deliberately drop .obs and keep only the
    # cell-type columns here.
    composition = pd.DataFrame(
        data_scanpy_1.X,
        index=data_scanpy_1.obs.index,
        columns=data_scanpy_1.var.index,
    )
    composition.index.name = None  # keep an unnamed index so the CSV header reads "X" in R

    if "leiden" in cluster:
        composition.columns = ['cluster' + str(col) for col in data_scanpy_1.var.index]
    else:
        composition.columns = list(data_scanpy_1.var.index)

    composition_aligned = composition.loc[common_samples]
    composition_aligned.to_csv(output_file2)
    print(f"Saved {output_file2}")
print("Processing completed.")