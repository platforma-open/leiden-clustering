import polars as pl
import pandas as pd
import scanpy as sc
import argparse
import numpy as np


def process_pca_embeddings(input_csv, n_neighbors):
    """
    Process PCA embeddings from CSV and construct AnnData object.
    
    Parameters:
        input_csv (str): Path to PCA embeddings CSV file.
        n_neighbors (int): Number of neighbors for graph construction.

    Returns:
        adata (AnnData): Annotated data object with PCA embeddings.
    """
    # Load PCA embeddings
    df = pl.read_csv(input_csv)

    # Identify which PC value column to use
    if "Principal Component Value" in df.columns:
        pc_value_column = "Principal Component Value"
    elif "Principal Component Value - Harmony corrected" in df.columns:
        pc_value_column = "Principal Component Value - Harmony corrected"
    else:
        raise ValueError("Input CSV must contain either 'Principal Component Value' or 'Principal Component Value - Harmony corrected'.")

    # Validate and normalize column headers to support both legacy and new names
    base_required = {"Sample", "Principal Component Number", pc_value_column}
    cell_headers = {"Cell Barcode", "Cell ID"}
    missing_base = base_required - set(df.columns)
    has_cell_header = any(h in df.columns for h in cell_headers)
    if missing_base or not has_cell_header:
        expected_desc = f"{sorted(base_required)} and one of {sorted(cell_headers)}"
        raise KeyError(f"PCA CSV must contain columns: {expected_desc}. Found: {list(df.columns)}")

    # Normalize to legacy internal name 'Cell Barcode'
    if "Cell ID" in df.columns and "Cell Barcode" not in df.columns:
        df = df.rename({"Cell ID": "Cell Barcode"})

    # Create a unique identifier: SampleId + CellId
    df = df.with_columns(
        (pl.col("Sample").cast(pl.Utf8) + "_" + pl.col("Cell Barcode").cast(pl.Utf8)).alias("UniqueCellId")
    )

    # Sort by cell and PC number to ensure correct order for reshaping
    df_sorted = df.sort("UniqueCellId", "Principal Component Number")

    # Get dimensions for reshaping
    n_cells = df_sorted.get_column("UniqueCellId").n_unique()
    n_pcs = df_sorted.get_column("Principal Component Number").n_unique()

    # Reshape values into a dense matrix
    pca_matrix_values = df_sorted.get_column(pc_value_column).to_numpy().reshape((n_cells, n_pcs))

    # Get cell and PC names for AnnData object creation
    cell_names = df_sorted.get_column("UniqueCellId").unique().to_list()
    pc_names = df_sorted.get_column("Principal Component Number").unique().to_list()

    # Create AnnData object
    adata = sc.AnnData(
        X=pca_matrix_values,
        obs=pd.DataFrame(index=cell_names),
        var=pd.DataFrame(index=pc_names)
    )

    # Compute the neighborhood graph
    sc.pp.neighbors(adata, use_rep="X", n_neighbors=n_neighbors)

    return adata


def perform_clustering(adata, leiden_resolution):
    """
    Perform Leiden clustering on an AnnData object.
    
    Parameters:
        adata (AnnData): Annotated data object.
        leiden_resolution (float): Resolution parameter for Leiden clustering.

    Returns:
        cluster_assignments (DataFrame): Dataframe with cell clustering assignments.
    """
    # Perform Leiden clustering
    sc.tl.leiden(adata, resolution=leiden_resolution, flavor='igraph', n_iterations=2, directed=False)
    
    # Create a Polars Series for cluster labels and prepend "CL-" efficiently
    cluster_series = "CL-" + pl.Series(adata.obs["leiden"]).cast(pl.Utf8)

    # Extract cluster assignments into a Polars DataFrame
    cluster_assignments = pl.DataFrame({
        "UniqueCellId": adata.obs_names.to_list(),
        "Cluster": cluster_series
    })

    # Split UniqueCellId into SampleId and CellId
    cluster_assignments = cluster_assignments.with_columns([
        pl.col("UniqueCellId").str.split_exact("_", 1).struct.field("field_0").alias("SampleId"),
        pl.col("UniqueCellId").str.split_exact("_", 1).struct.field("field_1").alias("CellId"),
    ])
    
    # Reorder columns
    cluster_assignments = cluster_assignments.select(["UniqueCellId", "SampleId", "CellId", "Cluster"])

    return cluster_assignments


def main():
    parser = argparse.ArgumentParser(description="Run Leiden clustering on PCA embeddings with duplicate CellIds across samples.")
    parser.add_argument("--input_csv", type=str, required=True, help="Path to the PCA embeddings CSV file.")
    parser.add_argument("--output_csv", type=str, required=True, help="Path to save cluster assignments.")
    parser.add_argument("--linker_csv", type=str, default="leiden_linker.csv", help="Path to save linker data for pFrame construction.")
    parser.add_argument("--n_neighbors", type=int, default=15, help="Number of neighbors for the graph (default: 15).")
    parser.add_argument("--leiden_resolution", type=float, default=1.0, help="Resolution for Leiden clustering (default: 1.0).")

    args = parser.parse_args()

    # Process PCA embeddings
    adata = process_pca_embeddings(args.input_csv, args.n_neighbors)

    # Perform clustering
    cluster_assignments = perform_clustering(adata, args.leiden_resolution)

    # Create linker data: SampleId,CellId,Cluster,Link
    linker_data = cluster_assignments.select(["SampleId", "CellId", "Cluster"])
    linker_data = linker_data.with_columns(pl.lit(1).alias("Link"))

    # Save outputs
    cluster_assignments.write_csv(args.output_csv)
    linker_data.write_csv(args.linker_csv)


if __name__ == "__main__":
    main()
