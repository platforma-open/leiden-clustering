import pandas as pd
import scanpy as sc
import argparse

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
    df = pd.read_csv(input_csv)

    # Identify which PC value column to use
    if "Principal Component Value" in df.columns:
        pc_value_column = "Principal Component Value"
    elif "Principal Component Value - Harmony corrected" in df.columns:
        pc_value_column = "Principal Component Value - Harmony corrected"
    else:
        raise ValueError("Input CSV must contain either 'Principal Component Value' or 'Principal Component Value - Harmony corrected'.")
    
    # Create a unique identifier: SampleId + CellId
    df["UniqueCellId"] = df["Sample"] + "_" + df["Cell Barcode"]
    
    # Pivot data to have cells as rows, PCs as columns
    pca_matrix = df.pivot(index="UniqueCellId", columns="Principal Component Number", values=pc_value_column)
    
    # Create AnnData object
    adata = sc.AnnData(pca_matrix)
    
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
    
    # Extract cluster assignments
    cluster_assignments = pd.DataFrame({
        "UniqueCellId": adata.obs_names,  # Unique ID as first column
        "Cluster": adata.obs["leiden"]
    })

    # Split SampleId and CellId
    cluster_assignments[["SampleId", "CellId"]] = cluster_assignments["UniqueCellId"].str.split("_", n=1, expand=True)
    
    # Reorder columns to have UniqueCellId first
    cluster_assignments = cluster_assignments[["UniqueCellId", "SampleId", "CellId", "Cluster"]]

    return cluster_assignments

def main():
    parser = argparse.ArgumentParser(description="Run Leiden clustering on PCA embeddings with duplicate CellIds across samples.")
    parser.add_argument("--input_csv", type=str, required=True, help="Path to the PCA embeddings CSV file.")
    parser.add_argument("--output_csv", type=str, required=True, help="Path to save cluster assignments.")
    parser.add_argument("--linker_csv", type=str, required=True, help="Path to save linker data for pFrame construction.")
    parser.add_argument("--n_neighbors", type=int, default=15, help="Number of neighbors for the graph (default: 15).")
    parser.add_argument("--leiden_resolution", type=float, default=1.0, help="Resolution for Leiden clustering (default: 1.0).")

    args = parser.parse_args()

    # Process PCA embeddings
    adata = process_pca_embeddings(args.input_csv, args.n_neighbors)

    # Perform clustering
    cluster_assignments = perform_clustering(adata, args.leiden_resolution)

    # Create linker data: SampleId,CellId,Cluster,Link
    linker_data = cluster_assignments[["SampleId", "CellId", "Cluster"]].copy()
    linker_data["Link"] = 1

    # Save outputs
    cluster_assignments.to_csv(args.output_csv, index=False)
    linker_data.to_csv(args.linker_csv, index=False)

if __name__ == "__main__":
    main()
