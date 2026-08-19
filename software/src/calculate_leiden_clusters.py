import argparse

import polars as pl
import scanpy as sc

# The exported frame labels its value column after the upstream p-column, so a
# batch-corrected run arrives under a different header than a plain one.
PC_VALUE_COLUMNS = (
    "Principal Component Value",
    "Principal Component Value - Harmony corrected",
)
CELL_HEADERS = ("Cell Barcode", "Cell ID")
SAMPLE_HEADER = "Sample"
PC_NUMBER_HEADER = "Principal Component Number"


"""Ordering note.

Leiden is order-dependent: its local-moving phase walks nodes by index, so two
runs over the *same* neighbour graph give different partitions if the rows were
fed in a different order. The pre-Parquet implementation pivoted with pandas,
which sorts both the index and the columns of its result; polars preserves
first-appearance order instead. Both axes are therefore pinned here to the order
pandas produced, so this version reproduces the partitions the block emitted
before the migration.

Row order is what actually matters. Column order is provably irrelevant to the
graph (permuting features leaves every pairwise distance unchanged), but it is
pinned too so the whole pipeline is reproducible by inspection rather than by
argument.
"""


def load_embeddings(input_parquet):
    """Read the long-format PCA embeddings and pivot them to a cell x PC matrix.

    Returns the wide frame (sample id, cell id, one column per PC), the cell-id
    header that was actually present, and the ordered PC column names.
    """
    scan = pl.scan_parquet(input_parquet)
    schema = scan.collect_schema()
    column_names = set(schema.names())

    pc_value_column = next((c for c in PC_VALUE_COLUMNS if c in column_names), None)
    if pc_value_column is None:
        raise ValueError(
            f"Embeddings Parquet must contain one of {sorted(PC_VALUE_COLUMNS)}. "
            f"Found: {sorted(column_names)}"
        )

    # Support both the legacy 'Cell Barcode' header and the current 'Cell ID'.
    cell_column = next((h for h in CELL_HEADERS if h in column_names), None)
    base_required = {SAMPLE_HEADER, PC_NUMBER_HEADER}
    missing_base = base_required - column_names
    if missing_base or cell_column is None:
        expected_desc = f"{sorted(base_required)} and one of {sorted(CELL_HEADERS)}"
        raise KeyError(
            f"Embeddings Parquet must contain columns: {expected_desc}. "
            f"Found: {sorted(column_names)}"
        )

    # Sample and cell ids repeat once per principal component, so holding them as
    # categoricals keeps the long frame small. Parquet carries its own dtypes, so
    # the cast rides in the scan plan instead of a read-time schema override.
    long_df = scan.with_columns(
        [pl.col(c).cast(pl.Categorical) for c in (SAMPLE_HEADER, cell_column)]
    ).collect()

    # Pivot on the (sample, cell) pair rather than on a concatenated key, so the
    # two ids stay in their own columns and never have to be split apart again.
    wide = long_df.pivot(
        on=PC_NUMBER_HEADER,
        index=[SAMPLE_HEADER, cell_column],
        values=pc_value_column,
    )

    # Lexicographic, matching pandas — see the ordering note at the top.
    pc_columns = sorted(c for c in wide.columns if c not in (SAMPLE_HEADER, cell_column))

    # Restore the dtypes the input carried. Downstream pfconv reads these back
    # against the axis specs of the source p-column, so the Parquet we emit has
    # to hold the same physical types the Parquet we read did.
    wide = wide.with_columns(
        pl.col(SAMPLE_HEADER).cast(schema[SAMPLE_HEADER]),
        pl.col(cell_column).cast(schema[cell_column]),
    )

    # The row key doubles as the sort key and as the emitted UniqueCellId. Sorting
    # on the joined string rather than on (sample, cell) reproduces pandas' order
    # exactly, including where one sample name is a prefix of another.
    wide = wide.with_columns(
        pl.concat_str(
            [
                pl.col(SAMPLE_HEADER).cast(pl.String),
                pl.lit("_"),
                pl.col(cell_column).cast(pl.String),
            ]
        ).alias("UniqueCellId")
    ).sort("UniqueCellId")

    return wide, cell_column, pc_columns


def build_neighbour_graph(wide, pc_columns, n_neighbors):
    """Construct an AnnData object over the PC matrix and compute its kNN graph."""
    adata = sc.AnnData(wide.select(pc_columns).to_numpy())
    sc.pp.neighbors(adata, use_rep="X", n_neighbors=n_neighbors)
    return adata


def perform_clustering(adata, wide, cell_column, leiden_resolution):
    """Run Leiden clustering and attach the assignments to the cell identities.

    AnnData preserves row order, so the i-th observation corresponds to the i-th
    row of `wide` and the two can be zipped without a join.
    """
    sc.tl.leiden(
        adata,
        resolution=leiden_resolution,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )

    clusters = pl.Series("Cluster", ("CL-" + adata.obs["leiden"].astype(str)).to_numpy())

    return wide.select(
        "UniqueCellId",
        pl.col(SAMPLE_HEADER).alias("SampleId"),
        pl.col(cell_column).alias("CellId"),
    ).with_columns(clusters)


def main():
    parser = argparse.ArgumentParser(
        description="Run Leiden clustering on PCA embeddings with duplicate CellIds across samples."
    )
    parser.add_argument(
        "--input_parquet",
        type=str,
        required=True,
        help="Path to the long-format PCA embeddings Parquet file.",
    )
    parser.add_argument(
        "--output_parquet",
        type=str,
        required=True,
        help="Path to save cluster assignments as Parquet.",
    )
    parser.add_argument(
        "--linker_parquet",
        type=str,
        default="leiden_linker.parquet",
        help="Path to save linker data for pFrame construction.",
    )
    parser.add_argument(
        "--n_neighbors",
        type=int,
        default=15,
        help="Number of neighbors for the graph (default: 15).",
    )
    parser.add_argument(
        "--leiden_resolution",
        type=float,
        default=1.0,
        help="Resolution for Leiden clustering (default: 1.0).",
    )

    args = parser.parse_args()

    wide, cell_column, pc_columns = load_embeddings(args.input_parquet)
    adata = build_neighbour_graph(wide, pc_columns, args.n_neighbors)
    cluster_assignments = perform_clustering(
        adata, wide, cell_column, args.leiden_resolution
    )

    # Linker column: [SampleId][CellId][Cluster] -> 1. Int32 matches the "Int"
    # valueType the workflow declares for the linker p-column.
    linker_data = cluster_assignments.select("SampleId", "CellId", "Cluster").with_columns(
        pl.lit(1, dtype=pl.Int32).alias("Link")
    )

    cluster_assignments.write_parquet(args.output_parquet)
    linker_data.write_parquet(args.linker_parquet)


if __name__ == "__main__":
    main()
