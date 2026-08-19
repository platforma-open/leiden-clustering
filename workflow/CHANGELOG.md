# @platforma-open/milaboratories.leiden-clustering.workflow

## 1.7.1

### Patch Changes

- aa3241e: Migrate block onto the structurer (block-tools 2.14.2) — full SDK upgrade: model 1.82.0, ui-vue 1.82.1, workflow-tengo 6.8.2, tengo-builder 4.0.23, package-builder 3.15.0, test 1.82.3. Adopts the canonical tool-managed layout (oxlint/oxfmt, tsconfig, turbo, CI workflows, managed package.json + catalog), the mandatory sibling `kind/` package declaring the block's init-params contract, and the slim facade for the root block package. Author-code fixes for the SDK majors: explicit type argument on the `isPColumn` filters feeding `createPFrame`, removal of the retired `@platforma-sdk/ui-vue/styles` and `@milaboratories/graph-maker/styles` imports, `UMAPPf`/`tSNEPf` moved to `outputWithStatus` for GraphMaker's `p-frame` prop, and the model export renamed `model` -> `platforma` for the facade. Two dead test files were removed and the dependency graph deduped to single copies of `@platforma-sdk/model`, `@platforma-sdk/ui-vue`, `zod` and `software-small-binaries`.

  Exchange the PCA embeddings and cluster assignments as Parquet instead of CSV, and give the embeddings conversion its own budget.

  The `xsv.exportFrame` step that materialises the long-format `pcvalue` p-column for the Python tool
  is the heaviest conversion in the block — it holds one row per (cell, principal component) pair,
  reaching tens of millions of rows on real datasets, while running with a hardcoded 16 GiB / 1 CPU.
  Two of its four columns are high-cardinality strings repeated once per principal component, which is
  the worst case for CSV text encoding; the upstream block already stores the column as Parquet, so
  the export was inflating columnar data into text for no benefit.

  - workflow: `xsv.exportFrame([embeddings], "parquet", ...)` with a dedicated 32 GiB / 4 CPU budget.
    The cluster and linker imports also move to Parquet but keep the 16 GiB / 1 CPU default — they
    carry one row per cell. The intermediate is renamed `csvEmbeddings` -> `embeddingsParquet` and
    `*.csv` -> `*.parquet` end-to-end, along with the child template's output keys.
  - software: `calculate_leiden_clusters.py` moves from pandas to polars, reading via
    `pl.scan_parquet` with the repeated sample/cell columns cast to `Categorical` inside the scan
    plan, and writing both outputs with `write_parquet`. Input dtypes are restored before writing so
    the emitted Parquet matches the axis specs pfconv reads it back against.
  - Fixes a latent id bug: the pivot now keys on the (sample, cell) pair instead of a
    `Sample + "_" + CellId` string that was split back apart on the first underscore, so a sample
    whose name contains an underscore no longer has its SampleId and CellId recovered incorrectly.

  The pandas -> polars rewrite itself does not change cluster assignments: the old and new scripts
  were run against the same 24592-cell input at several resolutions and produced byte-identical
  partitions (ARI 1.0000). The pivot ordering is pinned to the order pandas produced, so that stays
  true regardless of how the exchange file happens to be ordered.

  Cluster assignments are also preserved, but that took a second fix. The block pinned five Python
  packages and left `pynndescent`, `numba`, `llvmlite`, `umap-learn`, `scikit-learn`, `igraph` and
  `anndata` floating — all of which sit on the neighbours/Leiden path. Two builds at different dates
  resolved different versions and produced different clusterings with no code change at all. Measured
  on a real 24592-cell input at one resolution, three environments differing _only_ in these versions
  gave 20 / 18 / 18 clusters with ARI between them as low as 0.52.

  The full dependency set is now pinned to the versions the **previously published block actually
  ran**, recovered from its released container image, so this release reproduces the clusterings
  users already have rather than silently re-partitioning their data.

- Updated dependencies [aa3241e]
  - @platforma-open/milaboratories.leiden-clustering.software@1.6.4

## 1.7.0

### Minor Changes

- 71fb4c2: Enable block deduplication, fix labelling and update metadata.

## 1.6.5

### Patch Changes

- f6789a7: Update to Parquet

## 1.6.4

### Patch Changes

- ec60dad: technical release
- dcbd64b: technical release
- 69e135c: technical release
- ba6a3fe: technical release
- Updated dependencies [ec60dad]
- Updated dependencies [dcbd64b]
- Updated dependencies [69e135c]
- Updated dependencies [ba6a3fe]
  - @platforma-open/milaboratories.leiden-clustering.software@1.6.3

## 1.6.3

### Patch Changes

- Updated dependencies [9c10a43]
  - @platforma-open/milaboratories.leiden-clustering.software@1.6.2

## 1.6.2

### Patch Changes

- e80cf06: Updated trace

## 1.6.1

### Patch Changes

- Updated dependencies [b7a8ac6]
  - @platforma-open/milaboratories.leiden-clustering.software@1.6.1

## 1.6.0

### Minor Changes

- b7c9bfa: Add batch support

## 1.5.1

### Patch Changes

- Updated dependencies [7407462]
  - @platforma-open/milaboratories.leiden-clustering.software@1.6.0

## 1.5.0

### Minor Changes

- 213044c: Minor change to labels of clusters (from 1, 2, 3 to CL-1, CL-2, CL-3)

### Patch Changes

- Updated dependencies [213044c]
  - @platforma-open/milaboratories.leiden-clustering.software@1.5.0

## 1.4.0

### Minor Changes

- 925fc38: Update trace and importance

## 1.3.0

### Minor Changes

- 9e5428c: Add linker pcol for export

### Patch Changes

- Updated dependencies [9e5428c]
  - @platforma-open/milaboratories.leiden-clustering.software@1.4.0

## 1.2.0

### Minor Changes

- be07e24: Improved UI

## 1.1.3

### Patch Changes

- e322eef: chore: update deps

## 1.1.2

### Patch Changes

- Updated dependencies [de3457f]
  - @platforma-open/milaboratories.leiden-clustering.software@1.3.0

## 1.1.1

### Patch Changes

- ffcfae0: Updated dependencies
- Updated dependencies [ffcfae0]
  - @platforma-open/milaboratories.leiden-clustering.software@1.2.1

## 1.1.0

### Minor Changes

- 62bf70f: First block version

### Patch Changes

- Updated dependencies [62bf70f]
  - @platforma-open/milaboratories.leiden-clustering.software@1.2.0
