---
'@platforma-open/milaboratories.leiden-clustering.model': patch
'@platforma-open/milaboratories.leiden-clustering.ui': patch
'@platforma-open/milaboratories.leiden-clustering.workflow': patch
'@platforma-open/milaboratories.leiden-clustering.software': patch
'@platforma-open/milaboratories.leiden-clustering': patch
---

Migrate block onto the structurer (block-tools 2.13.0) — full SDK upgrade: model/ui-vue 1.81.1, workflow-tengo 6.8.2, tengo-builder 4.0.22, package-builder 3.15.0, test 1.81.3. Adopts the canonical tool-managed layout (oxlint/oxfmt, tsconfig, turbo, CI workflows, managed package.json + catalog) and the slim facade for the root block package. Author-code fixes for the SDK majors: explicit type argument on the `isPColumn` filters feeding `createPFrame`, removal of the retired `@platforma-sdk/ui-vue/styles` and `@milaboratories/graph-maker/styles` imports, `UMAPPf`/`tSNEPf` moved to `outputWithStatus` for GraphMaker's `p-frame` prop, and the model export renamed `model` -> `platforma` for the facade. Two dead test files were removed and the dependency graph deduped to single copies of `@platforma-sdk/model`, `@platforma-sdk/ui-vue`, `zod` and `software-small-binaries`.

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
