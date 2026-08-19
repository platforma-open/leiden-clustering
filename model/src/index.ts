import type { GraphMakerState } from "@milaboratories/graph-maker";
import type {
  InferOutputsType,
  PColumnIdAndSpec,
  PFrameHandle,
  PlRef,
  TreeNodeAccessor,
} from "@platforma-sdk/model";
import { BlockModel, isPColumn, isPColumnSpec } from "@platforma-sdk/model";

export type UiState = {
  graphStateUMAP: GraphMakerState;
  graphStateTSNE: GraphMakerState;
  anchorColumn?: PlRef;
};

export type BlockArgs = {
  principalComponentsRef?: PlRef;
  resolution: number;
  title?: string;
};

export const platforma = BlockModel.create()

  .withArgs<BlockArgs>({
    resolution: 0.5,
  })

  .withUiState<UiState>({
    graphStateUMAP: {
      title: "UMAP",
      template: "dots",
      currentTab: "settings",
    },
    graphStateTSNE: {
      title: "tSNE",
      template: "dots",
      currentTab: null,
    },
  })

  .argsValid((ctx) => ctx.args.principalComponentsRef !== undefined)

  .output("embeddingOptions", (ctx) =>
    ctx.resultPool.getOptions(
      (spec) => isPColumnSpec(spec) && spec.name === "pl7.app/rna-seq/pcvalue",
      { includeNativeLabel: false, addLabelAsSuffix: true },
    ),
  )

  .outputWithStatus("UMAPPf", (ctx): PFrameHandle | undefined => {
    // This block's own clusters first. Until the workflow has produced them there is
    // nothing to plot, and returning here keeps the lambda clear of the result pool.
    // `getData()` resolves data for every entry in the pool and marks the whole lambda
    // unstable while that is in flight; GraphMaker renders an unstable `undefined` as
    // "running". Reading the pool before this guard therefore made the chart claim the
    // block was running the moment a dataset was picked. Guarding first yields a stable
    // `undefined`, which GraphMaker renders as its idle "no pFrame" screen.
    const clusters = ctx.outputs?.resolve("leidenClusters")?.getPColumns();
    if (clusters === undefined) return undefined;

    // Get input data, to discern batch corrected or not
    if (!ctx.uiState?.anchorColumn) return undefined;
    const anchorSpec = ctx.resultPool.getPColumnSpecByRef(ctx.uiState?.anchorColumn);
    if (!anchorSpec) return undefined;

    const pCols = ctx.resultPool
      .getData()
      .entries.map((c) => c.obj)
      .filter(isPColumn<TreeNodeAccessor>)
      .filter((col) => {
        return (
          (col.spec.name === "pl7.app/rna-seq/umap1" ||
            col.spec.name === "pl7.app/rna-seq/umap2" ||
            col.spec.name === "pl7.app/rna-seq/umap3") &&
          col.spec.domain?.["pl7.app/blockId"] === anchorSpec.domain?.["pl7.app/blockId"]
        );
      });

    // Return batch corrected UMAP if present
    let finalPcols = pCols.filter(
      (col) => col.spec.domain?.["pl7.app/rna-seq/batch-corrected"] === "true",
    );
    if (finalPcols.length === 0) {
      finalPcols = pCols.filter(
        (col) => col.spec.domain?.["pl7.app/rna-seq/batch-corrected"] === "false",
      );
    }

    return ctx.createPFrame([...finalPcols, ...clusters]);
  })

  .outputWithStatus("tSNEPf", (ctx): PFrameHandle | undefined => {
    // See UMAPPf: guard on this block's own clusters before touching the result pool.
    const clusters = ctx.outputs?.resolve("leidenClusters")?.getPColumns();
    if (clusters === undefined) return undefined;

    // Get input data, to discern batch corrected or not
    if (!ctx.uiState?.anchorColumn) return undefined;
    const anchorSpec = ctx.resultPool.getPColumnSpecByRef(ctx.uiState?.anchorColumn);
    if (!anchorSpec) return undefined;

    const pCols = ctx.resultPool
      .getData()
      .entries.map((c) => c.obj)
      .filter(isPColumn<TreeNodeAccessor>)
      .filter((col) => {
        return (
          (col.spec.name === "pl7.app/rna-seq/tsne1" ||
            col.spec.name === "pl7.app/rna-seq/tsne2" ||
            col.spec.name === "pl7.app/rna-seq/tsne3") &&
          col.spec.domain?.["pl7.app/blockId"] === anchorSpec.domain?.["pl7.app/blockId"]
        );
      });

    // Return batch corrected UMAP if present
    let finalPcols = pCols.filter(
      (col) => col.spec.domain?.["pl7.app/rna-seq/batch-corrected"] === "true",
    );
    if (finalPcols.length === 0) {
      finalPcols = pCols.filter(
        (col) => col.spec.domain?.["pl7.app/rna-seq/batch-corrected"] === "false",
      );
    }

    return ctx.createPFrame([...finalPcols, ...clusters]);
  })

  .output("plotPcols", (ctx) => {
    // See UMAPPf: guard on this block's own clusters before touching the result pool.
    const clusters = ctx.outputs?.resolve("leidenClusters")?.getPColumns();
    if (clusters === undefined) return undefined;

    // Get input data, to discern batch corrected or not
    if (!ctx.uiState?.anchorColumn) return undefined;
    const anchorSpec = ctx.resultPool.getPColumnSpecByRef(ctx.uiState?.anchorColumn);
    if (!anchorSpec) return undefined;

    const pCols = ctx.resultPool
      .getData()
      .entries.map((c) => c.obj)
      .filter(isPColumn<TreeNodeAccessor>)
      .filter((col) => {
        return (
          (col.spec.name.slice(0, -1) === "pl7.app/rna-seq/tsne" ||
            col.spec.name.slice(0, -1) === "pl7.app/rna-seq/umap") &&
          col.spec.domain?.["pl7.app/blockId"] === anchorSpec.domain?.["pl7.app/blockId"]
        );
      });

    // Return batch corrected UMAP if present
    let finalPcols = pCols.filter(
      (col) => col.spec.domain?.["pl7.app/rna-seq/batch-corrected"] === "true",
    );
    if (finalPcols.length === 0) {
      finalPcols = pCols.filter(
        (col) => col.spec.domain?.["pl7.app/rna-seq/batch-corrected"] === "false",
      );
    }

    return [...finalPcols, ...clusters].map(
      (c) =>
        ({
          columnId: c.id,
          spec: c.spec,
        }) satisfies PColumnIdAndSpec,
    );
  })

  .output("isRunning", (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .sections((_ctx) => [{ type: "link", href: "/", label: "Main" }])

  .title((ctx) => (ctx.args.title ? `Leiden Clustering - ${ctx.args.title}` : "Leiden Clustering"))

  .done(2);

export type BlockOutputs = InferOutputsType<typeof platforma>;
