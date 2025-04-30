import type { GraphMakerState } from '@milaboratories/graph-maker';
import type {
  InferOutputsType,
  PColumnIdAndSpec,
  PFrameHandle,
  PlRef } from '@platforma-sdk/model';
import {
  BlockModel,
  isPColumn,
  isPColumnSpec,
} from '@platforma-sdk/model';

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

export const model = BlockModel.create()

  .withArgs<BlockArgs>({
    resolution: 1,
  })

  .withUiState<UiState>({
    graphStateUMAP: {
      title: 'UMAP',
      template: 'dots',
      currentTab: null,
    },
    graphStateTSNE: {
      title: 'tSNE',
      template: 'dots',
      currentTab: null,
    },
  })

  .argsValid((ctx) => ctx.args.principalComponentsRef !== undefined)

  .output('embeddingOptions', (ctx) =>
    ctx.resultPool.getOptions((spec) => isPColumnSpec(spec)
      && spec.name === 'pl7.app/rna-seq/pcvalue'
    , { includeNativeLabel: true, addLabelAsSuffix: true }),
  )

  .output('UMAPPf', (ctx): PFrameHandle | undefined => {
    // Get input data, to discern batch corrected or not
    if (!ctx.uiState?.anchorColumn) return undefined;
    const anchorSpec = ctx.resultPool.getPColumnSpecByRef(ctx.uiState?.anchorColumn);
    if (!anchorSpec) return undefined;

    const pCols
      = ctx.resultPool
        .getData()
        .entries.map((c) => c.obj)
        .filter(isPColumn)
        .filter((col) => {
          return (col.spec.name === 'pl7.app/rna-seq/umap1'
            || col.spec.name === 'pl7.app/rna-seq/umap2'
            || col.spec.name === 'pl7.app/rna-seq/umap3')
          && col.spec.domain?.['pl7.app/rna-seq/batch-corrected'] === anchorSpec.domain?.['pl7.app/rna-seq/batch-corrected'];
        });

    // enriching with leiden clusters data
    const upstream
      = ctx.outputs?.resolve('leidenClusters')?.getPColumns();

    if (upstream === undefined) {
      return undefined;
    }

    return ctx.createPFrame([...pCols, ...upstream]);
  })

  .output('tSNEPf', (ctx): PFrameHandle | undefined => {
    // Get input data, to discern batch corrected or not
    if (!ctx.uiState?.anchorColumn) return undefined;
    const anchorSpec = ctx.resultPool.getPColumnSpecByRef(ctx.uiState?.anchorColumn);
    if (!anchorSpec) return undefined;

    const pCols
      = ctx.resultPool
        .getData()
        .entries.map((c) => c.obj)
        .filter(isPColumn)
        .filter((col) => {
          return (col.spec.name === 'pl7.app/rna-seq/tsne1'
            || col.spec.name === 'pl7.app/rna-seq/tsne2'
            || col.spec.name === 'pl7.app/rna-seq/tsne3')
          && col.spec.domain?.['pl7.app/rna-seq/batch-corrected'] === anchorSpec.domain?.['pl7.app/rna-seq/batch-corrected'];
        });

    // enriching with leiden clusters data
    const upstream
      = ctx.outputs?.resolve('leidenClusters')?.getPColumns();

    if (upstream === undefined) {
      return undefined;
    }

    return ctx.createPFrame([...pCols, ...upstream]);
  })

  .output('plotPcols', (ctx) => {
    // Get input data, to discern batch corrected or not
    if (!ctx.uiState?.anchorColumn) return undefined;
    const anchorSpec = ctx.resultPool.getPColumnSpecByRef(ctx.uiState?.anchorColumn);
    if (!anchorSpec) return undefined;

    const pCols
      = ctx.resultPool
        .getData()
        .entries.map((c) => c.obj)
        .filter(isPColumn)
        .filter((col) => {
          return ((col.spec.name.slice(0, -1) === 'pl7.app/rna-seq/tsne'
            || col.spec.name.slice(0, -1) === 'pl7.app/rna-seq/umap')
          && col.spec.domain?.['pl7.app/rna-seq/batch-corrected'] === anchorSpec.domain?.['pl7.app/rna-seq/batch-corrected']);
        });

    // enriching with leiden clusters data
    const upstream
      = ctx.outputs?.resolve('leidenClusters')?.getPColumns();

    if (upstream === undefined) {
      return undefined;
    }

    return [...pCols, ...upstream].map(
      (c) =>
        ({
          columnId: c.id,
          spec: c.spec,
        } satisfies PColumnIdAndSpec),
    );
  })

  .output('isRunning', (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .sections((_ctx) => ([
    { type: 'link', href: '/', label: 'Main' },
    { type: 'link', href: '/umap', label: 'UMAP' },
    { type: 'link', href: '/tsne', label: 'tSNE' },
  ]))

  .title((ctx) =>
    ctx.args.title
      ? `Leiden Clustering - ${ctx.args.title}`
      : 'Leiden Clustering',
  )

  .done();

export type BlockOutputs = InferOutputsType<typeof model>;
