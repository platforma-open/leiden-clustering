<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlBlockPage } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';

import type { PredefinedGraphOption } from '@milaboratories/graph-maker';
import { GraphMaker } from '@milaboratories/graph-maker';
import type { PColumnIdAndSpec } from '@platforma-sdk/model';
import { ref, watch } from 'vue';

const app = useApp();

function getDefaultOptions(plotPcols?: PColumnIdAndSpec[]) {
  if (!plotPcols) {
    return undefined;
  }

  function getIndex(name: string, pcols: PColumnIdAndSpec[]): number {
    return pcols.findIndex((p) => p.spec.name === name,
    );
  }
  const defaults: PredefinedGraphOption<'scatterplot-umap'>[] = [
    {
      inputName: 'x',
      selectedSource: plotPcols[getIndex('pl7.app/rna-seq/umap1',
        plotPcols)].spec,
    },
    {
      inputName: 'y',
      selectedSource: plotPcols[getIndex('pl7.app/rna-seq/umap2',
        plotPcols)].spec,
    },
    {
      inputName: 'grouping',
      selectedSource: plotPcols[getIndex('pl7.app/rna-seq/leidencluster',
        plotPcols)].spec,
    },
  ];
  return defaults;
};

// Steps needed to reset graph maker after changing input table
const defaultOptions = ref(getDefaultOptions(app.model.outputs.plotPcols));
const key = ref(defaultOptions.value ? JSON.stringify(defaultOptions.value) : '');
// Reset graph maker state to allow new selection of defaults
watch(() => app.model.args.principalComponentsRef, (_) => {
  delete app.model.ui.graphStateUMAP.optionsState;
  defaultOptions.value = getDefaultOptions(app.model.outputs.plotPcols);
  key.value = defaultOptions.value ? JSON.stringify(defaultOptions.value) : '';
},
);

</script>

<template>
  <PlBlockPage>
    <GraphMaker
      :key="key"
      v-model="app.model.ui.graphStateUMAP" chartType="scatterplot-umap" :p-frame="app.model.outputs.UMAPPf"
      :default-options="defaultOptions"
    />
  </PlBlockPage>
</template>
