<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlAlert, PlBlockPage, PlDropdownRef, PlNumberField, PlRow } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';

import type { PredefinedGraphOption } from '@milaboratories/graph-maker';
import { GraphMaker } from '@milaboratories/graph-maker';
import type { PColumnIdAndSpec } from '@platforma-sdk/model';
import { plRefsEqual, type PlRef } from '@platforma-sdk/model';
import { ref, watch } from 'vue';

const app = useApp();
const settingsOpen = ref(true);

function setInput(inputRef?: PlRef) {
  app.model.ui.anchorColumn = inputRef;
  if (inputRef)
    app.model.args.title = app.model.outputs.embeddingOptions?.find((o) => plRefsEqual(o.ref, inputRef))?.label;
  else
    app.model.args.title = undefined;
}

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
      v-model="app.model.ui.graphStateUMAP"
      chartType="scatterplot-umap"
      :p-frame="app.model.outputs.UMAPPf"
      :default-options="defaultOptions"
      @run="settingsOpen = false"
    >
      <template v-if="settingsOpen" #settingsSlot>
        <PlDropdownRef
          v-model="app.model.args.principalComponentsRef"
          :options="app.model.outputs.embeddingOptions"
          :style="{ width: '320px' }"
          label="Select dataset"
          clearable
          required
          @update:model-value="setInput"
        />
        <PlRow>
          <PlNumberField
            v-model="app.model.args.resolution"
            :style="{ width: '320px' }"
            label="Resolution"
            :minValue="0.1"
            :step="0.1"
          >
            <template #tooltip>
              <div>
                <strong>Leiden Clustering Resolution</strong><br/>
                Controls the granularity of clustering. Higher values result in more clusters.<br/><br/>
                <strong>Recommended ranges:</strong><br/>
                • 0.4-0.8: Optimal for most datasets<br/>
                • 0.1-0.4: For fewer, larger clusters<br/>
                • 0.8-1.2: For more, smaller clusters<br/><br/>
                <strong>Effect:</strong> Higher resolution creates more granular clusters, lower resolution creates broader clusters.
              </div>
            </template>
          </PlNumberField>
        </PlRow>
        <!-- Add warnings if selected parameters are out of most commonly used bounds -->
        <PlAlert v-if="app.model.args.resolution > 1" type="warn" :style="{ width: '320px' }">
          <template #title>High Resolution</template>
          The selected resolution ({{ app.model.args.resolution }}) is above the recommended range (0.4-1.0).
          This may result in over-clustering with too many small clusters.
        </PlAlert>
        <PlAlert v-if="app.model.args.resolution < 0.4" type="warn" :style="{ width: '320px' }">
          <template #title>Low Resolution</template>
          The selected resolution ({{ app.model.args.resolution }}) is below the recommended range (0.4-1.0).
          This may result in under-clustering with too few large clusters.
        </PlAlert>
      </template>
    </GraphMaker>
  </PlBlockPage>
</template>
