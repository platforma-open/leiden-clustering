<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlAlert, PlBlockPage, PlDropdownRef, PlNumberField, PlRow, PlTabs } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';

import type { PredefinedGraphOption } from '@milaboratories/graph-maker';
import { GraphMaker } from '@milaboratories/graph-maker';
import type { PColumnIdAndSpec, PlRef } from '@platforma-sdk/model';
import { plRefsEqual } from '@platforma-sdk/model';
import { computed, reactive } from 'vue';

const app = useApp();

const data = reactive({
  currentTab: 'umap',
});

const tabOptions = [
  { label: 'UMAP', value: 'umap' },
  { label: 't-SNE', value: 'tsne' },
];

function setInput(inputRef?: PlRef) {
  app.model.ui.anchorColumn = inputRef;
  if (inputRef)
    app.model.args.title = app.model.outputs.embeddingOptions?.find((o) => plRefsEqual(o.ref, inputRef))?.label;
  else
    app.model.args.title = undefined;
}

function getIndex(name: string, pcols: PColumnIdAndSpec[]): number {
  return pcols.findIndex((p) => (p.spec.name === name
  ));
}

/* Function to create default options according to the selected tab */
function createDefaultOptions(
  pcols: PColumnIdAndSpec[] | undefined,
  coord1Name: string,
  coord2Name: string,
): PredefinedGraphOption<'scatterplot-umap'>[] | undefined {
  if (!pcols || pcols.length === 0)
    return undefined;

  const coord1Index = getIndex(coord1Name, pcols);
  const coord2Index = getIndex(coord2Name, pcols);
  const leidenIndex = getIndex('pl7.app/rna-seq/leidencluster', pcols);

  if (coord1Index === -1 || coord2Index === -1)
    return undefined;

  const defaults: PredefinedGraphOption<'scatterplot-umap'>[] = [
    {
      inputName: 'x',
      selectedSource: pcols[coord1Index].spec,
    },
    {
      inputName: 'y',
      selectedSource: pcols[coord2Index].spec,
    },
    {
      inputName: 'grouping',
      selectedSource: pcols[leidenIndex].spec,
    },
  ];

  return defaults;
}

const defaultOptions = computed((): PredefinedGraphOption<'scatterplot-umap'>[] | undefined => {
  if (data.currentTab === 'umap') {
    return createDefaultOptions(
      app.model.outputs.plotPcols,
      'pl7.app/rna-seq/umap1',
      'pl7.app/rna-seq/umap2',
    );
  }
  if (data.currentTab === 'tsne') {
    return createDefaultOptions(
      app.model.outputs.plotPcols,
      'pl7.app/rna-seq/tsne1',
      'pl7.app/rna-seq/tsne2',
    );
  }
  return undefined;
});

/* Modify graph state, pframe and default options based on the selected tab */
const graphState = computed({
  get: () => data.currentTab === 'umap' ? app.model.ui.graphStateUMAP : app.model.ui.graphStateTSNE,
  set: (value) => {
    if (data.currentTab === 'umap')
      app.model.ui.graphStateUMAP = value;
    else
      app.model.ui.graphStateTSNE = value;
  },
});

const pFrame = computed(() => data.currentTab === 'umap' ? app.model.outputs.UMAPPf : app.model.outputs.tSNEPf);

/* Use both currentTab and pFrame in :key to force re-render the graph when either args (which changes the pFrame) or the tab changes */

</script>

<template>
  <PlBlockPage>
    <GraphMaker
      :key="`${data.currentTab}-${pFrame}`"
      v-model="graphState"
      chartType="scatterplot-umap"
      :p-frame="pFrame"
      :default-options="defaultOptions"
    >
      <template #titleLineSlot>
        <PlTabs v-model="data.currentTab" :options="tabOptions" :style="{ display: 'flex', justifyContent: 'flex-end' }"/>
      </template>
      <template #settingsSlot>
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
