<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlAlert, PlBlockPage, PlDropdownRef, PlNumberField } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';
import type { PlRef } from '@platforma-sdk/model';
import { plRefsEqual } from '@platforma-sdk/model';

const app = useApp();

function setInput(inputRef?: PlRef) {
  app.model.args.principalComponentsRef = inputRef;
  if (inputRef)
    app.model.args.title = app.model.outputs.embeddingOptions?.find((o) => plRefsEqual(o.ref, inputRef))?.label;
  else
    app.model.args.title = undefined;
}

</script>

<template>
  <PlBlockPage>
    <template #title>Settings</template>
    <PlDropdownRef
      v-model="app.model.args.principalComponentsRef" :options="app.model.outputs.embeddingOptions"
      :style="{ width: '320px' }"
      label="Select dataset"
      clearable @update:model-value="setInput"
    />
    <PlNumberField
      v-model="app.model.args.resolution"
      :style="{ width: '320px' }"
      label="Resolution" :minValue="0.1" :step="0.1"
    >
      <template #tooltip>
        Select resolution for Leiden clustering. The bigger the resolution, the more clusters will be found.
      </template>
    </PlNumberField>
    <PlAlert v-if="app.model.args.resolution > 1" type="warn" :style="{ width: '320px' }">
      {{ "Warning: The selected resolution is over commonly used range (0.4 - 1.0)" }}
    </PlAlert>
    <PlAlert v-if="app.model.args.resolution < 0.4" type="warn" :style="{ width: '320px' }">
      {{ "Warning: The selected resolution is under commonly used range (0.4 - 1.0)" }}
    </PlAlert>
  </PlBlockPage>
</template>
