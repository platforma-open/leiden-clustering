<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlAlert, PlBlockPage, PlBtnGhost, PlDropdownRef, PlMaskIcon24, PlNumberField, PlSlideModal } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';
import { reactive } from 'vue';
import type { PlRef } from '@platforma-sdk/model';
import { plRefsEqual } from '@platforma-sdk/model';

const app = useApp();

const data = reactive<{
  settingsOpen: boolean;
}>({
  settingsOpen: app.model.args.principalComponentsRef === undefined,
});

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
    <template #title>Leiden Clustering</template>
    <template #append>
      <PlBtnGhost @click.stop="() => data.settingsOpen = true">
        Settings
        <template #append>
          <PlMaskIcon24 name="settings" />
        </template>
      </PlBtnGhost>
    </template>

    <PlSlideModal v-model="data.settingsOpen">
      <template #title>Settings</template>
      <PlDropdownRef
        v-model="app.model.args.principalComponentsRef" :options="app.model.outputs.embeddingOptions"
        label="Select dataset"
        clearable @update:model-value="setInput"
      />
      <PlNumberField
        v-model="app.model.args.resolution"
        label="Resolution" :minValue="0.1" :step="0.1"
      >
        <template #tooltip>
          Select resolution for Leiden clustering. The bigger the resolution, the more clusters will be found.
        </template>
      </PlNumberField>
      <PlAlert v-if="app.model.args.resolution > 1" type="warn">
        {{ "Warning: The selected resolution is over commonly used range (0.4 - 1.0)" }}
      </PlAlert>
      <PlAlert v-if="app.model.args.resolution < 0.4" type="warn">
        {{ "Warning: The selected resolution is under commonly used range (0.4 - 1.0)" }}
      </PlAlert>
    </PlSlideModal>
  </PlBlockPage>
</template>
