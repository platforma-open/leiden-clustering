<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlAlert, PlBlockPage, PlBtnGhost, PlDropdownRef, PlMaskIcon24, PlNumberField, PlSlideModal } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';
import { ref } from 'vue';

const app = useApp();

// const settingsAreShown = ref(app.model.outputs.UMAPPf === undefined)
const settingsAreShown = ref(true);
const showSettings = () => {
  settingsAreShown.value = true;
};

</script>

<template>
  <PlBlockPage>
    <template #title>Leiden Clustering</template>
    <template #append>
      <PlBtnGhost @click.stop="showSettings">
        Settings
        <template #append>
          <PlMaskIcon24 name="settings" />
        </template>
      </PlBtnGhost>
    </template>

    <PlSlideModal v-model="settingsAreShown">
      <template #title>Settings</template>
      <PlDropdownRef
        v-model="app.model.args.principalComponentsRef" :options="app.model.outputs.embeddingOptions"
        label="Select dataset"
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
