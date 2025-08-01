import { model } from '@platforma-open/milaboratories.leiden-clustering.model';
import { defineApp } from '@platforma-sdk/ui-vue';
import MainPage from './pages/MainPage.vue';
import tSNE from './pages/tSNE.vue';

export const sdkPlugin = defineApp(model, () => {
  return {
    routes: {
      '/': () => MainPage,
      '/tsne': () => tSNE,
    },
  };
});

export const useApp = sdkPlugin.useApp;
