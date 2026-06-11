import { defineConfig } from 'vite';
import wasm from 'vite-plugin-wasm';
import path from 'path';

export default defineConfig({
  plugins: [
    wasm()
  ],
  build: {
    target: 'esnext'
  },
  resolve: {
    alias: {
      'brep': path.resolve(__dirname, '../../src')
    }
  },
  server: {
    fs: {
      allow: [
        '.',
        path.resolve(__dirname, '../../node_modules')
      ]
    }
  },
  optimizeDeps: {
    exclude: ['occt-wasm', 'brepjs-opencascade', 'brepjs'],
    rolldownOptions: {
      external: ['brepjs-opencascade']
    }
  }
});
