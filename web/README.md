# DeKnot in the browser

A static web page that runs DeKnot client-side via WebAssembly. Nothing is uploaded
anywhere: the FASTQ is read into the browser, assembled in a Web Worker, and the
walks are offered as a download.

## Files

| File | Purpose |
|---|---|
| `wasm_main.cpp` | C entry point (`dk_wasm_assemble`) called from JavaScript; writes `final_walks.fa` and a read→haplotype TSV into Emscripten's in-memory filesystem |
| `index.html` | The page: file picker, parameters, results table, downloads |
| `worker.js` | Loads `deknot.js` in a Web Worker so the page stays responsive |
| `test_wasm.js` | Node script to run the wasm build on a FASTQ and compare with the native binary |
| `deknot.js`, `deknot.wasm` | Build outputs (git-ignored) |

## Build

Install [Emscripten](https://emscripten.org/docs/getting_started/downloads.html)
(`brew install emscripten` on macOS), then from the repo root:

```
make wasm
```

## Run locally

Browsers refuse to load `.wasm` from `file://`, so serve the directory:

```
cd web && python3 -m http.server 8000
```

and open http://localhost:8000.

## Check against the native binary

```
./deknot locus.fq out_native
node web/test_wasm.js locus.fq out_wasm_walks.fa
diff out_native/final_walks.fa out_wasm_walks.fa
```

## Deploy

Pushing to the `web` branch runs `.github/workflows/pages.yml`, which builds with
Emscripten and publishes `web/` to GitHub Pages. Enable it once under
Settings → Pages → Source: "GitHub Actions".
