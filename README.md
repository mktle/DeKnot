# Table of Contents
* [Usage](#usage)
* [Parameters](#parameters)
* [Browser version](#browser)

<a name="usage"></a>
## Usage
```
deknot [input.fq] [output directory] <flags>
```

<a name="parameters"></a>
## Parameters
### General
* `-v`: verbose mode with logging statements
### Sketching
* `-s`: s value when selecting syncmers
### Graph growth
* `--k-prime`: initial number of syncmers per k-syncmer
* `--freq-thresh`: unitigs under this frequency will be dropped during graph compaction
* `--max-iterations`: maximum number of growth iterations before forcing DeKnot to terminate

<a name="browser"></a>
## Browser version
DeKnot can also run entirely in a web browser (compiled to WebAssembly), with no upload: see [web/README.md](web/README.md).
Build it with `make wasm` after installing Emscripten.
