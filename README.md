# Table of Contents
* [Usage](#usage)
* [Parameters](#parameters)

<a name="usage"></a>
## Usage
```
deknot [input.fq] [output directory] <flags>
```

<a name="parameters"></a>
## Parameters
### Sketching
* `-s`: s value when selecting syncmers
### Graph growth
* `--k-prime`: initial number of syncmers per k-syncmer
* `--freq-thresh`: unitigs under this frequency will be dropped during graph compaction
* `--max-iterations`: maximum number of growth iterations before forcing DeKnot to terminate
