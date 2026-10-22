// Web Worker that runs the DeKnot WebAssembly module off the main thread.
// Messages in:  { type: 'run', fastq: ArrayBuffer, name: string, params: {...} }
// Messages out: { type: 'log', text }  |  { type: 'done', ... }  |  { type: 'error', message }

importScripts('deknot.js');

self.onmessage = async (ev) => {
  const msg = ev.data;
  if (msg.type !== 'run') return;
  const { fastq, params } = msg;
  const log = (text) => self.postMessage({ type: 'log', text });

  let mod;
  try {
    // A fresh module instance per run: no state carries over and memory is released.
    mod = await createDeKnot({ print: log, printErr: log });
  } catch (e) {
    self.postMessage({ type: 'error', message: 'Could not load deknot.wasm: ' + (e && e.message ? e.message : e) });
    return;
  }

  const t0 = performance.now();
  try {
    mod.FS.writeFile('/input.fq', new Uint8Array(fastq));
    self.postMessage({ type: 'status', text: 'Assembling…' });
    const n = mod.ccall(
      'dk_wasm_assemble', 'number',
      ['string', 'string', 'string', 'number', 'number', 'number', 'number', 'number', 'number'],
      ['/input.fq', '/final_walks.fa', '/read_assignments.tsv',
       params.s, params.kPrime, params.freqThresh, params.maxIterations,
       params.checkStrand ? 1 : 0, params.verbose ? 1 : 0]
    );
    if (n === -1) throw new Error('No reads were parsed from the input. Is it an uncompressed FASTQ?');
    if (n < 0) throw new Error('Assembly failed (code ' + n + ').');
    const walks = mod.FS.readFile('/final_walks.fa', { encoding: 'utf8' });
    const assignments = mod.FS.readFile('/read_assignments.tsv', { encoding: 'utf8' });
    self.postMessage({
      type: 'done',
      nWalks: n,
      walks,
      assignments,
      seconds: (performance.now() - t0) / 1000,
    });
  } catch (e) {
    self.postMessage({ type: 'error', message: (e && e.message) ? e.message : String(e) });
  }
};
