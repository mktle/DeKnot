// Runs the wasm build on a FASTQ under Node and writes the walks FASTA, so it can be
// diffed against the native binary's final_walks.fa:
//   node web/test_wasm.js locus.fq out_walks.fa [out_assignments.tsv]
const fs = require('fs');
const path = require('path');
const createDeKnot = require(path.join(__dirname, 'deknot.js'));

(async () => {
  const [fq, outWalks, outAssign] = process.argv.slice(2);
  if (!fq || !outWalks) { console.error('usage: node test_wasm.js in.fq out_walks.fa [out_assign.tsv]'); process.exit(2); }
  const mod = await createDeKnot({ print: (t) => process.stdout.write(t + '\n'), printErr: (t) => process.stderr.write(t + '\n') });
  mod.FS.writeFile('/input.fq', fs.readFileSync(fq));
  const t0 = Date.now();
  const n = mod.ccall('dk_wasm_assemble', 'number',
    ['string', 'string', 'string', 'number', 'number', 'number', 'number', 'number', 'number'],
    ['/input.fq', '/walks.fa', '/assign.tsv', 9, 2, 4, 2000, 0, process.env.DK_VERBOSE ? 1 : 0]);
  console.error(`walks=${n} time=${((Date.now() - t0) / 1000).toFixed(1)}s`);
  if (n < 0) process.exit(1);
  fs.writeFileSync(outWalks, mod.FS.readFile('/walks.fa'));
  if (outAssign) fs.writeFileSync(outAssign, mod.FS.readFile('/assign.tsv'));
})();
