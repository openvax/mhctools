// Regenerate descriptor references from the checksum-verified original app.
// Prediction and routine tests do not require Node or execute this source.
const fs = require('fs');
const path = require('path');
const vm = require('vm');
const crypto = require('crypto');

if (process.argv.length !== 4) {
  throw Error('Usage: node scripts/record_cavaco_reference.js APP_SOURCE_DIRECTORY OUTPUT_JSON');
}
const root = path.resolve(process.argv[2]);
const hashes = {
  'calc.js': 'e780b0e90e8438fc700ecef1a8b681d97704949eb9b6edf8be630e6c18f19822',
  'amino.json': 'e0806f668cc9d98ac8d2163d0516392e3fb3cf349397018698801b1468e39a24',
  'codes.json': 'f02e2fbce52b00e6916830ed2a23c58582a9b8c76dbd224926dd440e0b767f3d',
  'package.json': 'e00f4372343f8a00eb8880d894ddf34657a2f6f652cd1be2e5d4bd4881e8b5a0',
};
const files = {};
for (const [name, expected] of Object.entries(hashes)) {
  const raw = fs.readFileSync(path.join(root, name));
  if (crypto.createHash('sha256').update(raw).digest('hex') !== expected) {
    throw Error(`Source checksum mismatch: ${name}`);
  }
  files[name] = raw.toString('utf8');
}
const manifest = JSON.parse(files['package.json']);
if (manifest.license !== 'CC0-1.0' || manifest.version !== '1.0.0') {
  throw Error('Unexpected app identity/license');
}
const amino = JSON.parse(files['amino.json']);
const elements = {};
function element() {
  return {children: [{children: [{}]}, {style: {}}], style: {},
    appendChild() {}, getElementsByClassName() { return []; }};
}
const document = {
  addEventListener() {}, createElement: element,
  getElementById(id) { return elements[id] ||= element(); },
};
const sandbox = {__dirname: root, document, amino, output: {},
  require(name) {
    if (name !== 'fs') throw Error(`Forbidden module: ${name}`);
    return {readFileSync(filename) {
      if (filename !== path.join(root, 'data', 'codes.json')) {
        throw Error(`Forbidden read: ${filename}`);
      }
      return files['codes.json'];
    }};
  },
};
vm.createContext(sandbox);
// Capture numeric results without editing the checksum-verified app source.
vm.runInContext(`const nativeFixed = Number.prototype.toFixed;
  Number.prototype.toFixed = function(digits) {
    if (digits === 1) output.pi = Number(this);
    return nativeFixed.call(this, digits);
  };
  const nativeExp = Math.exp;
  Math.exp = function(x) { output.ln_app = x; return nativeExp(x); };`, sandbox);
vm.runInContext(files['calc.js'], sandbox, {timeout: 1000});
vm.runInContext('aminos = amino.list; bind = amino.bind;', sandbox);
const reference = JSON.parse(fs.readFileSync(
  path.join(__dirname, '..', 'tests', 'data', 'cavaco_source_reference.json'), 'utf8'));
reference.records = reference.records.map(row => {
  sandbox.peptide = row.peptide;
  sandbox.output = {};
  vm.runInContext('runSequence(peptide, 1, "none")', sandbox, {timeout: 1000});
  const np = Array.from(row.peptide).filter(
    aa => amino.list[aa].composition_code === 3).length / row.peptide.length * 100;
  return {...row, pi: sandbox.output.pi, np_percent: np,
    ln_app: sandbox.output.ln_app,
    app_minutes: elements['row-hflf'].children[1].textContent,
    ln_published: 2.226 + 0.053 * np - 1.515 * Number(row.peptide.includes('W')) +
      1.290 * Number(row.peptide.split('Y').length - 1 >= 2) -
      1.052 * Number(sandbox.output.pi >= 10)};
});
fs.writeFileSync(process.argv[3], JSON.stringify(reference, null, 2) + '\n');
console.log(`Verified original app descriptors for ${reference.records.length} sequences`);
