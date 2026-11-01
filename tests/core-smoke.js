const fs = require('node:fs');
const vm = require('node:vm');
const assert = require('node:assert/strict');

const html = fs.readFileSync('index.html', 'utf8').replace(/\r\n/g, '\n');
const scripts = [...html.matchAll(/<script(?:\s[^>]*)?>([\s\S]*?)<\/script>/g)].map(m => m[1]);
const coreMarker = '/* =========================================================\n   三、渲染与交互';
const core = scripts.find(script => script.includes(coreMarker));
assert.ok(core, 'core script should be present');

const context = {SeqCore:require('../assets/js/sequence-core.js')};
vm.runInNewContext(core.slice(0, core.indexOf(coreMarker)), context);

assert.equal(JSON.stringify(context.parseInput('>demo\nATN U')), JSON.stringify({seq: 'ATNT', removed: 0, hadFasta: true}));
assert.equal(context.reverseComplement('ATNG'), 'CNAT');
assert.equal(JSON.stringify(context.baseCounts('ATNG')), JSON.stringify({A: 1, T: 1, C: 0, G: 1, N: 1}));
assert.equal(context.translateFrame('ATGTAA', 0), 'M*');
assert.equal(JSON.stringify(context.findORFs('ATGTAA', 0)[0]), JSON.stringify({start: 0, end: 6, len: 2, prot: 'M', stopCodon: 'TAA'}));
assert.equal(JSON.stringify(context.findORFs('ATGATGTAA', 0).map(orf => orf.start)), JSON.stringify([0, 3]));
assert.equal(context.findORFs('ATGTAA', 0)[0].stopCodon, 'TAA');
assert.equal(context.gcWindows('ATG').pts.length, 1);

const emptyProps = context.calcProtProps('');
assert.equal(emptyProps.pI, null);
assert.equal(emptyProps.gravy, null);
assert.equal(emptyProps.aromatic, null);

console.log('core smoke tests passed');
