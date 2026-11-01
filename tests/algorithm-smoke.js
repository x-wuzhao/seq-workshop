const fs = require('node:fs');
const vm = require('node:vm');
const assert = require('node:assert/strict');

const html = fs.readFileSync('index.html', 'utf8');
const scripts = [...html.matchAll(/<script(?:\s[^>]*)?>([\s\S]*?)<\/script>/g)].map(m => m[1]);

function makeElement() {
  return {
    value: '', innerHTML: '', textContent: '', hidden: false, disabled: false,
    parentNode: {replaceChild() {}},
    style: {}, className: '', children: [],
    classList: {add() {}, remove() {}, toggle() {}, contains() { return false; }},
    setAttribute() {}, getAttribute() { return null; }, addEventListener() {},
    appendChild(child) { this.children.push(child); }, removeChild() {},
    querySelector() { return makeElement(); }, querySelectorAll() { return []; }, cloneNode() { return makeElement(); },
    scrollIntoView() {}, closest() { return null; }, focus() {}, click() {},
    getBoundingClientRect() { return {top: 0, bottom: 1}; }
  };
}

const elements = {};
const document = {
  documentElement: makeElement(), body: makeElement(),
  getElementById(id) { return elements[id] || (elements[id] = makeElement()); },
  querySelector() { return makeElement(); }, querySelectorAll() { return []; },
  createElement() { return makeElement(); }, addEventListener() {}
};
const window = {localStorage: {getItem() { return null; }, setItem() {}}, addEventListener() {}, scrollY: 0, scrollTo() {}};
const context = {
  SeqCore:require('../assets/js/sequence-core.js'), document, window, console, navigator: {serviceWorker: null},
  location: {protocol: 'http:', hostname: 'localhost'},
  setTimeout() { return 0; }, clearTimeout() {}, setInterval() { return 0; }, clearInterval() {},
  requestAnimationFrame() { return 0; }, cancelAnimationFrame() {},
  Blob: function Blob() {}, URL: {createObjectURL() { return ''; }}, alert() {},
  matchMedia() { return {matches: false}; }
};
context.globalThis = context;
vm.createContext(context);
scripts.forEach(script => vm.runInContext(script, context));

const optimized = context.optimizeCodons('ATGGCTTAA', 'ecoli');
assert.equal(optimized.changes, 1);
assert.equal(optimized.seq, 'ATGGCGTAA');
assert.doesNotThrow(() => context.optimizeCodons('ATGNNNTAA', 'ecoli'));
assert.equal(context.optimizeCodons('ATGNNNTAA', 'ecoli').seq, 'ATGNNNTAA');

const mutation = context.classifyMutation('ATGAAATAA', 'ATGTAATAA', 'point', 3);
assert.equal(mutation.type, 'nonsense');

console.log('algorithm smoke tests passed');
