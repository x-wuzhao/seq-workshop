/*
 * v3 生信增强的冒烟测试：把 index.html 里的 v3 代码块抽出来，在 Node 的 vm 里跑纯函数。
 * 覆盖：最近邻 Tm、引物二聚体 ΔG、线状/环状酶切、黏性末端、局部比对、GenBank 生成、线粒体密码表。
 */
const fs = require('node:fs');
const vm = require('node:vm');
const assert = require('node:assert/strict');

const html = fs.readFileSync('index.html', 'utf8').replace(/\r\n/g, '\n');
const start = html.indexOf('十三、v3 生信增强');
assert.ok(start > 0, 'v3 block should exist');
const endMark = 'var _v3Ready = true;';
const end = html.indexOf(endMark) + endMark.length;
assert.ok(end > start, 'v3 block should end with marker');
const block = html.slice(html.lastIndexOf('/*', start), end);

const prelude = `
var CODON=(function(){var c={};var aas='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
var bases='TCAG';for(var i=0;i<64;i++){var b1=bases[Math.floor(i/16)],b2=bases[Math.floor(i/4)%4],b3=bases[i%4];c[b1+b2+b3]=aas[i];}return c;})();
function reverseComplement(seq){var m={A:'T',T:'A',C:'G',G:'C',N:'N'},o='';for(var i=seq.length-1;i>=0;i--)o+=m[seq.charAt(i)];return o;}
var ENZYMES=[['EcoRI','GAATTC',1],['BamHI','GGATCC',1],['HindIII','AAGCTT',1],['SmaI','CCCGGG',3],['PstI','CTGCAG',5],['ScaI','AGTACT',3],['SalI','GTCGAC',1]];
var DEGEN={R:['A','G'],Y:['C','T'],S:['G','C'],W:['A','T'],K:['G','T'],M:['A','C'],B:['C','G','T'],D:['A','G','T'],H:['A','C','T'],V:['A','C','G'],N:['A','C','G','T']};
var DEG_COMP={R:'Y',Y:'R',S:'S',W:'W',K:'M',M:'K',B:'V',V:'B',D:'H',H:'D',N:'N',A:'T',T:'A',C:'G',G:'C'};
function baseMatch(p,ch){if(p===ch)return true;var d=DEGEN[p];if(d){for(var i=0;i<d.length;i++)if(d[i]===ch)return true;}return false;}
function siteAt(seq,pos,site){for(var i=0;i<site.length;i++){if(pos+i>=seq.length||!baseMatch(site.charAt(i),seq.charAt(pos+i)))return false;}return true;}
function rcPattern(p){var out='';for(var i=p.length-1;i>=0;i--)out+=(DEG_COMP[p.charAt(i)]||'N');return out;}
function scanEnzyme(seq,enz){var site=enz[1],L=site.length,c1=enz[2],c2=L-c1,rcs=rcPattern(site),hits=[];
for(var i=0;i+L<=seq.length;i++){if(siteAt(seq,i,site))hits.push({start:i,cut:i+c1});
if(rcs!==site&&siteAt(seq,i,rcs))hits.push({start:i,cut:i+c2,rc:true});}hits.sort(function(a,b){return a.cut-b.cut});return hits;}
var lastSeq='',lastEnz=null;
var document={getElementById:function(){return null;},querySelectorAll:function(){return [];},createElement:function(){return {};},documentElement:{getAttribute:function(){return null;}}};
var window={};
function analyze(){}
`;

const ctx = {SeqCore:require('../assets/js/sequence-core.js')};
vm.createContext(ctx);
vm.runInContext(prelude + block, ctx);

/* 1. 最近邻 Tm：20–22 mer 应落在 54–64 °C（50 mM Na+、0.25 µM 引物） */
const tms = ['ATGGCCATTGTAATGGGCCG', 'GTCGACCTGCAGCCCGGGAT', 'CGCGAATTCGCGATGGCCATTG'].map(s => ctx.tmNN(s));
tms.forEach((t, i) => {
  assert.ok(Number.isFinite(t), `Tm ${i} should be finite`);
  assert.ok(t > 52 && t < 66, `Tm ${i}=${t.toFixed(1)} out of range`);
});

/* 2. 完全互补的 20 mer 二聚体 ΔG 应明显为负；随机序列应接近 0 */
assert.ok(ctx.pairDG('ATGGCCATTGTAATGGGCCG', 'CGGCCCATTTACATGGCCAT') < -15, 'perfect duplex dG should be strongly negative');
assert.ok(ctx.selfDimerDG('ATGGCCATTGTAATGGGCCG') < 0, 'self dimer dG should be negative');

/* 3. 线状酶切：片段长度之和 = 序列长度 */
const seq = 'ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGCGTACGTAGCTAGCTAGCATCGATCGTAGCTAGCTAGCGGATCCGAATTCAAAAAAGCTTGGTACCGAGCTCGCGGCCGCTCTAGACTCGAGGTCGACCTGCAGCCCGGGATCGATGGTACCATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG';
const eco = ctx.scanEnzyme(seq, ['EcoRI', 'GAATTC', 1]);
const lin = ctx.digestFragments(seq.length, eco, false);
const sumLin = lin.frags.reduce((a, f) => a + f.len, 0);
assert.equal(sumLin, seq.length, 'linear fragments should sum to sequence length');
assert.equal(lin.frags.length, eco.length + 1, 'linear fragment count = cuts + 1');

/* 4. 环状酶切：片段数 = 切点数，长度之和仍等于序列长度 */
const circ = ctx.digestFragments(seq.length, eco, true);
assert.equal(circ.frags.length, eco.length, 'circular fragment count = cuts');
assert.equal(circ.frags.reduce((a, f) => a + f.len, 0), seq.length, 'circular fragments should sum to length');

/* 5. 黏性末端 */
assert.equal(ctx.enzymeOverhang('GAATTC', 1).type, '5′ 突出');
assert.equal(ctx.enzymeOverhang('GAATTC', 1).seq, 'AATT');
assert.equal(ctx.enzymeOverhang('CTGCAG', 5).type, '3′ 突出');
assert.equal(ctx.enzymeOverhang('CCCGGG', 3).type, '平端');

/* 6. 局部比对：两条含突变的同源序列应得到高一致性 */
const b1 = seq;
const b2 = seq.split('').map((c, i) => (i % 20 === 5) ? ({ A: 'C', C: 'G', G: 'T', T: 'A' })[c] : c).join('');
const sw = ctx.swAlign(b1, b2, 1, -1, -2);
assert.ok(sw.score > 0, 'local alignment score should be positive');
assert.ok(sw.match / sw.len > 0.9, 'local alignment identity should exceed 90%');

/* 7. GenBank 导出结构 */
const gb = ctx.buildGenBank(seq, { name: 'demo', circular: false });
assert.ok(gb.startsWith('LOCUS'), 'GenBank should start with LOCUS');
assert.ok(gb.trimEnd().endsWith('//'), 'GenBank should end with //');
assert.ok(/^FEATURES/m.test(gb), 'GenBank should contain FEATURES');
assert.ok(/^ORIGIN/m.test(gb), 'GenBank should contain ORIGIN');
assert.ok(gb.includes('misc_feature'), 'GenBank should contain feature annotations');

/* 8. 线粒体遗传密码表 */
ctx.activeCodonTable = 'vertebrate_mito';
assert.equal(ctx.translateFrame('TGA', 0), 'W');
ctx.activeCodonTable = 'standard';
assert.equal(ctx.translateFrame('TGA', 0), '*');

console.log('v3 bio smoke tests passed');
