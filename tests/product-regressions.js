const assert=require('node:assert/strict');
const core=require('../assets/js/sequence-core');
const format=require('../assets/js/workspace-format');
const app=require('./load-app')();
const fixture=require('./fixtures/biopython-1.86.json');
let checks=0;
function test(name,fn){fn();checks++;console.log('PASS',name);}
test('strict single-record FASTA across entry points',()=>{
  const raw='>one\nATGTAA\n>two\nATGTGA';
  assert.ok(core.parseDNA(raw).error);assert.ok(app.enhancedParseInput(raw).error);assert.ok(app.parseInput(raw).error);
  assert.ok(core.parseDNA('>header-only').error);assert.ok(core.parseDNA('ATGXYZ').error);
  assert.equal(core.parseDNA('>one\naug nntaa').seq,'ATGNNTAA');
});
test('unknown base composition has an explicit denominator',()=>{
  assert.equal(core.reverseComplement('ATGNNTAA'),'TTANNCAT');
  const r=core.composition('ATGNNTAA');assert.equal(r.atPct,62.5);assert.equal(r.gcPct,12.5);assert.equal(r.nPct,25);
  assert.equal(r.atPct+r.gcPct+r.nPct,100);
});
test('all 128 standard and mitochondrial codons match Biopython 1.86',()=>{
  for(const [id,name] of [['1','standard'],['2','vertebrate_mito']]){
    app.activeCodonTable=name;
    for(const [codon,expected] of Object.entries(fixture.tables[id].translation))assert.equal(app.translateFrame(codon,0),expected,name+' '+codon);
  }
  app.activeCodonTable='standard';
});
test('nearest-neighbor Tm agrees with independent reference within 0.01 C',()=>{
  for(const item of fixture.meltingTemperatures)assert.ok(Math.abs(app.tmNN(item.sequence)-item.temperature)<0.01,item.sequence+': '+app.tmNN(item.sequence)+' vs '+item.temperature);
  assert.ok(Number.isNaN(app.tmNN('ATGNNNTAA')));
});
test('mutation uses the same original ORF on both strands',()=>{
  const original=app.SAMPLE_DNA;
  const changed=app.applyMutation(original,'point',87,'A');
  const result=app.classifyMutation(original,changed.seq,'point',86,'+1');
  assert.equal(result.oldCodon,'TAT');assert.equal(result.newCodon,'TAA');assert.equal(result.type,'nonsense');
  assert.equal(result.originalProtein.length,129);assert.equal(result.mutantProtein.length,28);
  const reverse=core.reverseComplement(original),mutantReverse=core.reverseComplement(changed.seq);
  assert.equal(app.classifyMutation(reverse,mutantReverse,'point',original.length-1-86,'-1').type,'nonsense');
});
test('start loss, stop loss, synonymous and uncertain changes are distinguished',()=>{
  assert.equal(app.classifyMutation('ATGAAATAA','ACGAAATAA','point',1,'+1').type,'start-loss');
  assert.equal(app.classifyMutation('ATGAAATAA','ATGAAACAA','point',6,'+1').type,'stop-loss');
  assert.equal(app.classifyMutation('ATGAAATAA','ATGAAGTAA','point',5,'+1').type,'synonym');
  assert.equal(app.classifyMutation('CCCTAA','CCATAA','point',2,'+1').type,'unresolved');
  assert.ok(app.applyMutation('ATGTAA','point','','A').error);
});
test('preference rank is not mislabeled as CAI and counts actual substitutions',()=>{
  const r=app.optimizeCodons('ATGTTATAA','ecoli');assert.equal(r.seq,'ATGCTGTAA');assert.equal(r.changes,1);
  assert.ok(Math.abs(core.preferenceScore('ATGTTATAA',app.SPECIES_CODON_BIAS.ecoli,app.CODON)-Math.sqrt(0.5))<1e-12);
  assert.equal(core.preferenceScore('ATGNNNTAA',app.SPECIES_CODON_BIAS.ecoli,app.CODON),null);
  assert.equal(app.optimizeCodons('ATGGCTTAG','ecoli').seq.slice(-3),'TAG');
  app.activeCodonTable='vertebrate_mito';
  const mito='ATGATATGATAA',mitoResult=app.optimizeCodons(mito,'ecoli');
  assert.equal(app.translateFrame(mitoResult.seq,0),app.translateFrame(mito,0));
  app.activeCodonTable='standard';
});
test('portable backup round trip validates all data before import',()=>{
  const backup=format.create('ATGTAA',[{start:1,end:6,strand:'+',label:'<img src=x onerror=alert(1)>',note:'hello'}],{codonTable:'standard',topology:'circular'});
  assert.deepEqual(format.validate(JSON.parse(JSON.stringify(backup))),backup);
  assert.throws(()=>format.validate({...backup,schemaVersion:2}));
  assert.throws(()=>format.validate({...backup,features:[{start:1,end:99,strand:'+',label:'bad'}]}));
  app.mapFeatures=backup.features;app.renderMapFeatures();
  assert.ok(!app.document.getElementById('mapFeatures').innerHTML.includes('<img'));
  assert.ok(app.document.getElementById('mapFeatures').innerHTML.includes('&lt;img'));
});
test('new analysis invalidates every derived tool result',()=>{
  app.window._lastPrimers=[{old:true}];app.window._seqMapMutation={old:true};app.window._optimizedSeq='old';
  for(const id of ['mutResult','primerResult','alignResult','optResult'])app.document.getElementById(id).innerHTML='old';
  app.resetSequenceDependentPanels();
  assert.equal(app.window._lastPrimers.length,0);assert.equal(app.window._seqMapMutation,null);assert.equal(app.window._optimizedSeq,'');
  for(const id of ['mutResult','primerResult','alignResult','optResult'])assert.equal(app.document.getElementById(id).innerHTML,'');
});
test('long ORF-free sequence and nested ORFs preserve exact results',()=>{
  assert.equal(app.findORFs('ATG'.repeat(16000),0).length,0);
  assert.deepEqual(Array.from(app.findORFs('ATGATGTAA',0),o=>o.start),[0,3]);
});
test('GenBank complex locations are skipped, multiple records are rejected',()=>{
  assert.equal(app.parseGbLocation('join(1..3,7..9)'),null);
  assert.ok(app.parseGenBank('LOCUS one\nORIGIN\n 1 atgtaa\n//\nLOCUS two\nORIGIN\n 1 atgtaa\n//').error);
});
console.log(checks+' product regression groups passed; 128 codons and '+fixture.meltingTemperatures.length+' independent Tm cases.');
