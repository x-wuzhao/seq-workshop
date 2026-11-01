/* Sequence Workshop Next: deterministic, browser/Node shared calculations. */
(function(root, factory) {
  var api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.SeqCore = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function() {
  'use strict';
  var pairs = {A:'T',T:'A',C:'G',G:'C',N:'N'};
  function parseDNA(raw, options) {
    options = options || {};
    var lines = String(raw || '').trim().split(/\r?\n/), headers = 0, header = '', seq = '';
    for (var i=0; i<lines.length; i++) {
      var line = lines[i].trim();
      if (!line) continue;
      if (line[0] === '>') {
        headers++;
        if (headers > 1 || seq.length) return {error:'一次只能分析一条 FASTA 记录，请分别选择或粘贴；不会自动拼接多条序列。'};
        header = line.slice(1).trim();
      } else seq += line.replace(/\s/g, '');
    }
    seq = seq.toUpperCase();
    var invalid = seq.match(/[^ATCGUN]/g);
    if (invalid) return {error:'序列包含不支持的字符：' + Array.from(new Set(invalid)).slice(0,16).join('、') + '。支持 A/T/C/G/N，U 自动转为 T。'};
    seq = seq.replace(/U/g,'T');
    if (seq.length < (options.minLength === undefined ? 3 : options.minLength)) return {error:'请输入至少 3 个碱基的 DNA 序列。'};
    if (seq.length > (options.maxLength || 50000)) return {error:'本工具单条序列上限为 50,000 bp，请先截取目标区域。'};
    var warnings=[];
    if (seq.includes('N')) warnings.push('含 N（未知碱基）：单列未知比例，翻译显示 X；Tm、分子量和部分评分不予估算。');
    if (seq.length > 10000) warnings.push('序列较长，建议按目标区域开展后续分析。');
    return {seq:seq,header:header,warning:warnings.join('；') || null};
  }
  function complement(seq) { return seq.replace(/[ATCGN]/g, function(base) { return pairs[base]; }); }
  function reverseComplement(seq) { return complement(seq).split('').reverse().join(''); }
  function composition(seq) {
    var counts={A:0,T:0,C:0,G:0,N:0};
    for (var i=0;i<seq.length;i++) if (counts[seq[i]] !== undefined) counts[seq[i]]++;
    var divisor=seq.length || 1;
    return {counts:counts,gcPct:(counts.G+counts.C)/divisor*100,atPct:(counts.A+counts.T)/divisor*100,nPct:counts.N/divisor*100};
  }
  function translate(seq, frame, map) {
    var protein='';
    for(var i=frame;i+3<=seq.length;i+=3) protein+=map[seq.slice(i,i+3)] || 'X';
    return protein;
  }
  function findORFs(seq, frame, map, starts) {
    var protein=translate(seq,frame,map), nextStop=-1, result=[];
    for(var i=protein.length-1;i>=0;i--) {
      if(protein[i]==='*') { nextStop=i; continue; }
      var start=frame+i*3;
      if(nextStop>=0 && starts.indexOf(seq.slice(start,start+3))>=0) {
        var end=frame+(nextStop+1)*3;
        // Alternative initiation codons code for methionine at the CDS start.
        result.push({start:start,end:end,len:nextStop-i+1,prot:'M'+protein.slice(i+1,nextStop),stopCodon:seq.slice(end-3,end)});
      }
    }
    return result.reverse();
  }
  function translateCDS(seq,start,map) {
    var protein='',stop=null,end=start;
    for(var i=start;i+3<=seq.length;i+=3) {
      var codon=seq.slice(i,i+3),aa=map[codon] || 'X'; end=i+3;
      if(aa==='*') { stop=codon; break; }
      protein+=(i===start?'M':aa);
    }
    return {protein:protein,stopCodon:stop,end:end};
  }
  function mutationEffect(original,mutant,type,position,frameSpec,map,starts) {
    frameSpec=frameSpec || '+1';
    var reverse=frameSpec[0]==='-',frame=Number(frameSpec.slice(1))-1;
    var before=reverse?reverseComplement(original):original;
    var after=reverse?reverseComplement(mutant):mutant;
    var positionOnStrand=reverse?original.length-1-position:position;
    if(type==='insert' && reverse) positionOnStrand++;
    var orfs=findORFs(before,frame,map,starts).filter(function(orf) { return positionOnStrand>=orf.start && positionOnStrand<orf.end; });
    orfs.sort(function(a,b) { return b.len-a.len || a.start-b.start; });
    var ref=orfs[0];
    var result={frame:frameSpec,reference:ref || null,originalProtein:ref?ref.prot:'',mutantProtein:''};
    function done(kind,label,explanation) {
      result.type=kind;result.label=label;result.explain=explanation;return result;
    }
    if(!ref) return done('unresolved','无法判定编码影响','该位置不在所选读码框的完整 ORF 中。请确认链方向、读码框及编码区，不能据此判定同义或错义突变。');
    var translated=translateCDS(after,ref.start,map);result.mutantProtein=translated.protein;
    if(type!=='point') return done('frameshift','移码突变','所选编码区内插入或删除 1 个碱基，导致后续阅读框移动；前后蛋白沿同一原始 ORF 起点比较。');
    var codonStart=ref.start+Math.floor((positionOnStrand-ref.start)/3)*3;
    var oldCodon=before.slice(codonStart,codonStart+3),newCodon=after.slice(codonStart,codonStart+3);
    result.oldCodon=oldCodon;result.newCodon=newCodon;
    var oldAA=map[oldCodon],newAA=map[newCodon];
    if(!oldAA || !newAA) return done('unresolved','含未知碱基，无法判定','受影响密码子含未知碱基，不输出确定的突变类别。');
    if(codonStart===ref.start && starts.indexOf(newCodon)<0) return done('start-loss','起始密码子改变','原始 ORF 的起始密码子被改变；实际翻译起始需要进一步验证，不能直接当作普通错义突变。');
    if(newAA==='*' && oldAA!=='*') return done('nonsense','无义突变','在原始编码区内产生提前终止密码子。此处比较同一 ORF 的截短产物，不重新挑选其他最长 ORF。');
    if(oldAA==='*' && newAA!=='*') return done('stop-loss','终止密码子丢失','原终止密码子改变，翻译可能延长至后续终止位点；末端序列不足时结果是不完整的。');
    if(oldAA===newAA) return done('synonym','同义突变','所选 ORF 中该密码子编码的氨基酸保持一致。');
    return done('missense','错义突变','所选 ORF 中该密码子编码的氨基酸发生改变。');
  }
  function preferenceScore(seq,bias,map) {
    var sum=0,count=0;
    for(var i=0;i+3<=seq.length;i+=3) {
      var codon=seq.slice(i,i+3),aa=map[codon];
      if(!aa) return null;
      if(aa==='*') continue;
      var rank=(bias[aa] || []).indexOf(codon);
      if(rank<0) return null;
      sum+=Math.log(1/(rank+1));count++;
    }
    return count?Math.exp(sum/count):null;
  }
  function escapeHTML(value) {
    return String(value).replace(/[&<>"']/g,function(c) { return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]; });
  }
  return {version:'4.0.0',parseDNA:parseDNA,complement:complement,reverseComplement:reverseComplement,composition:composition,translate:translate,findORFs:findORFs,mutationEffect:mutationEffect,preferenceScore:preferenceScore,escapeHTML:escapeHTML};
});
