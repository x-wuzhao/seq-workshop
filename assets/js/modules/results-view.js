var charts={gc:null,aa:null};
var lastChartArgs=null;
function isDarkTheme(){return document.documentElement.getAttribute('data-theme')==='dark';}
function chartAxisColor(){return isDarkTheme()?'#9FB6B2':'#3E5C59';}
function chartSplitColor(){return isDarkTheme()?'rgba(167,224,222,.13)':'rgba(13,43,42,.09)';}

function setMsg(cls,html){
  var m=document.getElementById('inputMsg');
  m.className='msg '+cls;m.textContent=html;
}
function clearMsg(){var m=document.getElementById('inputMsg');m.className='msg';m.innerHTML='';}

// 碱基序列 → 彩色 HTML（每 10 个碱基留间隔，便于数位）
function seqHTML(seq,isRNA){
  var h='';
  for(var i=0;i<seq.length;i++){
    if(i>0&&i%10===0)h+='<span class="sep"></span>';
    var ch=seq.charAt(i);
    var cls=isRNA?(ch==='U'?'bU':'b'+ch):'b'+ch;
    h+='<span class="'+cls+'" data-seq-index="'+i+'">'+ch+'</span>';
  }
  return h;
}

function aaChar(ch){
  if(ch==='*')return '<span class="aaStop">*</span>';
  if(ch==='M')return '<span class="aaM">M</span>';
  return ch;
}
function protPart(prot,a,b){
  var s='';
  for(var i=a;i<b;i++)s+=aaChar(prot.charAt(i));
  return s;
}

// 构建单个读码框数据
function buildFrame(s,sign,f){
  var prot=translateFrame(s,f);
  var orfs=findORFs(s,f),best=null;
  for(var i=0;i<orfs.length;i++){if(!best||orfs[i].len>best.len)best=orfs[i];}
  var html='';
  if(best){
    var si=(best.start-f)/3,ei=si+best.len;
    html+=protPart(prot,0,si);
    html+='<span class="orf-hi">'+protPart(prot,si,ei)+'</span>';
    html+=protPart(prot,ei,prot.length);
  }else{
    html+=protPart(prot,0,prot.length);
  }
  return {label:sign+(f+1),prot:prot,best:best,html:html,codons:prot.length};
}

// 长序列折叠
function setupCollapse(lineEl,btnEl){
  if(lineEl.scrollHeight>150){
    lineEl.classList.add('collapsed');btnEl.hidden=false;
  }else{
    lineEl.classList.remove('collapsed');btnEl.hidden=true;
  }
}

function renderResults(seq){
  var rev=reverseComplement(seq),mrna=transcribe(seq);
  var counts=baseCounts(seq);
  var composition=SeqCore.composition(seq),gcPct=composition.gcPct,atPct=composition.atPct;

  // 统计卡
  document.getElementById('stLen').innerHTML=seq.length+'<small> bp</small>';
  document.getElementById('stGC').innerHTML=gcPct.toFixed(1)+'<small> %</small>';
  var stATel=document.getElementById('stAT');
  if(stATel)stATel.innerHTML=atPct.toFixed(1)+'<small> %</small>';
  document.getElementById('stA').textContent=counts.A;
  document.getElementById('stT').textContent=counts.T;
  document.getElementById('stG').textContent=counts.G;
  document.getElementById('stC').textContent=counts.C;
  // v3.1 碱基组成条
  var totalBases=seq.length||1;
  [['cbA',counts.A],['cbT',counts.T],['cbG',counts.G],['cbC',counts.C]].forEach(function(pair){
    var el=document.getElementById(pair[0]);
    if(el)el.style.width=(pair[1]/totalBases*100).toFixed(2)+'%';
  });
  var compLabel=document.getElementById('smCompLabel');
  if(compLabel)compLabel.textContent='GC '+gcPct.toFixed(1)+'% · AT '+atPct.toFixed(1)+'% · 未知 N '+composition.nPct.toFixed(1)+'%（均按总长度计算）';
  var tm=calcTm(counts,seq.length),mwDa=calcMW(counts);
  document.getElementById('stTm').innerHTML=counts.N?'未知（含 N）':Number(tm).toFixed(1)+'<small> °C</small>';
  document.getElementById('stMW').innerHTML=counts.N?'未知（含 N）':(mwDa/1000).toFixed(1)+'<small> kDa</small>';

  // 序列展示
  document.getElementById('lenInput').textContent='共 '+seq.length+' bp';
  var line1=document.getElementById('seqInputShow');
  var line2=document.getElementById('seqRevShow');
  var line3=document.getElementById('seqMrnaShow');
  line1.innerHTML=seqHTML(seq,false);
  line2.innerHTML=seqHTML(rev,false);
  line3.innerHTML=seqHTML(mrna,true);
  setupCollapse(line1,line1.nextElementSibling);
  setupCollapse(line2,line2.nextElementSibling);
  setupCollapse(line3,line3.nextElementSibling);

  // 六框
  var defs=[
    {s:seq,sign:'+',f:0},{s:seq,sign:'+',f:1},{s:seq,sign:'+',f:2},
    {s:rev,sign:'-',f:0},{s:rev,sign:'-',f:1},{s:rev,sign:'-',f:2}
  ];
  var frames=[];
  for(var i=0;i<defs.length;i++){frames.push(buildFrame(defs[i].s,defs[i].sign,defs[i].f));}
  var gBest=null,gIdx=-1;
  for(var j=0;j<frames.length;j++){
    if(frames[j].best&&(!gBest||frames[j].best.len>gBest.len)){
      gBest=frames[j].best;
      gBest.frameSeq=defs[j].s;
      gBest.frameSign=defs[j].sign;
      gBest.frame=defs[j].f;
      gIdx=j;
    }
  }
  var box=document.getElementById('framesBox');box.innerHTML='';
  for(var k=0;k<frames.length;k++){
    var fr=frames[k];
    var row=document.createElement('div');
    row.className='frame-row'+(k===gIdx?' best':'');
    var head='<div class="frame-head"><span class="frame-name">'+fr.label+'</span>';
    head+='<span class="frame-meta">本框共翻译 '+fr.codons+' 个密码子</span>';
    if(fr.best){
      head+='<span class="orf-badge">'+(k===gIdx?'最长 ORF · ':'ORF · ')+(fr.best.len-1)+' 个氨基酸</span>';
    }else{
      head+='<span class="frame-meta">未检测到完整 ORF（ATG→终止）</span>';
    }
    head+='</div>';
    row.innerHTML=head+'<div class="prot-line">'+fr.html+'</div>';
    box.appendChild(row);
  }

  // 氨基酸组成：优先最长 ORF，否则取最长框翻译
  var protForAA;
  if(gBest){protForAA=gBest.prot;}
  else{
    protForAA=frames[0].prot;
    for(var m=1;m<frames.length;m++){if(frames[m].prot.length>protForAA.length)protForAA=frames[m].prot;}
  }
  var aaCounts={};
  for(var z=0;z<AA_ORDER.length;z++){aaCounts[AA_ORDER[z]]=0;}
  for(var q=0;q<protForAA.length;q++){
    var c2=protForAA.charAt(q);
    if(aaCounts[c2]!==undefined)aaCounts[c2]++;
  }
  document.getElementById('aaHint').textContent=gBest?
    '基于最长 ORF（'+gBest.prot.length+' 个氨基酸）翻译产物统计。':
    '未发现完整 ORF，基于最长读码框的翻译产物统计。';

  // 蛋白质理化性质
  var props=calcProtProps(protForAA);
  document.getElementById('ppPI').textContent=props.pI===null?'—':props.pI.toFixed(2);
  document.getElementById('ppMW').innerHTML=props.len?((props.mw/1000).toFixed(2)+'<small> kDa</small>'):'—';
  document.getElementById('ppGRAVY').textContent=props.gravy===null?'—':props.gravy.toFixed(3);
  document.getElementById('ppArom').innerHTML=props.aromatic===null?'—':(props.aromatic*100).toFixed(1)+'<small> %</small>';

  // 酶切位点 + 密码子频率
  lastSeq=seq;lastEnz=analyzeEnzymes(seq);
  renderEnzMap(seq,lastEnz);
  renderEnzTable(seq,lastEnz,document.getElementById('enzFilter').value);
  setupGel(seq,lastEnz);
  renderCodonGrid(seq,gBest);

  document.getElementById('results').hidden=false;
  initSeqMap();
  renderSeqMap(seq);
  renderMapFeatures();
  if (typeof scheduleWorkspaceSave === 'function') scheduleWorkspaceSave();
  if(window.updateSectionSpy)window.updateSectionSpy();
  renderStructure3D();
  renderCharts(seq,gcPct,aaCounts);
}

var lastSeq='',lastEnz=null;
