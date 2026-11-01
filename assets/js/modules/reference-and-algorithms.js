/* =========================================================
   一、数据：遗传密码表 / 氨基酸名 / 示例序列
   ========================================================= */

// 标准遗传密码表（DNA 密码子 → 氨基酸单字母；* 为终止密码子）
var CODON = {
  TTT:'F',TTC:'F',TTA:'L',TTG:'L',
  TCT:'S',TCC:'S',TCA:'S',TCG:'S',
  TAT:'Y',TAC:'Y',TAA:'*',TAG:'*',
  TGT:'C',TGC:'C',TGA:'*',TGG:'W',
  CTT:'L',CTC:'L',CTA:'L',CTG:'L',
  CCT:'P',CCC:'P',CCA:'P',CCG:'P',
  CAT:'H',CAC:'H',CAA:'Q',CAG:'Q',
  CGT:'R',CGC:'R',CGA:'R',CGG:'R',
  ATT:'I',ATC:'I',ATA:'I',ATG:'M',
  ACT:'T',ACC:'T',ACA:'T',ACG:'T',
  AAT:'N',AAC:'N',AAA:'K',AAG:'K',
  AGT:'S',AGC:'S',AGA:'R',AGG:'R',
  GTT:'V',GTC:'V',GTA:'V',GTG:'V',
  GCT:'A',GCC:'A',GCA:'A',GCG:'A',
  GAT:'D',GAC:'D',GAA:'E',GAG:'E',
  GGT:'G',GGC:'G',GGA:'G',GGG:'G'
};

// 氨基酸单字母 → 中文名
var AA_NAME = {
  A:'丙氨酸',R:'精氨酸',N:'天冬酰胺',D:'天冬氨酸',C:'半胱氨酸',
  E:'谷氨酸',Q:'谷氨酰胺',G:'甘氨酸',H:'组氨酸',I:'异亮氨酸',
  L:'亮氨酸',K:'赖氨酸',M:'甲硫氨酸',F:'苯丙氨酸',P:'脯氨酸',
  S:'丝氨酸',T:'苏氨酸',W:'色氨酸',Y:'酪氨酸',V:'缬氨酸'
};
var AA_ORDER = ['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V'];

// 示例蛋白（刻意覆盖全部 20 种氨基酸），再反推出一段教学用"示例基因"
var SAMPLE_PROTEIN =
  'MKTAYIAKQRQISFVKNDCEGHLMPSTWY'+
  'VDERKHNQLGAIVWCFYSTPNDEMQRKLH'+
  'GAVIFWCYSTPDENQRKLMHGAIVWCFYS'+
  'TPDENQRKLMHGAIVWCFYSTPDENQRK'+
  'LMHGAIVWCFYSTP';
var AA_TO_CODON = {
  A:['GCA','GCC'],R:['CGT'],N:['AAC'],D:['GAT','GAC'],C:['TGC'],
  E:['GAA','GAG'],Q:['CAA'],G:['GGT','GGC'],H:['CAT'],I:['ATT'],
  L:['CTG','TTG'],K:['AAA'],M:['ATG'],F:['TTT'],P:['CCA','CCG'],
  S:['TCT','AGC'],T:['ACA','ACC'],W:['TGG'],Y:['TAT'],V:['GTT','GTG']
};
function buildSampleDNA(){
  var prot=SAMPLE_PROTEIN,codons=[];
  for(var i=0;i<prot.length;i++){
    var opts=AA_TO_CODON[prot.charAt(i)];
    codons.push(opts[i%opts.length]);
  }
  // 在不改变氨基酸的前提下，用相邻密码子组合嵌入常见酶切位点：
  // S(AGC)+F(TTT)=AGCT(AluI)；V(GTC)+D(GAC)=GTCGAC(SalI)；S(AGT)+T(ACT)=AGTACT(ScaI)
  for(var j=0;j+1<prot.length;j++){
    var a=prot.charAt(j),b=prot.charAt(j+1);
    if(a==='S'&&b==='F'){codons[j]='AGC';codons[j+1]='TTT';}
    else if(a==='V'&&b==='D'){codons[j]='GTC';codons[j+1]='GAC';}
    else if(a==='S'&&b==='T'){codons[j]='AGT';codons[j+1]='ACT';}
  }
  return codons.join('')+'TAA';
}
var SAMPLE_DNA = buildSampleDNA();

/* =========================================================
   二、生信核心算法（纯函数）
   ========================================================= */

// 解析输入：剥离 FASTA 标题行与空白，统计被忽略的非法字符
function parseInput(raw){
  var parsed=SeqCore.parseDNA(raw,{minLength:0});
  if(parsed.error)return {seq:'',removed:0,hadFasta:/>/.test(raw),error:parsed.error};
  return {seq:parsed.seq,removed:0,hadFasta:/>/.test(raw)};
}

function baseCounts(seq){
  var c={A:0,T:0,C:0,G:0,N:0};
  for(var i=0;i<seq.length;i++){
    var base=seq.charAt(i);
    if(c[base]!==undefined)c[base]++;
  }
  return c;
}

// 反向互补：互补配对后整体反转
function reverseComplement(seq){return SeqCore.reverseComplement(seq);}
function complement(seq){return SeqCore.complement(seq);}

// 转录：编码链 T → U
function transcribe(seq){return seq.replace(/T/g,'U');}

// 按指定读码框（0/1/2）翻译
function translateFrame(seq,frame){
  var p='';
  for(var i=frame;i+3<=seq.length;i+=3){
    p+=(CODON[seq.substr(i,3)]||'X');
  }
  return p;
}

// 在指定读码框内寻找全部 ORF（ATG → 终止密码子）
function findORFs(seq,frame){return SeqCore.findORFs(seq,frame,CODON,['ATG']);}

// GC 含量滑窗
function gcWindows(seq){
  var n=seq.length;
  // 短序列也保留至少一个窗口，避免 GC 图表出现空白。
  var win=n<20?n:(n<200?Math.max(20,Math.floor(n/10)):100);
  var step=Math.max(1,Math.floor(win/2));
  var pts=[],totalGC=0;
  for(var i=0;i+win<=n;i+=step){
    var gc=0;
    for(var j=i;j<i+win;j++){var ch=seq.charAt(j);if(ch==='G'||ch==='C')gc++;}
    pts.push({pos:i+Math.floor(win/2)+1,gc:+(gc/win*100).toFixed(1)});
    totalGC+=gc;
  }
  return {pts:pts,win:win,step:step};
}

/* =========================================================
   二（补）、限制性酶切 / 序列比对 / 物理参数
   ========================================================= */

// 限制酶：[名称, 识别序列(支持简并碱基), 编码链切点(在第 n 个碱基后切)]
var ENZYMES = [
  ['EcoRI','GAATTC',1],['BamHI','GGATCC',1],['HindIII','AAGCTT',1],
  ['XbaI','TCTAGA',1],['SalI','GTCGAC',1],['PstI','CTGCAG',5],
  ['SmaI','CCCGGG',3],['KpnI','GGTACC',5],['SacI','GAGCTC',5],
  ['NotI','GCGGCCGC',2],['XhoI','CTCGAG',1],['NcoI','CCATGG',1],
  ['NdeI','CATATG',2],['BglII','AGATCT',1],['SpeI','ACTAGT',1],
  ['ApaI','GGGCCC',5],['HaeIII','GGCC',2],['AluI','AGCT',2],
  ['TaqI','TCGA',1],['Sau3AI','GATC',1],['DraI','TTTAAA',3],
  ['EcoRV','GATATC',3],['AflII','CTTAAG',1],['AscI','GGCGCGCC',2],
  ['SphI','GCATGC',5],['NheI','GCTAGC',1],['StuI','AGGCCT',3],
  ['BspHI','TCATGA',1],['MspI','CCGG',1],['HpaI','GTTAAC',3],
  ['ScaI','AGTACT',3],['AatII','GACGTC',5]
];

// 简并碱基及其互补
var DEGEN = {R:['A','G'],Y:['C','T'],S:['G','C'],W:['A','T'],K:['G','T'],M:['A','C'],
  B:['C','G','T'],D:['A','G','T'],H:['A','C','T'],V:['A','C','G'],N:['A','C','G','T']};
var DEG_COMP = {R:'Y',Y:'R',S:'S',W:'W',K:'M',M:'K',B:'V',V:'B',D:'H',H:'D',N:'N',
  A:'T',T:'A',C:'G',G:'C'};

function baseMatch(p,ch){
  if(p===ch)return true;
  var d=DEGEN[p];
  if(d){for(var i=0;i<d.length;i++)if(d[i]===ch)return true;}
  return false;
}
function siteAt(seq,pos,site){
  for(var i=0;i<site.length;i++){
    if(pos+i>=seq.length||!baseMatch(site.charAt(i),seq.charAt(pos+i)))return false;
  }
  return true;
}
function rcPattern(p){
  var out='';
  for(var i=p.length-1;i>=0;i--)out+=(DEG_COMP[p.charAt(i)]||'N');
  return out;
}

// 扫描一种酶的全部识别位点与切点
function scanEnzyme(seq,enz){
  var site=enz[1],L=site.length,c1=enz[2],c2=L-c1,rcs=rcPattern(site),hits=[];
  for(var i=0;i+L<=seq.length;i++){
    if(siteAt(seq,i,site))hits.push({start:i,cut:i+c1});
    if(rcs!==site&&siteAt(seq,i,rcs))hits.push({start:i,cut:i+c2,rc:true});
  }
  hits.sort(function(a,b){return a.cut-b.cut});
  return hits;
}
function analyzeEnzymes(seq){
  var res=[];
  for(var e=0;e<ENZYMES.length;e++){
    res.push({enz:ENZYMES[e],hits:scanEnzyme(seq,ENZYMES[e])});
  }
  res.sort(function(a,b){return b.hits.length-a.hits.length;});
  return res;
}

// Needleman–Wunsch 全局比对
function nwAlign(a,b,sc,mm,gap){
  var n=a.length,m=b.length;
  if(!n||!m)return null;
  var F=[],T=[];
  for(var i=0;i<=n;i++){F.push(new Array(m+1));T.push(new Array(m+1));}
  F[0][0]=0;
  for(var x=1;x<=n;x++){F[x][0]=x*gap;T[x][0]=1;}
  for(var y=1;y<=m;y++){F[0][y]=y*gap;T[0][y]=2;}
  for(var ii=1;ii<=n;ii++){
    for(var jj=1;jj<=m;jj++){
      var s=(a.charAt(ii-1)===b.charAt(jj-1))?sc:mm;
      var d=F[ii-1][jj-1]+s,u=F[ii-1][jj]+gap,l=F[ii][jj-1]+gap;
      var best=d,bt=0;
      if(u>best){best=u;bt=1;}
      if(l>best){best=l;bt=2;}
      F[ii][jj]=best;T[ii][jj]=bt;
    }
  }
  var aa='',bb='',mid='',ci=n,cj=m;
  while(ci>0||cj>0){
    var t=T[ci][cj];
    if(t===0){
      var ca=a.charAt(ci-1),cb=b.charAt(cj-1);
      aa=ca+aa;bb=cb+bb;mid=(ca===cb?'|':'·')+mid;ci--;cj--;
    }else if(t===1){
      aa=a.charAt(ci-1)+aa;bb='-'+bb;mid=' '+mid;ci--;
    }else{
      aa='-'+aa;bb=b.charAt(cj-1)+bb;mid=' '+mid;cj--;
    }
  }
  var matchN=0;
  for(var z=0;z<mid.length;z++)if(mid.charAt(z)==='|')matchN++;
  return {a:aa,b:bb,mid:mid,score:F[n][m],match:matchN,len:mid.length};
}

// 解链温度（短序列 Wallace 规则；长序列经验公式）
function calcTm(c,n){
  if(n<14)return 4*(c.G+c.C)+2*(c.A+c.T);
  return +(64.9+41*((c.G+c.C-16.4)/n)).toFixed(1);
}
// 单链分子量（Da，近似）
function calcMW(c){
  return c.A*313.2+c.T*304.2+c.G*329.2+c.C*289.2-61;
}

/* ===== 蛋白质理化性质（参照 ExPASy ProtParam）===== */
// 氨基酸残基平均分子量（Da）
var AA_RES_MW = {
  A:71.0788,R:156.1875,N:114.1038,D:115.0886,C:103.1388,
  E:129.1155,Q:128.1307,G:57.0519,H:137.1411,I:113.1594,
  L:113.1594,K:128.1741,M:131.1926,F:147.1766,P:97.1167,
  S:87.0782,T:101.1051,W:186.2132,Y:163.1760,V:99.1326
};
// Kyte–Doolittle 疏水指数
var AA_HYDRO = {
  A:1.8,R:-4.5,N:-3.5,D:-3.5,C:2.5,E:-3.5,Q:-3.5,G:-0.4,
  H:-3.2,I:4.5,L:3.8,K:-3.9,M:1.9,F:2.8,P:-1.6,S:-0.8,
  T:-0.7,W:-0.9,Y:-1.3,V:4.2
};
// 可电离基团 pKa（Bjellqvist / ProtParam 口径）
var PKA_POS = {Nterm:7.7,K:10.0,R:12.0,H:6.0};
var PKA_NEG = {Cterm:3.55,D:4.05,E:4.45,C:8.3,Y:10.0};

// 给定 pH 下蛋白质净电荷
function protCharge(pH, n){
  var pos = 1/(1+Math.pow(10,pH-PKA_POS.Nterm))
    + (n.K||0)/(1+Math.pow(10,pH-PKA_POS.K))
    + (n.R||0)/(1+Math.pow(10,pH-PKA_POS.R))
    + (n.H||0)/(1+Math.pow(10,pH-PKA_POS.H));
  var neg = 1/(1+Math.pow(10,PKA_NEG.Cterm-pH))
    + (n.D||0)/(1+Math.pow(10,PKA_NEG.D-pH))
    + (n.E||0)/(1+Math.pow(10,PKA_NEG.E-pH))
    + (n.C||0)/(1+Math.pow(10,PKA_NEG.C-pH))
    + (n.Y||0)/(1+Math.pow(10,PKA_NEG.Y-pH));
  return pos - neg;
}
// 二分法求等电点
function calcPI(n){
  var lo=0,hi=14;
  for(var i=0;i<60;i++){
    var mid=(lo+hi)/2;
    if(protCharge(mid,n)>0)lo=mid;else hi=mid;
  }
  return (lo+hi)/2;
}
// 蛋白质理化性质主函数，prot 为氨基酸序列（不含终止符）
function calcProtProps(prot){
  if(/X/.test(prot))return {pI:null,mw:0,gravy:null,aromatic:null,len:0};
  var n={},len=0,mw=0,hydro=0;
  for(var i=0;i<prot.length;i++){
    var a=prot.charAt(i);
    if(AA_RES_MW[a]===undefined)continue;  // 跳过终止符等
    n[a]=(n[a]||0)+1;len++;
    mw+=AA_RES_MW[a];
    hydro+=AA_HYDRO[a];
  }
  if(len===0)return {pI:null,mw:0,gravy:null,aromatic:null,len:0};
  mw+=18.015;  // 末端水
  var pI=calcPI(n);
  var aromatic=((n.F||0)+(n.Y||0)+(n.W||0))/len;
  return {pI:pI,mw:mw,gravy:hydro/len,aromatic:aromatic,len:len};
}

/* =========================================================
   三、渲染与交互
   ========================================================= */

// HERO 碱基滚动轨道
(function(){
  var base='ATCG',txt='';
  for(var i=0;i<500;i++){
    var ch=base.charAt(Math.floor(Math.random()*4));
    txt+='<span class="b'+ch+'">'+ch+'</span>';
  }
  var lane=document.getElementById('heroLane');
  lane.innerHTML=txt+txt;
})();
