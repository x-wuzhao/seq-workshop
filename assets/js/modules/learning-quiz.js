/* ============================================================
   高中题库系统
   ============================================================ */
var QUIZ_DATA = [
  {cat:'DNA结构', q:'下列关于DNA分子结构的叙述，正确的是？', options:['DNA分子由两条同向平行的脱氧核苷酸链组成','DNA分子中的脱氧核糖和磷酸交替连接，排列在外侧','碱基对通过肽键连接','A与T之间形成3个氢键，G与C之间形成2个氢键'], answer:1, explain:'DNA分子由两条反向平行的脱氧核苷酸链组成；脱氧核糖和磷酸交替连接排列在外侧，构成基本骨架；碱基对通过氢键连接；A与T之间形成2个氢键，G与C之间形成3个氢键。'},
  {cat:'DNA复制', q:'DNA分子复制时，解旋酶作用的部位是？', options:['磷酸二酯键','氢键','肽键','高能磷酸键'], answer:1, explain:'解旋酶的作用是断裂碱基对之间的氢键，使DNA双链解开。磷酸二酯键是连接脱氧核苷酸的化学键，由DNA聚合酶催化形成；肽键连接氨基酸；高能磷酸键存在于ATP中。'},
  {cat:'DNA复制', q:'一个DNA分子经过3次复制后，含有原始母链的DNA分子占总数的？', options:['1/2','1/4','1/8','1/16'], answer:1, explain:'DNA复制为半保留复制。1个DNA复制3次后产生2³=8个DNA分子，其中只有2个DNA分子含有原始母链（每个原始母链进入一个子代DNA），所以比例为2/8=1/4。'},
  {cat:'转录', q:'下列关于转录的叙述，错误的是？', options:['转录的场所主要是细胞核','转录需要RNA聚合酶催化','转录的原料是4种核糖核苷酸','转录时DNA两条链同时作为模板'], answer:3, explain:'转录时只有DNA的一条链作为模板（模板链），另一条链为编码链。转录的场所主要是细胞核（线粒体、叶绿体中也可发生）；需要RNA聚合酶催化；原料是4种核糖核苷酸（ATP、UTP、GTP、CTP）。'},
  {cat:'转录', q:'某DNA分子的模板链碱基序列为5′-ATGCAT-3′，其转录出的mRNA碱基序列为？', options:['5′-AUGCAU-3′','5′-UACGUA-3′','5′-ATGCAT-3′','5′-TACGTA-3′'], answer:0, explain:'转录时mRNA与模板链互补，且方向相反。模板链5′-ATGCAT-3′，互补链为3′-TACGTA-5′，按5′→3′书写后为5′-ATGCAT-3′，将T替换为U后得到5′-AUGCAU-3′。'},
  {cat:'翻译', q:'下列关于翻译过程的叙述，正确的是？', options:['翻译的场所是细胞核','一种tRNA只能识别并转运一种氨基酸','mRNA上的密码子与氨基酸一一对应','翻译时核糖体沿mRNA从3′端向5′端移动'], answer:1, explain:'翻译的场所是核糖体（细胞质中）；一种tRNA只能识别并转运一种氨基酸（但一种氨基酸可由多种tRNA转运）；mRNA上的密码子有64种，其中3种是终止密码子不编码氨基酸，所以不是一一对应；翻译时核糖体沿mRNA从5′端向3′端移动。'},
  {cat:'密码子', q:'下列关于密码子的叙述，错误的是？', options:['密码子是mRNA上决定一个氨基酸的3个相邻碱基','一种氨基酸可以由多种密码子决定','密码子具有通用性，所有生物共用一套密码子','密码子共有64种，其中编码氨基酸的有64种'], answer:3, explain:'密码子共有64种，其中3种（UAA、UAG、UGA）是终止密码子，不编码氨基酸，所以编码氨基酸的密码子有61种。一种氨基酸可以由一种或多种密码子决定（简并性）；密码子具有通用性（几乎所有生物共用一套密码子）。'},
  {cat:'密码子', q:'mRNA上的碱基序列为5′-AUGGGUUAA-3′，这段序列编码的氨基酸数为？', options:['2个','3个','4个','1个'], answer:0, explain:'每3个相邻碱基构成一个密码子。AUG（甲硫氨酸，起始密码子）、GGU（甘氨酸）、UAA（终止密码子，不编码氨基酸）。所以编码2个氨基酸。'},
  {cat:'中心法则', q:'下列关于中心法则的叙述，错误的是？', options:['中心法则包括DNA复制、转录、翻译等过程','逆转录是以RNA为模板合成DNA的过程','RNA复制只发生在RNA病毒中','中心法则中遗传信息只能从DNA流向RNA，不能反向'], answer:3, explain:'中心法则中遗传信息可以从DNA流向RNA（转录），也可以从RNA流向DNA（逆转录），逆转录过程发生在某些RNA病毒（如HIV）中。中心法则包括DNA复制、转录、翻译、逆转录、RNA复制等过程。'},
  {cat:'中心法则', q:'HIV病毒侵入人体后，其遗传信息的流动方向是？', options:['DNA→RNA→蛋白质','RNA→DNA→RNA→蛋白质','RNA→蛋白质','DNA→蛋白质'], answer:1, explain:'HIV是逆转录病毒，其遗传物质是RNA。侵入宿主细胞后，先在逆转录酶作用下以RNA为模板合成DNA（逆转录），然后DNA整合到宿主基因组中，再经转录和翻译合成病毒蛋白。所以流动方向是RNA→DNA→RNA→蛋白质。'},
  {cat:'基因突变', q:'镰刀型细胞贫血症的根本原因是？', options:['血红蛋白中一个氨基酸发生改变','DNA中一个碱基对发生替换','mRNA中一个碱基发生改变','染色体结构发生变异'], answer:1, explain:'镰刀型细胞贫血症的根本原因是DNA中一个碱基对发生替换（T-A替换为A-T），导致转录出的mRNA中一个密码子改变（GAG变为GUG），翻译出的血红蛋白中一个氨基酸改变（谷氨酸变为缬氨酸）。根本原因是基因突变（碱基对替换），直接原因是蛋白质结构改变。'},
  {cat:'基因突变', q:'下列关于基因突变的叙述，正确的是？', options:['基因突变一定会导致生物性状改变','基因突变只能发生在细胞分裂间期','基因突变是新基因产生的途径','基因突变对生物体都是有害的'], answer:2, explain:'基因突变不一定导致生物性状改变（如密码子简并性、隐性突变等）；基因突变可以发生在个体发育的任何时期，主要发生在细胞分裂间期DNA复制时；基因突变是新基因产生的途径，是生物变异的根本来源；基因突变对生物体大多是有害的，也有有利或中性的。'},
  {cat:'基因表达', q:'下列关于基因表达的叙述，错误的是？', options:['基因表达包括转录和翻译两个过程','基因表达具有选择性，不同细胞中表达的基因不同','基因表达的产物都是蛋白质','基因表达受环境因素影响'], answer:2, explain:'基因表达的产物不都是蛋白质，有些基因表达的产物是RNA（如tRNA、rRNA、microRNA等），这些基因称为非编码RNA基因。基因表达包括转录和翻译（编码蛋白质的基因）；具有选择性（细胞分化的实质是基因的选择性表达）；受环境因素影响。'},
  {cat:'基因表达', q:'胰岛B细胞能合成胰岛素，而胰岛A细胞不能，其根本原因是？', options:['两种细胞中含有的基因不同','两种细胞中含有的mRNA不同','两种细胞中核糖体的结构不同','两种细胞中tRNA的种类不同'], answer:1, explain:'胰岛B细胞和胰岛A细胞都来自同一个受精卵的有丝分裂，含有相同的基因。两种细胞功能不同的根本原因是基因的选择性表达，即转录出的mRNA不同，进而翻译出的蛋白质不同。核糖体结构和tRNA种类在不同细胞中是相同的。'},
  {cat:'DNA结构', q:'在一个DNA分子中，腺嘌呤与胸腺嘧啶之和占全部碱基的42%，则鸟嘌呤占全部碱基的？', options:['29%','21%','42%','58%'], answer:0, explain:'根据碱基互补配对原则，A=T，G=C。A+T=42%，所以G+C=100%-42%=58%，G=58%/2=29%。'},
  {cat:'DNA复制', q:'用¹⁵N标记的一个DNA分子放在含有¹⁴N的培养基中复制3次，则含¹⁵N的DNA分子占全部DNA分子的比例和含¹⁵N的单链占全部单链的比例分别是？', options:['1/4和1/8','1/8和1/16','1/2和1/4','1/4和1/4'], answer:0, explain:'复制3次产生8个DNA分子（16条单链）。根据半保留复制，含¹⁵N的DNA分子有2个（每个含一条¹⁵N链和一条¹⁴N链），比例为2/8=1/4。含¹⁵N的单链有2条，比例为2/16=1/8。'},
  {cat:'转录', q:'真核细胞中，mRNA从细胞核进入细胞质需要穿过几层生物膜？', options:['0层','1层','2层','4层'], answer:0, explain:'mRNA是大分子，通过核孔从细胞核进入细胞质，不穿过生物膜（核孔是核膜上的通道，物质通过核孔进出不穿过膜结构）。所以穿过0层生物膜。'},
  {cat:'翻译', q:'一条多肽链中有100个氨基酸，则作为合成该多肽链模板的mRNA分子和用来转录mRNA的DNA分子至少分别有碱基多少个？', options:['300个和600个','300个和300个','600个和600个','100个和200个'], answer:0, explain:'mRNA上3个相邻碱基决定一个氨基酸，所以100个氨基酸至少需要mRNA上100×3=300个碱基（不考虑终止密码子）。DNA是双链，转录时以一条链为模板，所以DNA分子至少有300×2=600个碱基。'},
  {cat:'基因突变', q:'某基因模板链的一个密码子由3′-TAC-5′变为3′-TCC-5′，该突变属于？', options:['同义突变','错义突变','无义突变','移码突变'], answer:1, explain:'原模板密码子3′-TAC-5′转录为mRNA密码子5′-AUG-3′，突变后3′-TCC-5′转录为5′-AGG-3′，氨基酸由甲硫氨酸变为精氨酸，因此属于错义突变。'},
  {cat:'中心法则', q:'下列过程中，不发生碱基互补配对的是？', options:['DNA复制','转录','翻译','蛋白质的加工'], answer:3, explain:'DNA复制时母链与子链之间碱基互补配对；转录时DNA模板链与mRNA之间碱基互补配对；翻译时mRNA密码子与tRNA反密码子之间碱基互补配对。蛋白质的加工（如折叠、糖基化等）不涉及碱基互补配对。'},
  {cat:'基因表达', q:'下列关于基因、蛋白质与性状关系的叙述，错误的是？', options:['基因可以通过控制酶的合成来控制代谢过程，进而控制生物性状','基因可以通过控制蛋白质的结构直接控制生物性状','生物性状只由基因决定，与环境无关','一个基因可以影响多个性状，一个性状也可以由多个基因控制'], answer:2, explain:'生物性状是基因与环境共同作用的结果，不是只由基因决定。例如，同一株水毛茛，裸露在空气中的叶和浸在水中的叶形态不同，说明环境影响性状表现。基因控制性状有两条途径：间接控制（通过控制酶的合成控制代谢）和直接控制（通过控制蛋白质结构）。基因与性状之间不是简单的一一对应关系。'},
  {cat:'密码子', q:'下列关于tRNA的叙述，正确的是？', options:['tRNA的合成场所是核糖体','一种tRNA可以转运多种氨基酸','tRNA上的反密码子与mRNA上的密码子互补配对','tRNA含有密码子'], answer:2, explain:'tRNA是在细胞核中通过转录合成的（不是核糖体）；一种tRNA只能识别并转运一种氨基酸（专一性）；tRNA一端的反密码子与mRNA上的密码子互补配对；密码子位于mRNA上，tRNA上含有的是反密码子。'},
  {cat:'DNA结构', q:'下列关于DNA分子特性的叙述，错误的是？', options:['DNA分子具有多样性','DNA分子具有特异性','DNA分子具有稳定性','DNA分子的多样性由碱基配对方式决定'], answer:3, explain:'DNA分子的多样性由碱基对的排列顺序决定，而不是由碱基配对方式决定（碱基配对方式只有A-T、G-C两种，是固定的）。DNA分子具有多样性（碱基排列顺序千变万化）、特异性（每个DNA分子有特定的碱基排列顺序）、稳定性（双螺旋结构和氢键维持稳定）。'},
  {cat:'转录', q:'下列关于RNA的叙述，错误的是？', options:['RNA分子通常是单链结构','RNA分子中不含胸腺嘧啶','mRNA是翻译的直接模板','RNA分子中不含有氢键'], answer:3, explain:'RNA分子通常是单链结构，但tRNA分子中存在局部双链区域，含有氢键（如tRNA的三叶草结构中，互补碱基之间形成氢键）。RNA分子中含有的碱基是A、U、G、C，不含T（胸腺嘧啶）；mRNA是翻译的直接模板。'}
];

var quizState = {
  current: 0,
  answers: {},  // {index: selectedOption}
  submitted: {},  // {index: true/false}
  filter: 'all'
};

function renderQuiz(){
  var container = document.getElementById('quizContainer');
  if(!container) return;

  var filtered = getFilteredQuestions();
  if(filtered.length === 0){
    container.innerHTML = '<div style="text-align:center;padding:40px;color:var(--ink-soft)">该分类暂无题目</div>';
    return;
  }

  var q = filtered[quizState.current];
  var qIndex = QUIZ_DATA.indexOf(q);
  var isSubmitted = quizState.submitted[qIndex];
  var selected = quizState.answers[qIndex];
  var total = filtered.length;
  var answered = Object.keys(quizState.submitted).length;
  var correct = 0;
  for(var k in quizState.submitted){
    if(quizState.submitted[k] && quizState.answers[k] === QUIZ_DATA[k].answer) correct++;
  }

  var html = '';

  // 统计栏
  html += '<div class="quiz-header">';
  html += '<div class="quiz-stats">';
  html += '<div class="quiz-stat"><div class="v">'+total+'</div><div class="k">总题数</div></div>';
  html += '<div class="quiz-stat"><div class="v">'+answered+'</div><div class="k">已作答</div></div>';
  html += '<div class="quiz-stat"><div class="v v-ok">'+correct+'</div><div class="k">正确</div></div>';
  html += '<div class="quiz-stat"><div class="v v-bad">'+(answered-correct)+'</div><div class="k">错误</div></div>';
  html += '</div>';
  html += '<button class="btn btn-ghost" onclick="resetQuiz()">重新开始</button>';
  html += '</div>';

  // 进度条
  html += '<div class="quiz-progress-bar"><div class="quiz-progress-fill" style="width:'+((quizState.current+1)/total*100)+'%"></div></div>';

  // 分类筛选
  var cats = getCategories();
  html += '<div class="quiz-filter">';
  html += '<button class="'+(quizState.filter==='all'?'active':'')+'" onclick="filterQuiz(\'all\')">全部</button>';
  for(var c=0;c<cats.length;c++){
    html += '<button class="'+(quizState.filter===cats[c]?'active':'')+'" onclick="filterQuiz(\''+cats[c]+'\')">'+cats[c]+'</button>';
  }
  html += '</div>';

  // 题目卡片
  html += '<div class="quiz-card">';
  html += '<div class="quiz-qmeta">';
  html += '<span class="quiz-qnum">第 '+(quizState.current+1)+' / '+total+' 题</span>';
  html += '<span class="quiz-qcat">'+q.cat+'</span>';
  if(isSubmitted){
    html += '<span class="quiz-qcat" style="background:'+(selected===q.answer?'rgba(42,122,42,.15)':'rgba(154,42,42,.15)')+';color:'+(selected===q.answer?'#2A7A2A':'#9A2A2A')+'">'+(selected===q.answer?'回答正确 ✓':'回答错误 ✗')+'</span>';
  }
  html += '</div>';
  html += '<div class="quiz-qtext">'+q.q+'</div>';
  html += '<div class="quiz-options">';
  for(var i=0;i<q.options.length;i++){
    var optClass = 'quiz-option';
    if(selected === i) optClass += ' selected';
    if(isSubmitted){
      optClass += ' disabled';
      if(i === q.answer) optClass += ' correct';
      if(i === selected && selected !== q.answer) optClass += ' wrong';
    }
    html += '<button type="button" class="'+optClass+'" '+(isSubmitted?'disabled ':'')+'aria-pressed="'+(selected===i?'true':'false')+'" onclick="'+(isSubmitted?'':'selectQuizOption('+i+')')+'">';
    html += '<span class="letter">'+String.fromCharCode(65+i)+'</span>';
    html += '<span class="text">'+q.options[i]+'</span>';
    html += '</button>';
  }
  html += '</div>';

  // 解析
  html += '<div class="quiz-explain '+(isSubmitted?'show':'')+'">';
  html += '<h4>答案与解析</h4>';
  html += '<p><strong>正确答案：</strong>'+String.fromCharCode(65+q.answer)+'。'+q.options[q.answer]+'</p>';
  html += '<p><strong>解析：</strong>'+q.explain+'</p>';
  html += '</div>';

  // 操作按钮
  html += '<div class="quiz-actions">';
  html += '<div class="left">';
  html += '<button class="btn btn-ghost" onclick="prevQuiz()" '+(quizState.current===0?'disabled':'')+'>上一题</button>';
  html += '<button class="btn btn-primary" onclick="nextQuiz()">'+(quizState.current===total-1?'查看成绩':'下一题')+'</button>';
  html += '</div>';
  if(!isSubmitted && selected !== undefined){
    html += '<button class="btn btn-primary" style="background:#2A7A2A" onclick="submitQuiz()">提交答案</button>';
  }
  html += '</div>';
  html += '</div>';

  // 题目导航
  html += '<div class="quiz-nav">';
  for(var n=0;n<total;n++){
    var navQ = filtered[n];
    var navIdx = QUIZ_DATA.indexOf(navQ);
    var navClass = '';
    if(n === quizState.current) navClass = 'current';
    else if(quizState.submitted[navIdx]){
      navClass = quizState.answers[navIdx] === QUIZ_DATA[navIdx].answer ? 'correct' : 'wrong';
    }
    html += '<button class="'+navClass+'" onclick="goQuiz('+n+')">'+(n+1)+'</button>';
  }
  html += '</div>';

  container.innerHTML = html;
}

function getCategories(){
  var cats = [];
  for(var i=0;i<QUIZ_DATA.length;i++){
    if(cats.indexOf(QUIZ_DATA[i].cat) === -1) cats.push(QUIZ_DATA[i].cat);
  }
  return cats;
}

function getFilteredQuestions(){
  if(quizState.filter === 'all') return QUIZ_DATA.slice();
  return QUIZ_DATA.filter(function(q){return q.cat === quizState.filter;});
}

function filterQuiz(cat){
  quizState.filter = cat;
  quizState.current = 0;
  renderQuiz();
}

function selectQuizOption(idx){
  var filtered = getFilteredQuestions();
  var q = filtered[quizState.current];
  var qIndex = QUIZ_DATA.indexOf(q);
  quizState.answers[qIndex] = idx;
  renderQuiz();
}

function submitQuiz(){
  var filtered = getFilteredQuestions();
  var q = filtered[quizState.current];
  var qIndex = QUIZ_DATA.indexOf(q);
  if(quizState.answers[qIndex] === undefined) return;
  quizState.submitted[qIndex] = true;
  renderQuiz();
}

function prevQuiz(){
  if(quizState.current > 0){
    quizState.current--;
    renderQuiz();
  }
}

function nextQuiz(){
  var filtered = getFilteredQuestions();
  if(quizState.current < filtered.length - 1){
    quizState.current++;
    renderQuiz();
  } else {
    showQuizResult();
  }
}

function goQuiz(n){
  quizState.current = n;
  renderQuiz();
}

function resetQuiz(){
  quizState.current = 0;
  quizState.answers = {};
  quizState.submitted = {};
  quizState.filter = 'all';
  renderQuiz();
}

function showQuizResult(){
  var container = document.getElementById('quizContainer');
  var total = QUIZ_DATA.length;
  var answered = Object.keys(quizState.submitted).length;
  var correct = 0;
  for(var k in quizState.submitted){
    if(quizState.submitted[k] && quizState.answers[k] === QUIZ_DATA[k].answer) correct++;
  }
  var percent = answered > 0 ? Math.round(correct / answered * 100) : 0;

  var html = '<div class="quiz-result">';
  html += '<div class="score">'+percent+'<span style="font-size:28px">分</span></div>';
  html += '<div class="score-label">答题正确率（已作答 '+answered+' / '+total+' 题）</div>';
  html += '<div class="detail">';
  html += '<div class="detail-item"><div class="n right">'+correct+'</div><div class="l">答对</div></div>';
  html += '<div class="detail-item"><div class="n wrong">'+(answered-correct)+'</div><div class="l">答错</div></div>';
  html += '<div class="detail-item"><div class="n">'+(total-answered)+'</div><div class="l">未作答</div></div>';
  html += '</div>';
  html += '<div style="margin-top:30px;display:flex;gap:12px;justify-content:center;flex-wrap:wrap">';
  html += '<button class="btn btn-primary" onclick="renderQuiz()">继续练习</button>';
  html += '<button class="btn btn-ghost" onclick="resetQuiz()">重新开始</button>';
  html += '<button class="btn btn-ghost" onclick="showWrongQuestions()">查看错题</button>';
  html += '</div>';
  html += '</div>';

  container.innerHTML = html;
}

function showWrongQuestions(){
  var wrongIndices = [];
  for(var k in quizState.submitted){
    if(quizState.submitted[k] && quizState.answers[k] !== QUIZ_DATA[k].answer){
      wrongIndices.push(parseInt(k));
    }
  }
  if(wrongIndices.length === 0){
    alert('没有错题，太棒了！');
    return;
  }
  // 只显示错题
  quizState.filter = 'wrong';
  quizState.current = 0;
  // 临时修改QUIZ_DATA为只含错题的方式：用filter特殊处理
  renderWrongQuiz(wrongIndices);
}

function renderWrongQuiz(wrongIndices){
  var container = document.getElementById('quizContainer');
  var html = '<div style="margin-bottom:16px;display:flex;justify-content:space-between;align-items:center">';
  html += '<h3 style="color:var(--teal)">错题回顾（共'+wrongIndices.length+'题）</h3>';
  html += '<button class="btn btn-ghost" onclick="resetQuiz()">返回全部题目</button>';
  html += '</div>';

  for(var i=0;i<wrongIndices.length;i++){
    var idx = wrongIndices[i];
    var q = QUIZ_DATA[idx];
    html += '<div class="quiz-card" style="margin-bottom:16px">';
    html += '<div class="quiz-qmeta"><span class="quiz-qnum">第 '+(idx+1)+' 题</span><span class="quiz-qcat">'+q.cat+'</span></div>';
    html += '<div class="quiz-qtext">'+q.q+'</div>';
    html += '<div class="quiz-options">';
    for(var j=0;j<q.options.length;j++){
      var cls = 'quiz-option disabled';
      if(j === q.answer) cls += ' correct';
      if(j === quizState.answers[idx]) cls += ' wrong';
      html += '<button type="button" class="'+cls+'" disabled><span class="letter">'+String.fromCharCode(65+j)+'</span><span class="text">'+q.options[j]+'</span></button>';
    }
    html += '</div>';
    html += '<div class="quiz-explain show"><h4>答案与解析</h4><p><strong>正确答案：</strong>'+String.fromCharCode(65+q.answer)+'。'+q.options[q.answer]+'</p><p><strong>你的答案：</strong>'+String.fromCharCode(65+quizState.answers[idx])+'</p><p><strong>解析：</strong>'+q.explain+'</p></div>';
    html += '</div>';
  }

  container.innerHTML = html;
}
