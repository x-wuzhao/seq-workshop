var CODON_TABLES = {
  standard: {
    name: '标准密码表（核基因）',
    map: CODON,
    starts: ['ATG']
  },
  vertebrate_mito: {
    name: '脊椎动物线粒体密码表',
    map: (function () {
      var m = {};
      for (var k in CODON) m[k] = CODON[k];
      m.TGA = 'W';      // 线粒体：TGA 不是终止，而是色氨酸
      m.AGA = '*';      // AGA/AGG 在线粒体中为终止
      m.AGG = '*';
      m.ATA = 'M';      // ATA 在线粒体中为甲硫氨酸
      return m;
    })(),
    starts: ['ATG', 'ATT', 'ATC', 'ATA', 'GTG']
  }
};
var activeCodonTable = 'standard';

function codonMap() { return CODON_TABLES[activeCodonTable].map; }

// 覆盖基础版：按当前遗传密码表翻译
function translateFrame(seq, frame) {
  var m = codonMap(), p = '';
  for (var i = frame; i + 3 <= seq.length; i += 3) p += (m[seq.substr(i, 3)] || 'X');
  return p;
}

// 覆盖基础版：起始密码子随密码表变化
function findORFs(seq, frame) {
  return SeqCore.findORFs(seq,frame,codonMap(),CODON_TABLES[activeCodonTable].starts);
}

/* ---------------- 2. 最近邻热力学（SantaLucia 1998） ---------------- */
var NN_H = {
  AA: -7.9, TT: -7.9, AT: -7.2, TA: -7.2, CA: -8.5, TG: -8.5, GT: -8.4, AC: -8.4,
  CT: -7.8, AG: -7.8, GA: -8.2, TC: -8.2, CG: -10.6, GC: -9.8, GG: -8.0, CC: -8.0
};
var NN_S = {
  AA: -22.2, TT: -22.2, AT: -20.4, TA: -21.3, CA: -22.7, TG: -22.7, GT: -22.4, AC: -22.4,
  CT: -21.0, AG: -21.0, GA: -22.2, TC: -22.2, CG: -27.2, GC: -24.4, GG: -19.9, CC: -19.9
};
var NN_G37 = {
  AA: -1.00, TT: -1.00, AT: -0.88, TA: -0.58, CA: -1.45, TG: -1.45, GT: -1.44, AC: -1.44,
  CT: -1.28, AG: -1.28, GA: -1.30, TC: -1.30, CG: -2.17, GC: -2.17, GG: -1.84, CC: -1.84
};
var R_GAS = 1.987;          // cal/(mol·K)
var PRIMER_CT = 0.25e-6;    // 引物总浓度 0.25 µM（Primer3 默认量级）
var NA_CONC = 0.05;         // 单价阳离子 50 mM

/* 最近邻 Tm：ΔH 千卡/摩尔、ΔS 卡/摩尔/K，含末端修正与盐校正 */
function tmNN(seq, opts) {
  opts = opts || {};
  var na = opts.na || NA_CONC;
  var ct = opts.ct || PRIMER_CT;
  if (!seq || seq.length < 2) return NaN;
  var dh = 0, ds = 0; // SantaLucia 1998 / Biopython DNA_NN3 initiation
  for (var i = 0; i + 1 < seq.length; i++) {
    var k = seq.substr(i, 2);
    if (NN_H[k] === undefined) return NaN;  // 含 N 等简并碱基时不计算
    dh += NN_H[k];
    ds += NN_S[k];
  }
  var selfComplementary = reverseComplement(seq) === seq;
  if (selfComplementary) ds -= 1.4;
  // 末端 A/T 罚分（末端 G/C 更稳定）
  [seq.charAt(0), seq.charAt(seq.length - 1)].forEach(function (b) {
    if (b === 'A' || b === 'T') { dh += 2.3; ds += 4.1; }
    else { dh += 0.1; ds -= 2.8; }
  });
  // 盐校正：ds += 0.368 × (N-1) × ln[Na+]
  ds += 0.368 * (seq.length - 1) * Math.log(na);
  var tm = (dh * 1000) / (ds + R_GAS * Math.log(ct / (selfComplementary ? 1 : 4))) - 273.15;
  return tm;
}

/* 两条序列反向平行配对的最小 ΔG37（用于二聚体 / 发夹近似） */
function pairDG(a, b) {
  var best = 0;
  for (var off = -b.length + 1; off < a.length; off++) {
    var sum = 0, run = 0;
    for (var i = 0; i < a.length; i++) {
      var j = i + off;
      if (j < 0 || j >= b.length) { run = 0; continue; }
      // a 的 5'→3' 与 b 的反向互补配对
      var x = a.charAt(i), y = b.charAt(b.length - 1 - j);
      var pairOK = (x === 'A' && y === 'T') || (x === 'T' && y === 'A') ||
                   (x === 'G' && y === 'C') || (x === 'C' && y === 'G');
      if (!pairOK) { run = 0; continue; }
      run++;
      if (run >= 2) {
        var d = NN_G37[a.charAt(i - 1) + x];
        if (d !== undefined) sum += d; else sum += -1.0;
      }
    }
    if (sum < best) best = sum;
  }
  return best;   // kcal/mol，越负越容易形成二聚体
}

function selfDimerDG(seq) { return pairDG(seq, seq); }

function hairpinDG(seq) {
  var best = 0;
  for (var stem = 3; stem <= Math.min(9, Math.floor(seq.length / 2) - 3); stem++) {
    for (var loop = 3; loop <= 9; loop++) {
      for (var s = 0; s + stem * 2 + loop <= seq.length; s++) {
        var left = seq.substr(s, stem);
        var right = seq.substr(s + stem + loop, stem);
        var rc = '';
        for (var i = right.length - 1; i >= 0; i--) {
          rc += { A: 'T', T: 'A', G: 'C', C: 'G', N: 'N' }[right.charAt(i)] || 'N';
        }
        if (left !== rc) continue;
        var dg = -1.4;
        for (var j = 0; j + 1 < stem; j++) {
          var d = NN_G37[left.charAt(j) + left.charAt(j + 1)];
          if (d !== undefined) dg += d;
        }
        dg += (loop <= 4 ? 4.5 : loop <= 6 ? 4.0 : 3.5);
        if (dg < best) best = dg;
      }
    }
  }
  return best;
}

/* ---------------- 3. 引物设计 v2 ---------------- */
function gcOf(s) {
  var g = 0;
  for (var i = 0; i < s.length; i++) { var c = s.charAt(i); if (c === 'G' || c === 'C') g++; }
  return s.length ? g / s.length * 100 : 0;
}

function countOccurrences(template, probe) {
  var rc = reverseComplement(probe), n = 0;
  for (var i = 0; i + probe.length <= template.length; i++) {
    var sub = template.substr(i, probe.length);
    if (sub === probe || sub === rc) n++;
  }
  return n;
}

function designPrimersV2() {
  var seq = window._currentSeq || '';
  var out = document.getElementById('primerResult');
  if (!seq) { out.innerHTML = '<div class="msg info">请先在上方输入序列并点击“开始分析”</div>'; return; }

  var productLen = parseInt(document.getElementById('primerProduct').value, 10) || 300;
  var offset = parseInt(document.getElementById('primerOffset').value, 10) || 0;
  var targetTm = parseFloat(document.getElementById('primerTm').value) || 60;

  if (seq.length < productLen + 40) {
    out.innerHTML = '<div class="msg err">序列长度不足：设计 ' + productLen + ' bp 产物至少需要 ' + (productLen + 40) + ' bp 模板。</div>';
    return;
  }

  var requestedLen=Number(document.getElementById('primerLen').value);
  if(!Number.isInteger(requestedLen)||requestedLen<15||requestedLen>30||!Number.isInteger(productLen)||productLen<100||productLen>2000||!Number.isInteger(offset)||offset<0||offset>=seq.length||!Number.isFinite(targetTm)||targetTm<50||targetTm>70){
    out.textContent='请检查参数：引物 15–30 bp，目标 Tm 50–70°C，产物 100–2000 bp，偏移须在模板范围内。';return;
  }
  var lenList = [requestedLen];
  var fMin = Math.max(0, offset - 20);
  var fMax = Math.min(seq.length - productLen - 10, offset + 40);
  var fCands = [], i, li;

  for (var fStart = fMin; fStart <= fMax; fStart++) {
    for (li = 0; li < lenList.length; li++) {
      var fLen = lenList[li];
      var fSeq = seq.substr(fStart, fLen);
      if (fSeq.length < fLen || fSeq.indexOf('N') >= 0) continue;
      var fTm = tmNN(fSeq);
      if (!isFinite(fTm) || Math.abs(fTm - targetTm) > 8) continue;
      var fGC = gcOf(fSeq);
      if (fGC < 30 || fGC > 70) continue;
      var lastF = fSeq.charAt(fLen - 1);
      fCands.push({ start: fStart, len: fLen, seq: fSeq, tm: fTm, gc: fGC, gcClamp: (lastF === 'G' || lastF === 'C') });
    }
  }
  if (!fCands.length) {
    out.innerHTML = '<div class="msg err">没有找到 Tm 与 GC 合适的正向引物，请调整目标 Tm 或产物长度。</div>';
    return;
  }

  // 轻量评分：先用便宜指标粗排，再对前若干名做二聚体/发夹热力学计算
  var pairs = [];
  fCands.forEach(function (fc) {
    for (var rEnd = fc.start + productLen - 12; rEnd <= fc.start + productLen + 12; rEnd++) {
      if (rEnd > seq.length) break;
      for (var rj = 0; rj < lenList.length; rj++) {
        var rLen = lenList[rj];
        if (rEnd - rLen < fc.start + fc.len + 1) continue;
        var rSeq = reverseComplement(seq.substr(rEnd - rLen, rLen));
        if (rSeq.indexOf('N') >= 0) continue;
        var rTm = tmNN(rSeq);
        if (!isFinite(rTm) || Math.abs(rTm - targetTm) > 8) continue;
        var rGC = gcOf(rSeq);
        if (rGC < 30 || rGC > 70) continue;
        var cheap = Math.abs(fc.tm - targetTm) + Math.abs(rTm - targetTm) + Math.abs(fc.tm - rTm) * 1.5
          + Math.abs(fc.gc - 50) * 0.12 + Math.abs(rGC - 50) * 0.12
          + (fc.gcClamp ? 0 : 0.6);
        pairs.push({
          f: fc, rSeq: rSeq, rTm: rTm, rGC: rGC, rStart: rEnd - rLen + 1, rEnd: rEnd,
          product: rEnd - fc.start, cheap: cheap
        });
      }
    }
  });
  if (!pairs.length) {
    out.innerHTML = '<div class="msg err">在当前参数下没有找到成对引物，请放宽目标 Tm（如 58–62 °C）或调整产物长度。</div>';
    return;
  }
  pairs.sort(function (a, b) { return a.cheap - b.cheap; });

  var top = pairs.slice(0, 40).map(function (p) {
    var hetero = pairDG(p.f.seq, p.rSeq);
    var fHair = hairpinDG(p.f.seq), rHair = hairpinDG(p.rSeq);
    var fSelf = selfDimerDG(p.f.seq), rSelf = selfDimerDG(p.rSeq);
    p.hetero = hetero; p.fHair = fHair; p.rHair = rHair; p.fSelf = fSelf; p.rSelf = rSelf;
    p.penalty = p.cheap + Math.max(0, -hetero - 4) * 2 + Math.max(0, -fHair - 2) * 1.5
      + Math.max(0, -rHair - 2) * 1.5 + Math.max(0, -fSelf - 6) + Math.max(0, -rSelf - 6);
    return p;
  }).sort(function (a, b) { return a.penalty - b.penalty; }).slice(0, 5);

  top.forEach(function (p) {
    p.fHits = countOccurrences(seq, p.f.seq);
    p.rHits = countOccurrences(seq, p.rSeq);
  });
  // 地图/拖动用的扁平结构（候选表仍用 top 里的 f/r 嵌套结构）
  window._lastPrimers = top.map(function (p) {
    return {
      fStart: p.f.start + 1, fLen: p.f.len, fSeq: p.f.seq, fTm: p.f.tm, fGC: p.f.gc,
      rStart: p.rStart, rEnd: p.rEnd, rLen: p.rSeq ? p.rSeq.length : p.rLen, rSeq: p.rSeq, rTm: p.rTm, rGC: p.rGC,
      product: p.product, hetero: p.hetero, penalty: p.penalty, fHits: p.fHits, rHits: p.rHits
    };
  });

  var h = '<div class="primer-hint">按最近邻热力学模型（SantaLucia 1998）评分：Tm 匹配、GC 含量、3′ 端 GC 夹、引物二聚体与发夹 ΔG，取前 5 对候选。</div>';
  h += '<div class="primer-table-wrap"><table class="primer-table"><thead><tr>' +
    '<th>#</th><th>正向引物（5′→3′）</th><th>Tm</th><th>位置</th>' +
    '<th>反向引物（5′→3′）</th><th>Tm</th><th>位置</th><th>产物</th><th>ΔG(二聚体)</th><th>评分</th>' +
    '</tr></thead><tbody>';
  top.forEach(function (p, idx) {
    h += '<tr class="' + (idx === 0 ? 'best' : '') + '">' +
      '<td>' + (idx + 1) + (idx === 0 ? ' ★' : '') + '</td>' +
      '<td class="mono">' + p.f.seq + '</td><td>' + p.f.tm.toFixed(1) + '°C</td><td>' + (p.f.start + 1) + '–' + (p.f.start + p.f.len) + '</td>' +
      '<td class="mono">' + p.rSeq + '</td><td>' + p.rTm.toFixed(1) + '°C</td><td>' + p.rStart + '–' + p.rEnd + '</td>' +
      '<td>' + p.product + ' bp</td><td>' + p.hetero.toFixed(1) + '</td><td>' + p.penalty.toFixed(1) + '</td></tr>';
  });
  h += '</tbody></table></div>';

  var p0 = top[0];
  var warn = [];
  if (Math.abs(p0.f.tm - p0.rTm) > 2) warn.push('正反向 Tm 相差 ' + Math.abs(p0.f.tm - p0.rTm).toFixed(1) + ' °C（建议 < 2 °C）');
  if (p0.fHits > 1 || p0.rHits > 1) warn.push('候选引物在模板中出现 ' + Math.max(p0.fHits, p0.rHits) + ' 次，可能存在非特异扩增');
  if (p0.hetero < -6) warn.push('引物二聚体 ΔG = ' + p0.hetero.toFixed(1) + ' kcal/mol，建议换一对');
  if (p0.fHair < -3 || p0.rHair < -3) warn.push('存在较强的发夹结构（ΔG < −3 kcal/mol）');
  if (warn.length) {
    h += '<div class="primer-card" style="border-color:#F0C9C9;background:#FDF5F5"><div class="primer-name primer-warn">⚠ 首选引物仍有的风险</div><ul style="margin:8px 0 0 18px;font-size:13px;color:var(--ink-soft)">' +
      warn.map(function (w) { return '<li>' + w + '</li>'; }).join('') + '</ul></div>';
  } else {
    h += '<div class="primer-card" style="border-color:#C4DEDC;background:#F0F8F7"><div class="primer-name" style="color:#2A7A2A">✓ 首选引物通过基础筛选</div><p style="font-size:13px;color:var(--ink-soft);margin-top:6px">二聚体/发夹 ΔG 与 Tm 匹配均在常用阈值内。实验前仍建议用 Primer-BLAST 检查全基因组特异性。</p></div>';
  }
  h += '<p class="primer-hint" style="margin-top:10px">提示：在“序列地图”的引物轨道上可直接拖动紫色箭头，实时查看 Tm 与产物长度变化；此处的候选表始终保留原始推荐结果。</p>';
  out.innerHTML = h;
  if (typeof renderSeqMap === 'function') renderSeqMap();
}

/* ---------------- 4. 酶切 v2：线性/环状 + 双酶切 ---------------- */
var ENZ_BY_NAME = {};
(function () {
  for (var i = 0; i < ENZYMES.length; i++) ENZ_BY_NAME[ENZYMES[i][0]] = ENZYMES[i];
})();

function enzymeOverhang(site, c1) {
  var L = site.length, c2 = L - c1;
  if (c1 === c2) return { type: '平端', len: 0 };
  return { type: c1 < c2 ? '5′ 突出' : '3′ 突出', len: Math.abs(c1 - c2), seq: site.substr(Math.min(c1, c2), Math.abs(c1 - c2)) };
}

/* 汇总若干酶的切点，返回片段表；circular 时首尾相接 */
function digestFragments(seqLen, hits, circular) {
  var cuts = [];
  for (var i = 0; i < hits.length; i++) {
    var h = hits[i];
    if (h.cut > 0 && h.cut < seqLen) cuts.push(h);
  }
  cuts.sort(function (a, b) { return a.cut - b.cut; });
  var frags = [], k;
  if (!cuts.length) {
    frags.push({ start: 1, end: seqLen, len: seqLen, left: null, right: null, circular: false });
    return { frags: frags, cuts: cuts };
  }
  if (!circular) {
    frags.push({ start: 1, end: cuts[0].cut, len: cuts[0].cut, left: null, right: cuts[0] });
    for (k = 0; k + 1 < cuts.length; k++) {
      frags.push({ start: cuts[k].cut + 1, end: cuts[k + 1].cut, len: cuts[k + 1].cut - cuts[k].cut, left: cuts[k], right: cuts[k + 1] });
    }
    frags.push({ start: cuts[cuts.length - 1].cut + 1, end: seqLen, len: seqLen - cuts[cuts.length - 1].cut, left: cuts[cuts.length - 1], right: null });
  } else {
    for (k = 0; k < cuts.length; k++) {
      var next = cuts[(k + 1) % cuts.length];
      var len = (next.cut - cuts[k].cut + seqLen) % seqLen || seqLen;
      frags.push({ start: cuts[k].cut + 1, end: next.cut, len: len, left: cuts[k], right: next });
    }
  }
  return { frags: frags, cuts: cuts };
}

function digestLabel(h) {
  if (!h) return '—';
  var e = h.enz[0];
  var oh = enzymeOverhang(h.enz[1], h.enz[2]);
  return e + (oh.len ? '（' + oh.type + ' ' + oh.len + ' nt）' : '（平端）');
}

function renderDigest(seq, res) {
  var sel1 = document.getElementById('gelEnzSelect');
  var sel2 = document.getElementById('digestEnz2');
  var topo = document.getElementById('digestTopo');
  var enzNames = [];
  if (sel1 && sel1.value) enzNames.push(sel1.value);
  if (sel2 && sel2.value && sel2.value !== sel1.value) enzNames.push(sel2.value);
  var circular = topo && topo.value === 'circular';

  var hits = [];
  enzNames.forEach(function (n) {
    for (var i = 0; i < res.length; i++) {
      if (res[i].enz[0] === n) hits = hits.concat(res[i].hits.map(function (h) { return { start: h.start, cut: h.cut, rc: h.rc, enz: res[i].enz }; }));
    }
  });

  var d = digestFragments(seq.length, hits, circular);
  var tableBox = document.getElementById('digestTable');
  var title = enzNames.length ? enzNames.join(' + ') + (circular ? ' · 环状' : ' · 线性') : '未选择酶';

  var h = '<div class="digest-head"><strong>' + title + '</strong><span>共 ' + d.frags.length + ' 个片段' +
    (circular ? '（环状 DNA 被切 ' + d.cuts.length + ' 次）' : '') + '</span></div>';
  if (!hits.length) {
    h += '<p class="empty-note">所选酶在此序列上没有切点，可换一种酶或改选“环状”。</p>';
    if (tableBox) tableBox.innerHTML = h;
    renderGel(seq, res, enzNames);
    return;
  }
  h += '<table class="enz-table digest-table"><thead><tr><th>片段</th><th>长度 (bp)</th><th>位置</th><th>左端</th><th>右端</th></tr></thead><tbody>';
  d.frags.sort(function (a, b) { return b.len - a.len; });
  d.frags.forEach(function (f, i) {
    h += '<tr><td>' + (i + 1) + '</td><td><b>' + f.len + '</b></td><td>' + f.start + '–' + f.end + '</td>' +
      '<td>' + digestLabel(f.left) + '</td><td>' + digestLabel(f.right) + '</td></tr>';
  });
  h += '</tbody></table>';
  if (tableBox) tableBox.innerHTML = h;
  renderGel(seq, res, enzNames);
}

function renderGel(seq, res, enzNames) {
  var box = document.getElementById('gelBox');
  if (!box) return;
  var sel1 = document.getElementById('gelEnzSelect');
  var sel2 = document.getElementById('digestEnz2');
  var topo = document.getElementById('digestTopo');
  var names = [];
  if (sel1 && sel1.value) names.push(sel1.value);
  if (sel2 && sel2.value && sel2.value !== sel1.value) names.push(sel2.value);
  var circular = topo && topo.value === 'circular';

  var lanes = [];
  names.forEach(function (n, idx) {
    var hits = [];
    for (var i = 0; i < res.length; i++) if (res[i].enz[0] === n) hits = res[i].hits;
    lanes.push({ label: names.length > 1 ? (idx === 0 ? n : n + '（双）') : n, hits: hits });
  });
  if (!lanes.length) { box.innerHTML = '<p class="empty-note">选择酶后显示虚拟凝胶。</p>'; return; }

  // 所有片段长度，用于确定对数刻度范围
  var allSizes = [];
  lanes.forEach(function (l) {
    var d = digestFragments(seq.length, l.hits, circular);
    d.frags.forEach(function (f) { allSizes.push(f.len); });
  });
  var maxBp = Math.max.apply(null, allSizes.concat([seq.length, 600]));
  var minBp = Math.max(5, Math.min.apply(null, allSizes.concat([60])) * 0.8);

  var laneW = 112, leftPad = 96;
  var W = leftPad + laneW * (lanes.length + 1) + 24;
  var H = 330, yTop = 52, yBot = 292;
  function Y(bp) {
    var b = Math.max(minBp, Math.min(maxBp, bp));
    return yTop + (Math.log10(maxBp) - Math.log10(b)) / (Math.log10(maxBp) - Math.log10(minBp)) * (yBot - yTop);
  }

  var markerX = leftPad + laneW * 0.5;
  var svg = '<svg width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + ' ' + H + '" xmlns="http://www.w3.org/2000/svg">';

  // 点样孔
  function well(x) {
    return '<rect x="' + (x - 30) + '" y="30" width="60" height="9" rx="2" fill="#0d1a16" stroke="rgba(127,232,180,.25)"/>';
  }
  svg += well(markerX);
  svg += '<text x="' + markerX + '" y="22" text-anchor="middle" fill="#BFE0DD" font-size="11">Marker</text>';

  var marker = [100, 200, 300, 500, 700, 1000, 1500, 2000, 3000, 5000].filter(function (b) {
    return b >= minBp && b <= maxBp;
  });
  marker.forEach(function (bp) {
    var y = Y(bp);
    svg += '<rect x="' + (markerX - 26) + '" y="' + (y - 2) + '" width="52" height="4" rx="2" fill="#A7E0DE" opacity=".85"/>';
    svg += '<text x="' + (markerX - 32) + '" y="' + (y + 3.5) + '" text-anchor="end" fill="#9FB6B2" font-size="10">' + bp + '</text>';
  });

  lanes.forEach(function (l, i) {
    var cx = leftPad + laneW * (i + 1.5);
    svg += well(cx);
    svg += '<text x="' + cx + '" y="22" text-anchor="middle" fill="#BFE0DD" font-size="11">' + l.label + '</text>';
    svg += '<line x1="' + cx + '" y1="' + yTop + '" x2="' + cx + '" y2="' + yBot + '" stroke="rgba(159,198,194,.18)"/>';
    var d = digestFragments(seq.length, l.hits, circular);
    // 相近长度合并（±5%），避免条带重叠
    var groups = [];
    d.frags.forEach(function (f) {
      var g = null;
      for (var k = 0; k < groups.length; k++) {
        if (Math.abs(groups[k].len - f.len) / Math.max(groups[k].len, f.len) < 0.05) { g = groups[k]; break; }
      }
      if (g) { g.n++; g.len = Math.round((g.len * (g.n - 1) + f.len) / g.n); }
      else groups.push({ len: f.len, n: 1 });
    });
    groups.sort(function (a, b) { return b.len - a.len; });
    groups.forEach(function (g) {
      var y = Y(g.len);
      svg += '<rect x="' + (cx - 26) + '" y="' + (y - 2) + '" width="52" height="' + (g.n > 1 ? 6 : 4) + '" rx="2" fill="#7FE8B4" opacity="' + (0.5 + Math.min(g.n, 3) * 0.15) + '"/>';
      svg += '<text x="' + (cx + 31) + '" y="' + (y + 3.5) + '" fill="#9FB6B2" font-size="9.5">' + g.len + (g.n > 1 ? ' ×' + g.n : '') + '</text>';
    });
  });
  svg += '<text x="' + (W - 10) + '" y="' + (H - 8) + '" text-anchor="end" fill="#6E8A86" font-size="9.5">' +
    '纵轴：片段长度（对数刻度）· ' + (circular ? '环状' : '线性') + '</text>';
  svg += '</svg>';
  box.innerHTML = svg;
}

function setupGel(seq, res) {
  var sel1 = document.getElementById('gelEnzSelect');
  var sel2 = document.getElementById('digestEnz2');
  if (!sel1 || !sel2) return;
  var found = res.filter(function (r) { return r.hits.length > 0; });
  var cur1 = sel1.value, cur2 = sel2.value;
  sel1.innerHTML = '';
  sel2.innerHTML = '<option value="">（不做双酶切）</option>';
  found.forEach(function (r) {
    var o1 = document.createElement('option');
    o1.value = r.enz[0];
    o1.textContent = r.enz[0] + '（' + r.hits.length + ' 切点）';
    sel1.appendChild(o1);
    var o2 = document.createElement('option');
    o2.value = r.enz[0];
    o2.textContent = r.enz[0] + '（' + r.hits.length + ' 切点）';
    sel2.appendChild(o2);
  });
  if (!found.length) {
    var o = document.createElement('option');
    o.value = ''; o.textContent = '未检测到可切的酶';
    sel1.appendChild(o);
  }
  if (cur1) sel1.value = cur1;
  if (cur2) sel2.value = cur2;
  renderDigest(seq, res);
}

/* ---------------- 5. GenBank 导出 ---------------- */
function gbDate() {
  var m = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'];
  var d = new Date();
  return ('0' + d.getDate()).slice(-2) + '-' + m[d.getMonth()] + '-' + d.getFullYear();
}

function wrapQualifier(key, value) {
  return '                     /' + key + '="' + value + '"';
}

function buildGenBank(seq, opts) {
  opts = opts || {};
  var name = opts.name || 'seq_workshop';
  var circular = !!opts.circular;
  var L = seq.length;
  var lines = [];
  var locusName = (name.length > 16 ? name.slice(0, 16) : name);
  var pad = locusName + new Array(Math.max(1, 16 - locusName.length + 1)).join(' ');
  lines.push('LOCUS       ' + pad + String(L).padStart(11, ' ') + ' bp    DNA     ' + (circular ? 'circular' : 'linear') + '   SYN ' + gbDate());
  lines.push('DEFINITION  序列工坊导出的 DNA 序列 ' + name + '.');
  lines.push('ACCESSION   .');
  lines.push('VERSION     .');
  lines.push('KEYWORDS    .');
  lines.push('SOURCE      synthetic construct');
  lines.push('  ORGANISM  synthetic construct');
  lines.push('            other sequences; artificial sequences.');
  lines.push('COMMENT     由序列工坊（Seq Workshop）在浏览器本地生成；坐标 1-based。');
  lines.push('FEATURES             Location/Qualifiers');
  lines.push('     source          1..' + L);
  lines.push(wrapQualifier('organism', 'synthetic construct'));
  lines.push(wrapQualifier('mol_type', 'genomic DNA'));

  // ORF 注释
  var orfRows = [];
  for (var f = 0; f < 3; f++) {
    var pos = findORFs(seq, f);
    var neg = findORFs(reverseComplement(seq), f);
    pos.forEach(function (o) {
      orfRows.push({ start: o.start + 1, end: o.end, strand: '+', frame: f + 1, aa: o.prot.length, prot: o.prot });
    });
    neg.forEach(function (o) {
      var s = L - o.end + 1, e = L - o.start;
      orfRows.push({ start: s, end: e, strand: '-', frame: f + 1, aa: o.prot.length });
    });
  }
  orfRows.sort(function (a, b) { return a.start - b.start; });
  orfRows.slice(0, 40).forEach(function (o) {
    lines.push('     misc_feature    ' + (o.strand === '-' ? 'complement(' : '') + o.start + '..' + o.end + (o.strand === '-' ? ')' : ''));
    lines.push(wrapQualifier('label', 'ORF ' + o.strand + o.frame + ' (' + o.aa + ' aa)'));
  });

  // 引物注释
  if (window._lastPrimers && window._lastPrimers.length) {
    var p = window._lastPrimers[0];
    lines.push('     primer_bind     ' + p.fStart + '..' + (p.fStart + p.fLen - 1));
    lines.push(wrapQualifier('label', 'forward primer'));
    lines.push(wrapQualifier('note', p.fSeq));
    lines.push('     primer_bind     complement(' + p.rStart + '..' + p.rEnd + ')');
    lines.push(wrapQualifier('label', 'reverse primer'));
    lines.push(wrapQualifier('note', p.rSeq));
  }

  // 常用酶切位点
  var siteCount = 0;
  for (var i = 0; i < ENZYMES.length && siteCount < 20; i++) {
    var hits = scanEnzyme(seq, ENZYMES[i]);
    if (hits.length && ENZYMES[i][0].length >= 5) {
      hits.slice(0, 2).forEach(function (h) {
        if (siteCount >= 20) return;
        siteCount++;
        lines.push('     misc_feature    ' + (h.start + 1) + '..' + (h.start + ENZYMES[i][1].length));
        lines.push(wrapQualifier('label', ENZYMES[i][0] + ' site'));
      });
    }
  }

  // 用户注释
  (mapFeatures || []).forEach(function (f) {
    var s0 = Math.min(f.start, f.end), e0 = Math.max(f.start, f.end);
    lines.push('     misc_feature    ' + (f.strand === '-' ? 'complement(' : '') + s0 + '..' + e0 + (f.strand === '-' ? ')' : ''));
    lines.push(wrapQualifier('label', f.label || 'annotation'));
    if (f.note) lines.push(wrapQualifier('note', f.note));
  });

  lines.push('ORIGIN');
  var lower = seq.toLowerCase();
  for (var p2 = 0; p2 < lower.length; p2 += 60) {
    var chunk = lower.substr(p2, 60);
    var blocks = [];
    for (var q = 0; q < chunk.length; q += 10) blocks.push(chunk.substr(q, 10));
    lines.push(String(p2 + 1).padStart(9, ' ') + ' ' + blocks.join(' '));
  }
  lines.push('//');
  return lines.join('\n') + '\n';
}

function exportGenBank() {
  var seq = window._currentSeq;
  if (!seq) { alert('请先输入并分析一段序列'); return; }
  var topo = document.getElementById('digestTopo');
  var gb = buildGenBank(seq, { name: 'seq_workshop', circular: !!(topo && topo.value === 'circular') });
  var blob = new Blob([gb], { type: 'text/plain;charset=utf-8' });
  var a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = 'seq_workshop.gb';
  a.click();
}

function addGenBankButton() {
  var bar = document.getElementById('exportBar');
  if (!bar || document.getElementById('btnGenBank')) return;
  var b = document.createElement('button');
  b.className = 'btn btn-ghost btn-sm';
  b.id = 'btnGenBank';
  b.textContent = '导出 GenBank';
  b.onclick = exportGenBank;
  bar.appendChild(b);
}

/* ---------------- 6. 局部比对 + 点阵图 ---------------- */
function swAlign(a, b, sc, mm, gap) {
  var n = a.length, m = b.length;
  var F = [], T = [];
  for (var i = 0; i <= n; i++) { F.push(new Array(m + 1).fill(0)); T.push(new Array(m + 1).fill(-1)); }
  var best = 0, bi = 0, bj = 0;
  for (var x = 1; x <= n; x++) {
    for (var y = 1; y <= m; y++) {
      var s = (a.charAt(x - 1) === b.charAt(y - 1)) ? sc : mm;
      var diag = F[x - 1][y - 1] + s;
      var up = F[x - 1][y] + gap;
      var left = F[x][y - 1] + gap;
      var v = Math.max(0, diag, up, left);
      F[x][y] = v;
      T[x][y] = (v === 0) ? -1 : (v === diag ? 0 : (v === up ? 1 : 2));
      if (v > best) { best = v; bi = x; bj = y; }
    }
  }
  var aa = '', bb = '', mid = '';
  var ci = bi, cj = bj;
  while (ci > 0 && cj > 0 && T[ci][cj] !== -1) {
    var t = T[ci][cj];
    if (t === 0) {
      var ca = a.charAt(ci - 1), cb = b.charAt(cj - 1);
      aa = ca + aa; bb = cb + bb; mid = (ca === cb ? '|' : '·') + mid; ci--; cj--;
    } else if (t === 1) { aa = a.charAt(ci - 1) + aa; bb = '-' + bb; mid = ' ' + mid; ci--; }
    else { aa = '-' + aa; bb = b.charAt(cj - 1) + bb; mid = ' ' + mid; cj--; }
  }
  var match = 0;
  for (var z = 0; z < mid.length; z++) if (mid.charAt(z) === '|') match++;
  return { a: aa, b: bb, mid: mid, score: best, match: match, len: mid.length, startA: ci, startB: cj };
}

function drawDotPlot(a, b) {
  var cv = document.getElementById('dotPlot');
  if (!cv) return;
  var ctx = cv.getContext('2d');
  var W = cv.width, H = cv.height, pad = 46;
  ctx.clearRect(0, 0, W, H);
  var dark = document.documentElement.getAttribute('data-theme') === 'dark';
  ctx.fillStyle = dark ? '#0E1615' : '#FBFDFD';
  ctx.fillRect(0, 0, W, H);
  var pw = W - pad - 16, ph = H - pad - 16;
  ctx.strokeStyle = dark ? '#26322F' : '#D7E3E1';
  ctx.lineWidth = 1;
  ctx.strokeRect(pad, 12, pw, ph);
  var k = Math.max(4, Math.min(10, Math.round(Math.min(a.length, b.length) / 40)));
  var n = 0;
  for (var i = 0; i + k <= a.length; i++) {
    var kmer = a.substr(i, k);
    for (var j = 0; j + k <= b.length; j++) {
      if (b.substr(j, k) === kmer) {
        var x = pad + (i / a.length) * pw;
        var y = 12 + ph - (j / b.length) * ph;
        ctx.fillStyle = dark ? 'rgba(127,232,180,.55)' : 'rgba(14,124,123,.45)';
        ctx.fillRect(x, y, 1.6, 1.6);
        n++;
        if (n > 60000) break;
      }
    }
    if (n > 60000) break;
  }
  ctx.fillStyle = dark ? '#8FA5A1' : '#3E5C59';
  ctx.font = '11px "JetBrains Mono", monospace';
  ctx.fillText('序列 1（' + a.length + ' bp）', pad, H - 14);
  ctx.save();
  ctx.translate(14, 12 + ph);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText('序列 2（' + b.length + ' bp）', 0, 0);
  ctx.restore();
  ctx.fillText('k=' + k + ' 精确匹配 · 对角线＝同源片段', pad + 6, 24);
}

/* ---------------- 初始化 ---------------- */
(function initV3() {
  var tbl = document.getElementById('codonTableSel');
  if (tbl) tbl.addEventListener('change', function () {
    activeCodonTable = this.value;
    if (window._currentSeq) analyze();
  });

  var primerBtn = document.getElementById('primerDesign');
  if (primerBtn) primerBtn.onclick = designPrimersV2;

  var topo = document.getElementById('digestTopo');
  if (topo) topo.addEventListener('change', function () { if (lastSeq && lastEnz) renderDigest(lastSeq, lastEnz); });
  var e2 = document.getElementById('digestEnz2');
  if (e2) e2.addEventListener('change', function () { if (lastSeq && lastEnz) renderDigest(lastSeq, lastEnz); });

  // 用克隆节点替换，清掉旧版绑定
  var oldGel = document.getElementById('gelEnzSelect');
  if (oldGel) {
    var g2 = oldGel.cloneNode(true);
    oldGel.parentNode.replaceChild(g2, oldGel);
    g2.addEventListener('change', function () { if (lastSeq && lastEnz) renderDigest(lastSeq, lastEnz); });
  }

  var oldAlign = document.getElementById('btnAlign');
  if (oldAlign) {
    var a2 = oldAlign.cloneNode(true);
    oldAlign.parentNode.replaceChild(a2, oldAlign);
    a2.addEventListener('click', function () {
      var out = document.getElementById('alignResult');
      var p1 = parseInput(document.getElementById('seq-input').value);
      var p2 = parseInput(document.getElementById('alignInput2').value);
      if(p1.error||p2.error){out.textContent=p1.error||p2.error;return;}
      if (!p1.seq || !p2.seq) {
        out.innerHTML = '<p class="empty-note">需要两条序列才能比对：请确认上方主序列与第二条序列都已填写。</p>';
        return;
      }
      if (p1.seq.length > 1200 || p2.seq.length > 1200) {
        out.innerHTML = '<p class="empty-note">为保证浏览器流畅，本工具比对长度上限为 1200 bp，当前序列超出。</p>';
        return;
      }
      var mode = document.getElementById('alignMode').value;
      var sc = +document.getElementById('pMatch').value;
      var mm = +document.getElementById('pMismatch').value;
      var gp = +document.getElementById('pGap').value;
      var r = mode === 'local' ? swAlign(p1.seq, p2.seq, sc, mm, gp) : nwAlign(p1.seq, p2.seq, sc, mm, gp);
      var identity = r.len ? (r.match / r.len * 100) : 0;
      var gaps = (r.a.match(/-/g) || []).length;
      var h = '<div class="align-stats">' +
        '<span>算法 <b>' + (mode === 'local' ? 'Smith–Waterman 局部' : 'Needleman–Wunsch 全局') + '</b></span>' +
        '<span>比对得分 <b>' + r.score + '</b></span>' +
        '<span>一致位点 <b>' + r.match + '</b></span>' +
        '<span>空位数 <b>' + gaps + '</b></span>' +
        '<span>一致性 <b>' + identity.toFixed(1) + '%</b></span></div>';
      h += '<div class="align-view"><div>' + colorDNA(r.a) + '</div><div>' + colorMid(r.mid) + '</div><div>' + colorDNA(r.b) + '</div></div>';
      out.innerHTML = h;
      drawDotPlot(p1.seq, p2.seq);
    });
  }

  addGenBankButton();
})();


/* =========================================================
   v3.1 界面交互：滚动高亮 / 轻提示 / 分析按钮状态
   ========================================================= */
