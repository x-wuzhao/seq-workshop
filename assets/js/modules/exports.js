function exportCSV(filename, rows){
  var csv = rows.map(function(r){return r.map(function(c){return '"'+String(c).replace(/"/g,'""')+'"';}).join(',');}).join('\n');
  var blob = new Blob(['\ufeff'+csv], {type:'text/csv;charset=utf-8'});
  var a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = filename;
  a.click();
}

function exportFasta(seq, name){
  var fasta = '>'+name+'\n';
  for(var i=0;i<seq.length;i+=60) fasta += seq.substring(i,i+60)+'\n';
  var blob = new Blob([fasta], {type:'text/plain;charset=utf-8'});
  var a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = name+'.fasta';
  a.click();
}

function addExportButtons(){
  // 在统计卡片区域添加导出按钮
  var statsPanel = document.querySelector('#nav-complement');
  if(statsPanel && !document.getElementById('exportBar')){
    var bar = document.createElement('div');
    bar.className = 'export-bar';
    bar.id = 'exportBar';
    bar.innerHTML = '<button class="btn btn-ghost btn-sm" onclick="exportFasta(window._currentSeq,\'序列工坊_输入序列\')">导出 FASTA</button>' +
      '<button class="btn btn-ghost btn-sm" onclick="exportResultsCSV()">导出分析结果 CSV</button>' +
      '<button class="btn btn-ghost btn-sm" onclick="exportProteinCSV()">导出蛋白序列 CSV</button>';
    statsPanel.appendChild(bar);
  }
}

function exportResultsCSV(){
  var seq = window._currentSeq;
  if(!seq) return;
  var c = baseCounts(seq);
  var rows = [
    ['指标','数值'],
    ['序列长度(bp)', seq.length],
    ['A含量', c.A],
    ['T含量', c.T],
    ['G含量', c.G],
    ['C含量', c.C],
    ['N未知碱基', c.N],
    ['GC含量(%)', ((c.G+c.C)/seq.length*100).toFixed(2)],
    ['AT含量(%)', ((c.A+c.T)/seq.length*100).toFixed(2)],
    ['N含量(%)', (c.N/seq.length*100).toFixed(2)],
    ['解链温度Tm(°C)', c.N?'未知（含 N）':calcTm(c, seq.length).toFixed(1)],
    ['单链分子量(kDa)', c.N?'未知（含 N）':(calcMW(c)/1000).toFixed(2)]
  ];
  exportCSV('序列工坊_分析结果.csv', rows);
}

function exportProteinCSV(){
  var seq = window._currentSeq;
  if(!seq) return;
  var rows = [['读码框','方向','ORF长度(aa)','蛋白序列']];
  for(var f=0;f<3;f++){
    var orfs = findORFs(seq,f);
    for(var i=0;i<orfs.length;i++){
      rows.push(['Frame '+(f+1),'正向',orfs[i].prot.length,orfs[i].prot]);
    }
    var rc = reverseComplement(seq);
    var rorfs = findORFs(rc,f);
    for(var j=0;j<rorfs.length;j++){
      rows.push(['Frame '+(f+1),'反向',rorfs[j].prot.length,rorfs[j].prot]);
    }
  }
  exportCSV('序列工坊_蛋白序列.csv', rows);
}

/* ============================================================
   八、示例序列菜单
   ============================================================ */
