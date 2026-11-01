function renderEnzMap(seq,res){
  var W=620,H=118,x0=46,x1=W-22,y=56;
  function X(pos){return x0+(pos/seq.length)*(x1-x0);}
  var svg='<svg viewBox="0 0 '+W+' '+H+'" xmlns="http://www.w3.org/2000/svg">';
  var ts=[0,.25,.5,.75,1],nHits=0;
  for(var i=0;i<ts.length;i++){
    var tpos=Math.round(ts[i]*seq.length),tx=X(tpos);
    svg+='<line x1="'+tx+'" y1="'+y+'" x2="'+tx+'" y2="'+(y+5)+'" stroke="#3E5C59"/>';
    svg+='<text x="'+tx+'" y="'+(y+19)+'" font-size="10" text-anchor="middle" fill="#3E5C59">'+tpos+'</text>';
  }
  svg+='<line x1="'+x0+'" y1="'+y+'" x2="'+x1+'" y2="'+y+'" stroke="#0E7C7B" stroke-width="2"/>';
  for(var r=0;r<res.length;r++){
    for(var h=0;h<res[r].hits.length;h++){
      var hit=res[r].hits[h],cx=X(hit.cut);
      svg+='<line x1="'+cx+'" y1="'+(y-11)+'" x2="'+cx+'" y2="'+(y+11)+'" stroke="#C43E3E"><title>'+res[r].enz[0]+' · 切点 '+hit.cut+' bp</title></line>';
      nHits++;
    }
  }
  svg+='</svg>';
  document.getElementById('enzMap').innerHTML=nHits?svg:'<p class="empty-note">未检测到任何酶切位点。</p>';
}

// 酶切结果表格
function renderEnzTable(seq,res,filter){
  var box=document.getElementById('enzTableBox');
  var rows=(filter==='found')?res.filter(function(r){return r.hits.length>0;}):res;
  if(filter==='found'&&!rows.length){box.innerHTML='<p class="empty-note">未检测到常见酶的切点。</p>';return;}
  var h='<table class="enz-table"><thead><tr>'+
    '<th style="width:11%">酶</th><th style="width:14%">识别序列</th><th style="width:9%">切点数</th>'+
    '<th style="width:36%">切点位置 (bp)</th><th style="width:30%">片段长度 (bp)</th></tr></thead><tbody>';
  for(var i=0;i<rows.length;i++){
    var r=rows[i],site=r.enz[1],c1=r.enz[2];
    var siteDisp=site.substr(0,c1)+'<span class="cutmark">^</span>'+site.substr(c1);
    var posHtml='',fragHtml='';
    if(!r.hits.length){posHtml='—';fragHtml='—';}
    else{
      var prev=0;
      for(var k=0;k<r.hits.length;k++){
        posHtml+='<span class="pos">'+r.hits[k].cut+'</span>';
        fragHtml+='<span class="pos">'+(r.hits[k].cut-prev)+'</span>';
        prev=r.hits[k].cut;
      }
      fragHtml+='<span class="pos">'+(seq.length-prev)+'</span>';
    }
    h+='<tr><td><b>'+r.enz[0]+'</b></td><td class="site">'+siteDisp+'</td><td>'+r.hits.length+
      '</td><td>'+posHtml+'</td><td>'+fragHtml+'</td></tr>';
  }
  box.innerHTML=h+'</tbody></table>';
}

/* ===== 虚拟凝胶电泳 ===== */
// 计算某酶切后的片段长度
function enzymeFragments(seqLen,hits){
  var cuts=[];
  for(var i=0;i<hits.length;i++)cuts.push(hits[i].cut);
  cuts.sort(function(a,b){return a-b;});
  var frags=[],prev=0;
  for(var j=0;j<cuts.length;j++){frags.push(cuts[j]-prev);prev=cuts[j];}
  frags.push(seqLen-prev);
  return frags;
}
var GEL_MARKER=[500,400,300,200,150,100,50];  // bp ladder

function renderGel(seq,res,enzName){
  var target=null;
  for(var i=0;i<res.length;i++){if(res[i].enz[0]===enzName){target=res[i];break;}}
  if(!target)return;
  var frags=enzymeFragments(seq.length,target.hits);

  var W=480,H=300,xMarker=140,xSample=330,laneW=92;
  var yTop=50,yBot=268;
  var bpMax=GEL_MARKER[0],bpMin=GEL_MARKER[GEL_MARKER.length-1];
  function Y(bp){
    var b=bp;if(b>bpMax)b=bpMax;if(b<bpMin)b=bpMin;
    return yTop+(Math.log(bpMax)-Math.log(b))/(Math.log(bpMax)-Math.log(bpMin))*(yBot-yTop);
  }

  var h='<svg viewBox="0 0 '+W+' '+H+'" xmlns="http://www.w3.org/2000/svg">';
  // 点样孔
  h+='<rect x="'+(xMarker-laneW/2+8)+'" y="32" width="'+(laneW-16)+'" height="10" rx="2" fill="#0d1a16" stroke="rgba(127,232,180,.25)"/>';
  h+='<rect x="'+(xSample-laneW/2+8)+'" y="32" width="'+(laneW-16)+'" height="10" rx="2" fill="#0d1a16" stroke="rgba(127,232,180,.25)"/>';
  // lane 边框线
  var edges=[xMarker,xSample];
  for(var e=0;e<2;e++){
    var xc=edges[e];
    h+='<line class="gel-lane-line" x1="'+(xc-laneW/2)+'" y1="'+yTop+'" x2="'+(xc-laneW/2)+'" y2="'+yBot+'"/>';
    h+='<line class="gel-lane-line" x1="'+(xc+laneW/2)+'" y1="'+yTop+'" x2="'+(xc+laneW/2)+'" y2="'+yBot+'"/>';
  }
  // lane 标签
  h+='<text class="gel-lane-label" x="'+xMarker+'" y="22" text-anchor="middle">Marker</text>';
  h+='<text class="gel-lane-label" x="'+xSample+'" y="22" text-anchor="middle">'+enzName+' 酶切</text>';
  // Marker 条带 + 长度标签
  for(var m=0;m<GEL_MARKER.length;m++){
    var bp=GEL_MARKER[m],y=Y(bp);
    h+='<rect class="gel-band bright" x="'+(xMarker-laneW/2+10)+'" y="'+(y-2.5)+'" width="'+(laneW-20)+'" height="5" rx="2"/>';
    h+='<text class="gel-size-tag" x="'+(xMarker-laneW/2-6)+'" y="'+(y+3.5)+'" text-anchor="end">'+bp+'</text>';
  }
  // 样本条带（相同长度合并，条带数越多越亮；标签自动避让）
  var seen={};
  for(var f=0;f<frags.length;f++){var fl=frags[f];seen[fl]=(seen[fl]||0)+1;}
  var bands=[];
  for(var key in seen){bands.push({size:+key,y:Y(+key),n:seen[key]});}
  bands.sort(function(a,b){return a.y-b.y;});
  var tagX=xSample+laneW/2+6,prevTagY=-99;
  for(var bi=0;bi<bands.length;bi++){
    var bd=bands[bi],op=0.65+Math.min(bd.n,3)*0.12;
    h+='<rect class="gel-band" x="'+(xSample-laneW/2+10)+'" y="'+(bd.y-2.5)+'" width="'+(laneW-20)+'" height="5" rx="2" opacity="'+op+'"/>';
    var ty=bd.y+3.5;
    if(ty-prevTagY<13)ty=prevTagY+13;
    prevTagY=ty;
    if(Math.abs(ty-(bd.y+3.5))>1.5){
      h+='<line x1="'+(tagX-2)+'" y1="'+bd.y+'" x2="'+(tagX-2)+'" y2="'+(ty-4)+'" stroke="rgba(159,198,194,.55)" stroke-width=".8"/>';
    }
    h+='<text class="gel-size-tag" x="'+tagX+'" y="'+ty+'" text-anchor="start">'+bd.size+'</text>';
  }
  h+='</svg>';
  document.getElementById('gelBox').innerHTML=h;
}

// 初始化凝胶下拉并渲染
function setupGel(seq,res){
  var sel=document.getElementById('gelEnzSelect');
  sel.innerHTML='';
  var found=res.filter(function(r){return r.hits.length>0;});
  for(var i=0;i<found.length;i++){
    var opt=document.createElement('option');
    opt.value=found[i].enz[0];opt.textContent=found[i].enz[0]+'（'+found[i].hits.length+' 切点）';
    sel.appendChild(opt);
  }
  if(found.length)renderGel(seq,res,found[0].enz[0]);
}

// 密码子使用频率（按氨基酸分组，热力底色，频率以‰计）
function renderCodonGrid(seq,gBest){
  var CODON=codonMap();
  var source=gBest&&gBest.frameSeq?gBest.frameSeq:seq;
  var region=gBest?source.substring(gBest.start,gBest.end):source.substring(0,source.length-(source.length%3));
  var counts={},total=0;
  for(var i=0;i+3<=region.length;i+=3){
    var c=region.substr(i,3);
    if(c.length===3&&/^[ATCG]{3}$/.test(c)){counts[c]=(counts[c]||0)+1;total++;}
  }
  var h='';
  for(var a=0;a<AA_ORDER.length;a++){
    var aa=AA_ORDER[a];
    h+='<div class="codon-aahead">'+aa+' '+AA_NAME[aa]+'</div>';
    for(var key in CODON){
      if(CODON[key]===aa){
        var n=counts[key]||0,per=total?(n/total*1000):0;
        var op=Math.min(.45,per/1000*1.5);
        var codonBg=isDarkTheme()?'rgba(43,163,161,':'rgba(14,124,123,';
        h+='<div class="codon-cell" style="background:'+codonBg+op.toFixed(2)+';border:1px solid '+(isDarkTheme()?'rgba(43,163,161,.2)':'rgba(14,124,123,.15)')+'">'+
          '<span class="c">'+key+'</span><span class="a">'+n+'</span><span class="f">'+per.toFixed(1)+'‰</span></div>';
      }
    }
  }
  document.getElementById('codonGrid').innerHTML=h;
}

function renderCharts(seq,avgGC,aaCounts){
  lastChartArgs={seq:seq,avgGC:avgGC,aaCounts:aaCounts};
  if(typeof echarts==='undefined'){
    var g=document.getElementById('chartGC'),a=document.getElementById('chartAA');
    g.innerHTML='<p style="padding:30px;color:#3E5C59;font-size:13px">图表库未加载；文字与序列结果仍可使用。请确认 assets/echarts.min.js 已随应用提供。</p>';
    a.innerHTML='<p style="padding:30px;color:#3E5C59;font-size:13px">图表库未加载；文字与序列结果仍可使用。请确认 assets/echarts.min.js 已随应用提供。</p>';
    return;
  }
  if(!document.getElementById('nav-charts').hidden)doChartsRender(seq,avgGC,aaCounts);
}

function doChartsRender(seq,avgGC,aaCounts){
  if(document.getElementById('nav-charts').hidden)return;
  var gcBox=document.getElementById('chartGC'),aaBox=document.getElementById('chartAA');
  if(charts.gc){charts.gc.dispose();}
  if(charts.aa){charts.aa.dispose();}
  charts.gc=echarts.init(gcBox);
  charts.aa=echarts.init(aaBox);
  var ac=chartAxisColor(),sc=chartSplitColor();

  var w=gcWindows(seq),lineData=[];
  for(var i=0;i<w.pts.length;i++){lineData.push([w.pts[i].pos,w.pts[i].gc]);}
  charts.gc.setOption({
    animation:false,
    tooltip:{trigger:'axis',triggerOn:'click',renderMode:'richText',confine:true,
      formatter:function(p){return '位置 '+p[0].data[0]+' bp'+String.fromCharCode(10)+'GC 含量 '+p[0].data[1]+'%';}},
    grid:{left:8,right:24,top:18,bottom:6,containLabel:true},
    xAxis:{type:'value',name:'位置 (bp)',min:0,nameLocation:'middle',nameGap:28,
      axisLabel:{color:ac},axisLine:{lineStyle:{color:ac}},nameTextStyle:{color:ac},
      splitLine:{lineStyle:{color:sc}}},
    yAxis:{type:'value',min:0,max:100,
      axisLabel:{color:ac},axisLine:{lineStyle:{color:ac}},
      splitLine:{lineStyle:{color:sc}}},
    series:[{
      type:'line',smooth:true,symbol:'circle',symbolSize:4,data:lineData,
      areaStyle:{opacity:.08},
      markLine:{silent:true,symbol:'none',lineStyle:{type:'dashed'},
        data:[{yAxis:+avgGC.toFixed(1),label:{formatter:'平均 '+avgGC.toFixed(1)+'%',position:'insideStartTop',color:ac}}]}
    }]
  });

  var rows=[];
  for(var k=0;k<AA_ORDER.length;k++){
    if(aaCounts[AA_ORDER[k]]>0)rows.push({aa:AA_ORDER[k],n:aaCounts[AA_ORDER[k]]});
  }
  rows.sort(function(a,b){return a.n-b.n;});
  var cats=[],vals=[];
  for(var r=0;r<rows.length;r++){cats.push(rows[r].aa);vals.push(rows[r].n);}
  charts.aa.setOption({
    animation:false,
    tooltip:{trigger:'axis',axisPointer:{type:'shadow'},triggerOn:'click',renderMode:'richText',confine:true,
      formatter:function(p){var aa=p[0].name;return aa+' '+AA_NAME[aa]+String.fromCharCode(10)+p[0].value+' 个';}},
    grid:{left:8,right:30,top:14,bottom:6,containLabel:true},
    xAxis:{type:'value',minInterval:1,
      axisLabel:{color:ac},axisLine:{lineStyle:{color:ac}},
      splitLine:{lineStyle:{color:sc}}},
    yAxis:{type:'category',data:cats,axisLabel:{fontSize:12,fontWeight:600,interval:0,color:ac},
      axisLine:{lineStyle:{color:ac}}},
    series:[{type:'bar',data:vals,barMaxWidth:15,
      label:{show:true,position:'right',formatter:'{c}',color:ac}}]
  });
}

window.addEventListener('resize',function(){
  if(charts.gc)charts.gc.resize();
  if(charts.aa)charts.aa.resize();
});

// 展开 / 收起
document.addEventListener('click',function(e){
  var b=e.target;
  if(b&&b.classList&&b.classList.contains('expand-link')){
    var t=document.getElementById(b.getAttribute('data-target'));
    var isC=t.classList.contains('collapsed');
    t.classList.toggle('collapsed');
    b.textContent=isC?'收起':'展开全部';
  }
});
