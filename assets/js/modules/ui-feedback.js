function showToast(msg){
  var t=document.getElementById('toast');
  if(!t){
    t=document.createElement('div');
    t.id='toast';t.className='toast';t.setAttribute('role','status');t.setAttribute('aria-live','polite');
    document.body.appendChild(t);
  }
  t.textContent=msg;
  t.classList.add('show');
  clearTimeout(t._hide);
  t._hide=setTimeout(function(){t.classList.remove('show');},1900);
}

// Task navigation is owned by module-router.js.

(function initAnalyzeButton(){
  var btn=document.getElementById('btnAnalyze');
  if(!btn||typeof analyze!=='function')return;
  var fresh=btn.cloneNode(true);
  btn.parentNode.replaceChild(fresh,btn);
  fresh.addEventListener('click',function(){
    var self=this;
    if(self.classList.contains('is-loading'))return;
    self.classList.add('is-loading');
    var label=self.textContent;
    self.textContent='分析中…';
    setTimeout(function(){
      try{
        analyze();
        var resultsEl=document.getElementById('results');
        if(resultsEl&&!resultsEl.hidden)showToast('分析完成');
      }finally{
        self.classList.remove('is-loading');
        self.textContent=label;
      }
    },30);
  });
})();

(function initExportToast(){
  ['exportCSV','exportFasta','exportGenBank'].forEach(function(fn){
    var original=window[fn];
    if(typeof original!=='function')return;
    window[fn]=function(){
      var r=original.apply(this,arguments);
      if(fn==='exportCSV')showToast('已导出 CSV');
      else if(fn==='exportFasta')showToast('已导出 FASTA');
      else showToast('已导出 GenBank（.gb）');
      return r;
    };
  });
})();


/* =========================================================
   十四、序列地图（横向注释轨道）
   ========================================================= */
