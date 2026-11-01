var WS_KEY = 'seqw_next_workspace_v4';
var wsTimer = 0;

function scheduleWorkspaceSave() {
  clearTimeout(wsTimer);
  wsTimer = setTimeout(saveWorkspaceNow, 400);
}

function saveWorkspaceNow() {
  if(window.SeqModuleRouter && !window.SeqModuleRouter.hasSequence())return;
  var seq = window._currentSeq || '';
  if (!seq) return;
  try {
    var data=SeqWorkspaceFormat.create(seq,mapFeatures || [],{codonTable:activeCodonTable,topology:seqMapState.circular?'circular':'linear'});
    localStorage.setItem(WS_KEY, JSON.stringify(data));
    document.getElementById('workspaceStatus').textContent='已在当前浏览器保存 '+seq.length+' bp 与 '+data.features.length+' 条注释；可导出备份换设备恢复。';
  } catch (e) { document.getElementById('workspaceStatus').textContent='浏览器本地保存不可用，请导出工作区备份以防丢失。'; }
}

function loadWorkspace() {
  try {
    var raw = localStorage.getItem(WS_KEY);
    if (!raw) return null;
    var clean=SeqWorkspaceFormat.validate(JSON.parse(raw));
    return {seq:clean.sequence,features:clean.features,primers:[],mutation:null,circular:clean.settings.topology==='circular',codonTable:clean.settings.codonTable};
  } catch (e) { return null; }
}

function clearWorkspace() {
  try { localStorage.removeItem(WS_KEY); } catch (e) {}
}

function restoreWorkspace() {
  var ws = loadWorkspace();
  if (!ws || !ws.seq) return;
  document.getElementById('seq-input').value = ws.seq;
  activeCodonTable=ws.codonTable;
  document.getElementById('codonTableSel').value=activeCodonTable;
  window._lastAnalyzedSeq=null;
  window._pendingFeatures = ws.features || [];
  analyze();
  window._lastPrimers = ws.primers || [];
  window._seqMapMutation = ws.mutation || null;
  if (ws.circular) {
    seqMapState.circular = true;
    var lin = document.getElementById('mapLinear');
    var cir = document.getElementById('mapCircular');
    if (lin) lin.classList.remove('active');
    if (cir) cir.classList.add('active');
  }
  renderSeqMap();
  renderMapFeatures();
  if (typeof showToast === 'function') showToast('已恢复上次的工作区');
}

function initWorkspace() {
  var bar = document.getElementById('workspaceBar');
  if (!bar) return;
  var ws = loadWorkspace();
  if (ws && ws.seq) {
    bar.hidden = false;
    document.getElementById('wsRestore').addEventListener('click', function () {
      restoreWorkspace();
      bar.hidden = true;
    });
    document.getElementById('wsClear').addEventListener('click', function () {
      clearWorkspace();
      bar.hidden = true;
    });
  }
  var upBtn = document.getElementById('btnUpload');
  var upInput = document.getElementById('uploadFile');
  if (upBtn && upInput) {
    upBtn.addEventListener('click', function () { upInput.click(); });
    upInput.addEventListener('change', function () {
      var f = this.files && this.files[0];
      if (!f) return;
      var fr = new FileReader();
      fr.onload = function () {
        var text = String(fr.result || '');
        if (/^\s*LOCUS/m.test(text)) {
          importGbText(text);
        } else {
          document.getElementById('seq-input').value = text;
          analyze();
          if (typeof showToast === 'function') showToast('已加载文件：' + f.name);
        }
      };
      fr.readAsText(f);
      this.value = '';
    });
  }
  window.addEventListener('beforeunload', saveWorkspaceNow);
}

/* 让 v3 功能在 analyze() 之后也能重新初始化按钮 */
var _v3Ready = true;
initWorkspace();
