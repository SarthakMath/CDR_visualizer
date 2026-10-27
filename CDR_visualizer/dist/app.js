import {preset,exactSolution,isTransient,velocity,physicalDiffusion,addedDiffusion} from './fem.js';
import {meshPlot,surfacePlot,crossPlot,matrixPlot,fieldValues,short} from './plots.js';
const $=id=>document.getElementById(id);
let state=preset(),result=null,solvedState=null,requestId=0,timer,toastTimer,selected=null,editor='matrix';
let camera={yaw:-.66,elevation:.67},view='surface',field='raw',showMesh=true,slice=.5,points=[],meshPoints=[],pending=false;
let trajectory=null,frameIndex=0,selectedTime=0,playing=false,playRequest=0,lastPlayback=0;
const worker=new Worker(new URL('./worker.js',import.meta.url),{type:'module'});
const methodLabel=()=>state.stabilization==='supg'?`SUPG × ${state.supgScale}`:state.stabilization==='none'?'no stabilization':`isotropic δ = ${short(state.diffusion)}`;
const format=v=>Math.abs(v)<1e-12?'0.000':Math.abs(v)<.001||Math.abs(v)>9999?v.toExponential(2):v.toFixed(3);
function toast(message){$('toast').textContent=message;$('toast').classList.remove('hidden');clearTimeout(toastTimer);toastTimer=setTimeout(()=>$('toast').classList.add('hidden'),4000);}
function failure(message){$('errorBanner').textContent=message+' The plots show the last successful solve.';$('errorBanner').classList.remove('hidden');$('solveState').textContent='Check inputs';document.querySelector('.plot-card').classList.add('stale');$('export').disabled=true;pending=false;}
function schedule(){stopPlayback();clearTimeout(timer);pending=true;$('errorBanner').classList.add('hidden');$('solveState').textContent=isTransient(state)?'Building trajectory…':'Solving…';$('export').disabled=true;for(const id of ['playTime','previousTime','nextTime','rewindTime','timeSlider'])$(id).disabled=true;document.querySelector('.plot-card').classList.add('stale');requestId++;const id=requestId;timer=setTimeout(()=>worker.postMessage({id,state:structuredClone(state)}),90);}
worker.onmessage=({data})=>{
  if(data.id!==requestId)return;
  if(data.progress!==undefined){$('solveState').textContent=`Time steps · ${Math.round(data.progress*100)}%`;return;}
  if(data.error){failure(data.error);return;}
  pending=false;solvedState=structuredClone(state);$('solveState').textContent='✓ Solved';$('export').disabled=false;document.querySelector('.plot-card').classList.remove('stale');
  if(data.result.frames){trajectory=data.result;const index=trajectory.frames.reduce((best,f,i)=>Math.abs(f.time-selectedTime)<Math.abs(trajectory.frames[best].time-selectedTime)?i:best,0);selectFrame(index);}
  else{trajectory=null;result=data.result;}
  if(field==='reference'&&!result.reference)setField('raw');
  updateSummary();updateEditor();render();
};
worker.onerror=()=>failure('The solver could not start. Open the app through its web server (npm start), then reload.');

function stopPlayback(){playing=false;cancelAnimationFrame(playRequest);if($('playTime')){$('playTime').textContent='▶ Play';$('playTime').setAttribute('aria-label','Play time evolution');}}
function selectFrame(index){
  if(!trajectory)return;frameIndex=Math.max(0,Math.min(trajectory.frames.length-1,index));
  const f=trajectory.frames[frameIndex];selectedTime=f.time;
  result={...trajectory,...f,rawB:Float64Array.from(trajectory.mesh.interior,id=>f.raw[id]),stableB:Float64Array.from(trajectory.mesh.interior,id=>f.stable[id]),
    timeMatrix:frameIndex===trajectory.steps?trajectory.lastTimeMatrix:trajectory.timeMatrix,
    stableTimeMatrix:frameIndex===trajectory.steps?trajectory.lastStableTimeMatrix:trajectory.stableTimeMatrix};
  $('timeSlider').max=trajectory.steps;$('timeSlider').value=frameIndex;$('currentTime').textContent=`t = ${format(f.time)}`;
  $('frameCounter').textContent=`Frame ${frameIndex} / ${trajectory.steps}`;
  $('timelineNote').textContent=`${state.scheme==='be'?'Backward Euler':'Crank–Nicolson'} · Δt = ${Number(f.dt.toPrecision(5))} · comparison: ${state.massType==='lumped'?'lumped':'consistent'} M · ${methodLabel()}`;
  for(const id of ['playTime','previousTime','nextTime','rewindTime','timeSlider'])$(id).disabled=pending;
  $('previousTime').disabled=pending||frameIndex===0;$('nextTime').disabled=pending||frameIndex===trajectory.steps;
  updateSummary();updateEditor();render();
}
function playbackTick(timestamp){
  if(!playing||pending||!trajectory)return;
  const interval=1000/Number($('playbackSpeed').value);
  if(timestamp-lastPlayback>=interval){lastPlayback=timestamp;selectFrame(frameIndex+1);if(frameIndex===trajectory.steps){stopPlayback();return;}}
  playRequest=requestAnimationFrame(playbackTick);
}
$('playTime').addEventListener('click',()=>{
  if(!trajectory||pending)return;if(playing){stopPlayback();return;}
  if(frameIndex===trajectory.steps)selectFrame(0);playing=true;lastPlayback=0;$('playTime').textContent='Ⅱ Pause';$('playTime').setAttribute('aria-label','Pause time evolution');playRequest=requestAnimationFrame(playbackTick);
});
$('timeSlider').addEventListener('input',()=>{stopPlayback();if(!pending)selectFrame(Number($('timeSlider').value));});
for(const[id,change]of [['previousTime',-1],['nextTime',1],['rewindTime',0]])$(id).addEventListener('click',()=>{stopPlayback();if(!pending)selectFrame(change===0?0:frameIndex+change);});

function syncControls(){
  for(const key of ['mode','epsilon','reaction','betaX','betaY','source','amplitude','pattern','diffusion','stabilization','supgScale'])$(key).value=state[key];
  for(const side of ['left','right','bottom','top'])$(`bc-${side}`).value=state.boundaries[side];
  $('meshSize').value=state.n;$('meshSizeLabel').textContent=`${state.n} × ${state.n}`;
  $('diffusionSlider').max=Math.max(.15,state.diffusion);$('diffusionSlider').value=state.diffusion;$('diffusionMax').textContent=short(+$('diffusionSlider').max);
  $('supgControls').classList.toggle('hidden',state.stabilization!=='supg');$('isotropicControls').classList.toggle('hidden',state.stabilization!=='isotropic');
  const temporal=isTransient(state),convective=Math.hypot(...velocity({...state,betaX:1,betaY:1}))>0;
  $('velocityFields').classList.toggle('hidden',!convective);$('epsilon').disabled=['advection','advection-reaction','reaction'].includes(state.mode);
  for(const id of ['timeControls','timeline','timeDiagnostics'])$(id).classList.toggle('hidden',!temporal);
  $('matrixPickerWrap').classList.remove('hidden');
  for(const option of $('matrixPicker').options)option.disabled=!temporal&&['ML','step','comparisonMass'].includes(option.value);
  for(const key of ['initial','dt','duration','scheme','massType'])$(key).value=state[key];
  $('pdeLabel').textContent=[temporal?'uₜ':'',!['advection','advection-reaction','reaction'].includes(state.mode)?'− εΔu':'',convective?'+ β·∇u':'',(state.reaction!==0||state.mode.includes('reaction')||state.mode.includes('adr'))?'+ cu':''].filter(Boolean).join(' ').replace(/^\+ /,'')+' = f';
  $('discreteLabel').textContent=temporal?'Mḃ + Ab = f':'Ab = f';$('inspectorEquation').textContent=temporal?'Mḃ + Ab = f':'Ab = f';
  $('originalFieldName').textContent=$('originalLegendName').textContent=temporal?'Consistent baseline':'Original';
  $('comparisonFieldName').textContent=$('comparisonLegendName').textContent=temporal?'Comparison':state.stabilization==='supg'?'With SUPG':state.stabilization==='none'?'Comparison':'With diffusion';
  document.querySelectorAll('[data-preset]').forEach(b=>{b.classList.toggle('selected',b.dataset.preset===state.preset);b.setAttribute('aria-pressed',b.dataset.preset===state.preset);});
}
function updateSummary(){
  const r=result;$('tauSummary').textContent=`Element τ: ${r.tauRange.map(v=>short(v)).join(' – ')}${Math.hypot(...velocity(state))===0?' · no flow, SUPG is zero':''}`;$('nodeCount').textContent=r.mesh.nodes.length;$('elementCount').textContent=r.mesh.triangles.length;$('dofCount').textContent=r.N;
  $('systemSize').textContent=`${r.N} unknowns`;$('matrixShape').textContent=`${r.N} × ${r.N}`;
  const count=Object.keys(state.matrixEdits).length+Object.keys(state.rhsEdits).length+Object.keys(state.nodeEdits).length;$('editCount').textContent=`${count} override${count===1?'':'s'}`;
  $('timing').textContent=`${r.mesh.triangles.length} elements · ${r.ms.toFixed(0)} ms · ${trajectory?`${r.steps} time steps · reused LU`:'partial-pivot solve'}`;
  document.querySelector('[data-field="reference"]').disabled=!r.reference;
  $('exactLegend').classList.toggle('hidden',!r.reference);
  $('fieldBadge').textContent={poisson:'POISSON',transport:'CONVECTION–DIFFUSION',heat:'DIFFUSION',advection:'CONVECTION',parabolic:'PARABOLIC',reaction:'REACTION','diffusion-reaction':'DIFFUSION–REACTION','advection-reaction':'ADVECTION–REACTION',adr:'ADVECTION–DIFFUSION–REACTION','steady-adr':'STEADY ADR'}[state.mode];
  $('fieldSubtitle').textContent=field==='reference'?'Exact solution sampled at mesh nodes':field==='difference'?'Comparison minus original · shared coordinates':trajectory?`t = ${format(r.time)} · ${field==='stable'&&state.massType==='lumped'?'lumped mass':'consistent mass'} · ${state.scheme==='be'?'backward Euler':'Crank–Nicolson'}${field==='stable'?` · ${methodLabel()}`:' · δ = 0'}`:'Piecewise-linear finite element solution';
  let note;
  if(count)note=`${count} algebraic or boundary override${count===1?' is':'s are'} active. The original PDE reference no longer applies. Changes to mesh or triangulation clear these edits.`;
  else if(trajectory&&r.initialMismatch)note=`Initial and boundary data disagree at ${r.initialMismatch} boundary nodes. Prescribed boundary values take precedence; the two runs share the same initial interior field.`;
  else if(state.stabilization==='supg')note=`SUPG adds residual stabilization along the flow (τ = ${r.tauRange.map(v=>short(v)).join(' – ')}). Compare extrema and L² error: oscillations can remain. Switch to isotropic diffusion to compare spreading across the flow; refine the mesh to resolve the layer.`;
  else if(trajectory){
    if(r.initialMismatch)note=`Initial and boundary data disagree at ${r.initialMismatch} boundary nodes. Prescribed boundary data take precedence. Playback compares the same initial interior field.`;
    else if(state.reaction!==0)note=`Reaction c = ${state.reaction}: ${state.reaction>0?'decay':'growth'}. Compare the exact profile at the same time. Added diffusion spreads the field; reaction changes its amplitude. Mass lumping changes only the time derivative, while cM remains consistent.`;
    else if(state.mode==='heat')note='Step forward from the narrow initial pulse. Lumped mass with backward Euler preserves nonnegativity here; compare negative nodes and L² error. Positivity and accuracy are different questions.';
    else if(state.mode==='advection')note='Watch the front move with the flow. Lumping alone can still produce oscillations; compare with added diffusion and switch between Crank–Nicolson and backward Euler to see temporal damping.';
    else note='Physical diffusion spreads the pulse as convection carries it. Compare the effect of extra δK and mass lumping against the exact profile at the same time.';
  }
  else if(state.reaction!==0)note='The spatial operator is A = εK + C + cM. Reaction changes the balance of the steady solution; added δK changes only diffusion.';
  else if(state.preset==='layer'&&state.mode==='transport')note=state.diffusion>0?`Added diffusion δ = ${Number(state.diffusion.toPrecision(5))} brings the coordinate mesh Péclet number from ${short(r.pe)} to ${short(r.effectivePe)}. Compare oscillations and layer spreading in the cross-section.`:`The layer is thin relative to the mesh (Peₕ = ${short(r.pe)}). Try “Apply mesh-based value” to compare the original solution with artificial diffusion.`;
  else if(state.preset==='sharp')note='A layer of width 0.025 near x = 1 challenges the mesh. Refinement improves resolution; added isotropic diffusion changes the amplitude and does not restore the sharp layer.';
  else note=state.diffusion>0?'Artificial diffusion changes the Poisson equation by increasing its diffusivity. Compare with the exact reference: extra smoothness can increase error.':'A smooth starting point. Refine the mesh to see the finite element solution approach the exact solution.';
  $('insightText').textContent=note;
  for(const id of ['matrixRow','matrixCol','rhsRow'])$(id).max=r.N-1;
  $('nodeId').max=r.mesh.nodes.length-1;
  updateMetrics();
}
function updateMetrics(){if(!result)return;
  const values=fieldValues(result,field),m=field==='stable'?result.metrics.stable:result.metrics.raw;
  $('minValue').textContent=format(Math.min(...values));$('maxValue').textContent=format(Math.max(...values));
  $('l2Value').textContent=field==='difference'||field==='reference'?'—':m.l2===null?'n/a':m.l2.toExponential(2);
  const pe=field==='stable'?result.effectivePe:result.pe;$('peValue').textContent=Number.isFinite(pe)?short(pe):'∞';
  $('residualValue').textContent=['difference','reference'].includes(field)||m.backwardError===null?'—':m.backwardError.toExponential(1);
  if(trajectory){
    $('integralValue').textContent=field==='reference'?'—':(field==='difference'?result.metrics.stable.integral-result.metrics.raw.integral:m.integral).toPrecision(5);
    $('negativeValue').textContent=Array.from(values).filter(v=>v < -1e-10).length;
    $('courantValue').textContent=short(Math.hypot(...velocity(state))*result.dt*state.n);
    $('diffusionNumberValue').textContent=short(physicalDiffusion(state)*result.dt*state.n**2);
  }
}
function render(){if(!result)return;meshPoints=meshPlot($('meshPreview'),result.mesh,selected);
  const p=surfacePlot($('solution'),result,{...camera,view,field,mesh:showMesh,selected,slice});points=p.points;
  $('legendMin').textContent=short(p.lo);$('legendMax').textContent=short(p.hi);
  crossPlot($('crossSection'),result,solvedState,slice);matrixPlot($('matrixPreview'),{...result,A:inspectedMatrix()},+$('matrixRow').value||0,+$('matrixCol').value||0);
}
function inspectedMatrix(){if(!result)return null;const key=$('matrixPicker').value;return key==='A'?result.A:key==='step'?result.stableTimeMatrix:result[key];}
function setField(next){field=next;document.querySelectorAll('[data-field]').forEach(b=>{b.classList.toggle('active',b.dataset.field===field);b.setAttribute('aria-pressed',b.dataset.field===field);});if(result){updateSummary();render();}}
function selectEditor(next){editor=next;document.querySelectorAll('[data-editor]').forEach(b=>{b.classList.toggle('active',b.dataset.editor===editor);b.setAttribute('aria-pressed',b.dataset.editor===editor);});for(const key of ['matrix','rhs','node'])$(`${key}Editor`).classList.toggle('hidden',key!==editor);$('symmetricLabel').classList.toggle('hidden',editor!=='matrix');updateEditor();}
function validIndex(input,max){const v=Number(input.value);return input.value.trim()!==''&&Number.isInteger(v)&&v>=0&&v<max?v:null;}
function updateEditor(){if(!result)return;
  const i=validIndex($('matrixRow'),result.N),j=validIndex($('matrixCol'),result.N),k=validIndex($('rhsRow'),result.N),node=validIndex($('nodeId'),result.mesh.nodes.length);
  const readOnly=$('matrixPicker').value!=='A',matrix=inspectedMatrix(),key=$('matrixPicker').value;
  $('matrixValue').readOnly=readOnly;$('matrixValueLabel').textContent=key==='A'?'Aᵢⱼ':key==='ML'?'Mᴸᵢⱼ':key==='step'?'Bᵢⱼ':`${key}ᵢⱼ`;
  $('applyEdit').disabled=editor==='matrix'&&readOnly;$('symmetricLabel').classList.toggle('hidden',editor!=='matrix'||readOnly);
  if(i!==null&&j!==null)$('matrixValue').value=Number(matrix[i*result.N+j].toPrecision(10));
  if(k!==null)$('rhsValue').value=Number(result.F[k].toPrecision(10));
  if(node!==null){const[x,y]=result.mesh.nodes[node],dof=result.mesh.dof[node],boundary=dof<0;
    $('nodeValue').value=Number(result.raw[node].toPrecision(10));$('nodeValue').disabled=!boundary;
    $('nodeInfo').textContent=`(${short(x)}, ${short(y)}) · ${boundary?'Dirichlet boundary; value is editable':`interior DOF ${dof}; edit its row in A or f`}`;
    if(editor==='node')$('matrixDetail').textContent=`Node ${node}: u = ${format(result.raw[node])}, uδ = ${format(result.stable[node])}`;
  }
  if(editor==='matrix')$('matrixDetail').textContent=i===null||j===null?`Enter row and column in 0…${result.N-1}.`:readOnly?`${key==='step'?'Comparison B':key}[${i}, ${j}] = ${matrix[i*result.N+j].toPrecision(6)} · read only`:`A[${i}, ${j}] = ${short(result.A[i*result.N+j])} · δK = ${short(addedDiffusion(state)*result.K[i*result.N+j])} · SUPG S = ${short(result.S[i*result.N+j])} · total = ${short(result.effectiveA[i*result.N+j])}`;
  if(editor==='rhs')$('matrixDetail').textContent=k===null?`Enter a DOF in 0…${result.N-1}.`:`f[${k}] = ${short(result.F[k])} · stabilized load = ${short(result.effectiveF[k])} · b[${k}] = ${short(result.rawB[k])}`;
  matrixPlot($('matrixPreview'),{...result,A:matrix},i||0,j||0);
}
function clearIndexed(){state.matrixEdits={};state.rhsEdits={};state.nodeEdits={};selected=null;for(const id of ['matrixRow','matrixCol','rhsRow','nodeId'])$(id).value=0;}
function loadPreset(name){stopPlayback();state=preset(name);selectedTime=0;trajectory=null;$('matrixPicker').value='A';clearIndexed();$('symmetric').checked=!Math.hypot(...velocity(state));field='raw';syncControls();setField('raw');schedule();}

document.querySelectorAll('[data-preset]').forEach(b=>b.addEventListener('click',()=>loadPreset(b.dataset.preset)));
$('reset').addEventListener('click',()=>{loadPreset('poisson');camera={yaw:-.66,elevation:.67};toast('Reset to the smooth Poisson experiment.');});
$('meshSize').addEventListener('input',()=>{state.n=+$('meshSize').value;clearIndexed();syncControls();schedule();});
$('pattern').addEventListener('change',()=>{state.pattern=$('pattern').value;clearIndexed();schedule();});
$('mode').addEventListener('change',()=>loadPreset({poisson:'poisson',transport:'layer',advection:'moving-front',heat:'heat-pulse',parabolic:'moving-pulse',reaction:'reaction-decay','diffusion-reaction':'diffusion-reaction','advection-reaction':'advection-reaction',adr:'adr-pulse','steady-adr':'steady-adr'}[$('mode').value]));
for(const key of ['epsilon','reaction','betaX','betaY','source','amplitude'])$(key).addEventListener('input',()=>{
  const v=key==='source'?$(key).value:$(key).value.trim()===''?NaN:Number($(key).value);state[key]=v;if(key==='reaction')syncControls();schedule();
});
for(const key of ['initial','dt','duration','scheme','massType'])$(key).addEventListener(['scheme','massType'].includes(key)?'change':'input',()=>{
  state[key]=['dt','duration'].includes(key)?$(key).value.trim()===''?NaN:Number($(key).value):$(key).value;
  if(key==='massType')setField('stable');schedule();
});
$('matrixPicker').addEventListener('change',()=>{selectEditor('matrix');render();});
for(const side of ['left','right','bottom','top'])$(`bc-${side}`).addEventListener('input',()=>{state.boundaries[side]=$(`bc-${side}`).value;schedule();});
$('stabilization').addEventListener('change',()=>{state.stabilization=$('stabilization').value;syncControls();setField('stable');schedule();});
$('supgScale').addEventListener('input',()=>{state.supgScale=$('supgScale').value.trim()===''?NaN:Number($('supgScale').value);schedule();});
$('diffusion').addEventListener('input',()=>{state.diffusion=$('diffusion').value.trim()===''?NaN:Number($('diffusion').value);if(Number.isFinite(state.diffusion)){$('diffusionSlider').max=Math.max(.15,state.diffusion);$('diffusionSlider').value=state.diffusion;$('diffusionMax').textContent=short(+$('diffusionSlider').max);}schedule();});
$('diffusionSlider').addEventListener('input',()=>{state.diffusion=+$('diffusionSlider').value;$('diffusion').value=state.diffusion;schedule();});
$('suggestDiffusion').addEventListener('click',()=>{
  state.stabilization='isotropic';
  if(Math.hypot(...velocity(state))===0){toast('There is no convection term: the mesh-based value is δ = 0. You can still add diffusion manually.');state.diffusion=0;}
  else{state.diffusion=Math.max(0,Math.hypot(...velocity(state))/(2*state.n)-physicalDiffusion(state));toast(`Applied δ = ${state.diffusion.toPrecision(5)}. Compare the cross-section below.`);}
  syncControls();setField('stable');schedule();
});
document.querySelectorAll('[data-field]').forEach(b=>b.addEventListener('click',()=>setField(b.dataset.field)));
document.querySelectorAll('[data-view]').forEach(b=>b.addEventListener('click',()=>{view=b.dataset.view;document.querySelectorAll('[data-view]').forEach(x=>{x.classList.toggle('active',x.dataset.view===view);x.setAttribute('aria-pressed',x.dataset.view===view);});$('plotHint').innerHTML=view==='surface'?'<span>↔</span> Drag to orbit <span>·</span> Click a node to inspect':'Click a node to inspect its value';render();}));
$('meshToggle').addEventListener('click',()=>{showMesh=!showMesh;$('meshToggle').classList.toggle('active',showMesh);$('meshToggle').setAttribute('aria-pressed',showMesh);render();});
$('resetView').addEventListener('click',()=>{camera={yaw:-.66,elevation:.67};render();});
$('sliceY').addEventListener('change',()=>{const value=Number($('sliceY').value);if($('sliceY').value===''||!Number.isFinite(value)||value<0||value>1){$('sliceY').value=slice;toast('Choose a horizontal slice from y = 0 to y = 1.');return;}slice=value;render();});
document.querySelectorAll('[data-editor]').forEach(b=>b.addEventListener('click',()=>selectEditor(b.dataset.editor)));
for(const id of ['matrixRow','matrixCol','rhsRow','nodeId'])$(id).addEventListener('input',updateEditor);
$('clearEdits').addEventListener('click',()=>{clearIndexed();schedule();toast('Cleared all matrix, load, and boundary node overrides.');});
$('editForm').addEventListener('submit',e=>{
  e.preventDefault();if(!result||pending){toast('Wait for the current solve before applying an edit.');return;}
  const valueId=editor==='matrix'?'matrixValue':editor==='rhs'?'rhsValue':'nodeValue',value=Number($(valueId).value);
  if($(valueId).value.trim()===''||!Number.isFinite(value)){toast('Enter a finite number.');return;}
  if(editor==='matrix'){
    if($('matrixPicker').value!=='A'){toast('This assembled matrix is read only. Select spatial operator A to edit it.');return;}
    const i=validIndex($('matrixRow'),result.N),j=validIndex($('matrixCol'),result.N);if(i===null||j===null){toast(`Choose indices from 0 to ${result.N-1}.`);return;}
    state.matrixEdits[`${i},${j}`]=value;if($('symmetric').checked)state.matrixEdits[`${j},${i}`]=value;
  }else if(editor==='rhs'){
    const i=validIndex($('rhsRow'),result.N);if(i===null){toast(`Choose a DOF from 0 to ${result.N-1}.`);return;}state.rhsEdits[i]=value;
  }else{
    const id=validIndex($('nodeId'),result.mesh.nodes.length);if(id===null||result.mesh.dof[id]>=0){toast('Select a boundary node to prescribe its value. Interior nodes are solution unknowns.');return;}state.nodeEdits[id]=value;
  }
  schedule();
});

function nearest(coords,x,y,radius=20){let best=null,dist=radius;coords.forEach((p,i)=>{const d=Math.hypot(x-p[0],y-p[1]);if(d<dist){dist=d;best=i;}});return best;}
function inspectNode(id){if(id===null||!result||pending)return;selected=id;$('nodeId').value=id;const dof=result.mesh.dof[id];if(dof>=0){$('matrixRow').value=dof;$('matrixCol').value=dof;$('rhsRow').value=dof;}selectEditor('node');render();}
$('meshPreview').addEventListener('click',e=>{const rect=e.currentTarget.getBoundingClientRect();inspectNode(nearest(meshPoints,e.clientX-rect.left,e.clientY-rect.top,18));});
let drag=null;
$('solution').addEventListener('pointerdown',e=>{if(pending)return;drag={x:e.clientX,y:e.clientY,yaw:camera.yaw,elevation:camera.elevation,moved:false};e.currentTarget.setPointerCapture(e.pointerId);});
$('solution').addEventListener('pointermove',e=>{if(!result||pending)return;const rect=e.currentTarget.getBoundingClientRect(),x=e.clientX-rect.left,y=e.clientY-rect.top;
  if(drag&&view==='surface'){
    if(Math.hypot(e.clientX-drag.x,e.clientY-drag.y)>3)drag.moved=true;
    if(drag.moved){camera.yaw=drag.yaw+(e.clientX-drag.x)*.008;camera.elevation=Math.max(.15,Math.min(1.35,drag.elevation+(e.clientY-drag.y)*.006));$('nodeTooltip').classList.add('hidden');render();return;}
  }
  const id=nearest(points,x,y,15);if(id===null){$('nodeTooltip').classList.add('hidden');return;}
  const p=result.mesh.nodes[id],val=fieldValues(result,field)[id];$('nodeTooltip').textContent=`Node ${id} · (${short(p[0])}, ${short(p[1])}) · u = ${short(val)}`;$('nodeTooltip').style.left=`${Math.min(x+14,Math.max(5,rect.width-270))}px`;$('nodeTooltip').style.top=`${Math.max(5,y-35)}px`;$('nodeTooltip').classList.remove('hidden');
});
$('solution').addEventListener('pointerup',e=>{if(drag&&!drag.moved){const rect=e.currentTarget.getBoundingClientRect();inspectNode(nearest(points,e.clientX-rect.left,e.clientY-rect.top,24));}drag=null;});
$('solution').addEventListener('pointercancel',()=>drag=null);
$('solution').addEventListener('pointerleave',()=>$('nodeTooltip').classList.add('hidden'));
for(const id of ['guideButton','insightGuide'])$(id).addEventListener('click',()=>$('guide').showModal());
$('closeGuide').addEventListener('click',()=>$('guide').close());
$('guide').addEventListener('click',e=>{if(e.target===$('guide')){const r=$('guide').getBoundingClientRect();if(e.clientX<r.left||e.clientX>r.right||e.clientY<r.top||e.clientY>r.bottom)$('guide').close();}});
$('export').addEventListener('click',()=>{
  if(!result||pending)return;const r=result;
  const sparse=matrix=>{const entries=[];for(let i=0;i<r.N;i++)for(let j=0;j<r.N;j++)if(matrix[i*r.N+j]!==0)entries.push([i,j,matrix[i*r.N+j]]);return entries;};
  const data={format:'elliptic-lab-v3',description:'Zero-based indices. A and f are the reduced edited base system. Time-dependent runs compare consistent-mass baseline with the selected mass and added diffusion. SUPG uses the physical PDE residual; S and H are its spatial and time corrections. Time-step loads include changing-boundary mass lifting including H. Pe excludes directional SUPG.',state:solvedState,
    nodes:r.mesh.nodes,triangles:r.mesh.triangles,interiorNodeIds:r.mesh.interior,solutions:{original:Array.from(r.raw),stabilized:Array.from(r.stable),exactAtNodes:r.reference?Array.from(r.reference):null},
    system:{size:r.N,A:sparse(r.A),K:sparse(r.K),M:sparse(r.M),supgSpatial:sparse(r.S),supgTime:sparse(r.H),tauRange:r.tauRange,lumpedMassDiagonal:Array.from(r.lumped),f:Array.from(r.F),stabilizedA:sparse(r.effectiveA),stabilizedF:Array.from(r.effectiveF)},metrics:r.metrics};
  if(trajectory){data.system.comparisonMass=sparse(r.comparisonMass);data.time=r.time;data.trajectory=trajectory.frames.map(f=>({time:f.time,dt:f.dt,original:Array.from(f.raw),comparison:Array.from(f.stable),metrics:f.metrics}));data.system.timeStepMatrix=sparse(r.timeMatrix);data.system.comparisonTimeStepMatrix=sparse(r.stableTimeMatrix);data.system.timeStepRhs=r.rawRhs?Array.from(r.rawRhs):null;data.system.comparisonTimeStepRhs=r.stableRhs?Array.from(r.stableRhs):null;}
  const url=URL.createObjectURL(new Blob([JSON.stringify(data,null,2)],{type:'application/json'})),a=document.createElement('a');a.href=url;a.download=`elliptic-lab-${state.preset}-n${state.n}.json`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);toast('Exported experiment, mesh, matrices, and solutions as JSON.');
});
new ResizeObserver(()=>requestAnimationFrame(render)).observe(document.querySelector('main'));
syncControls();schedule();
