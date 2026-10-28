import {get as matrixGet,entries as matrixEntries} from './matrix.js';
import {isAfc} from './afc.js';
import {experiments,families,equationFamily,timeSetting,examplesFor,defaultExample,compatibleMethods} from './experiments.js';
import {preset,exactSolution,isTransient,velocity,physicalDiffusion,addedDiffusion} from './fem.js';
import {meshPlot,surfacePlot,crossPlot,matrixPlot,fieldValues,short} from './plots.js';
const $=id=>document.getElementById(id);
const requestedExample=new URLSearchParams(location.search).get('example');
const knownExamples=experiments.map(e=>e.id);
let state=preset(knownExamples.includes(requestedExample)?requestedExample:'afc-layer'),result=null,solvedState=null,requestId=0,timer,toastTimer,selected=null,editor='matrix';
let camera={yaw:-.66,elevation:.67},view='surface',field='raw',showMesh=true,slice=.5,points=[],meshPoints=[],pending=false;
let trajectory=null,frameIndex=0,selectedTime=0,playing=false,playRequest=0,lastPlayback=0;
let worker;
function startWorker(){worker=new Worker(new URL('./worker.js',import.meta.url),{type:'module'});worker.onmessage=handleWorkerMessage;worker.onerror=()=>failure('The solver could not start. Open the app through its web server, then reload.');}
const methodLabel=()=>isAfc(state)?(state.stabilization==='afc'?(isTransient(state)?'Zalesak FCT':`BJK AFC · μ = ${state.afcMu}`):'low-order graph diffusion'):state.stabilization==='supg'?`SUPG × ${state.supgScale}`:state.stabilization==='none'?'no stabilization':`isotropic δ = ${short(state.diffusion)}`;
const format=v=>Math.abs(v)<1e-12?'0.000':Math.abs(v)<.001||Math.abs(v)>9999?v.toExponential(2):v.toFixed(3);
function toast(message){$('toast').textContent=message;$('toast').classList.remove('hidden');clearTimeout(toastTimer);toastTimer=setTimeout(()=>$('toast').classList.add('hidden'),4000);}
function failure(message){$('errorBanner').textContent=message+' The plots show the last successful solve.';$('errorBanner').classList.remove('hidden');$('solveState').textContent='Check inputs';document.querySelector('.plot-card').classList.add('stale');$('export').disabled=true;pending=false;}
function schedule(){stopPlayback();clearTimeout(timer);if(pending&&worker){worker.terminate();startWorker();}pending=true;$('errorBanner').classList.add('hidden');$('solveState').textContent=isTransient(state)?'Building trajectory…':'Solving…';$('export').disabled=true;for(const id of ['playTime','previousTime','nextTime','rewindTime','timeSlider'])$(id).disabled=true;document.querySelector('.plot-card').classList.add('stale');requestId++;const id=requestId;timer=setTimeout(()=>worker.postMessage({id,state:structuredClone(state)}),90);}
function handleWorkerMessage({data}){
  if(data.id!==requestId)return;
  if(data.progress!==undefined){$('solveState').textContent=`Time steps · ${Math.round(data.progress*100)}%`;return;}
  if(data.error){failure(data.error);return;}
  pending=false;solvedState=structuredClone(state);$('solveState').textContent='✓ Solved';$('export').disabled=false;document.querySelector('.plot-card').classList.remove('stale');
  if(data.result.afc?.converged===false)$('solveState').textContent='AFC iteration limit';
  if(data.result.frames){trajectory=data.result;const index=trajectory.frames.reduce((best,f,i)=>Math.abs(f.time-selectedTime)<Math.abs(trajectory.frames[best].time-selectedTime)?i:best,0);selectFrame(index);}
  else{trajectory=null;result=data.result;}
  if(field==='reference'&&!result.reference)setField('raw');
  updateSummary();updateEditor();render();
};
startWorker();

function stopPlayback(){playing=false;cancelAnimationFrame(playRequest);if($('playTime')){$('playTime').textContent='▶ Play';$('playTime').setAttribute('aria-label','Play time evolution');}}
function selectFrame(index){
  if(!trajectory)return;frameIndex=Math.max(0,Math.min(trajectory.frames.length-1,index));
  const f=trajectory.frames[frameIndex];selectedTime=f.time;
  result={...trajectory,...f,rawB:Float64Array.from(trajectory.mesh.interior,id=>f.raw[id]),stableB:Float64Array.from(trajectory.mesh.interior,id=>f.stable[id]),
    timeMatrix:frameIndex===trajectory.steps?trajectory.lastTimeMatrix:trajectory.timeMatrix,
    stableTimeMatrix:frameIndex===trajectory.steps?trajectory.lastStableTimeMatrix:trajectory.stableTimeMatrix};
  $('timeSlider').max=trajectory.steps;$('timeSlider').value=frameIndex;$('currentTime').textContent=`t = ${format(f.time)}`;
  $('frameCounter').textContent=`Frame ${frameIndex} / ${trajectory.steps}${trajectory.saveStride>1?` · ${trajectory.computedSteps} computed steps`:''}`;
  $('timelineNote').textContent=`${state.scheme==='be'?'Backward Euler':'Crank–Nicolson'} · Δt = ${Number(f.dt.toPrecision(5))} · comparison: ${isAfc(state)?'SSPRK2 · lumped':state.massType==='lumped'?'lumped':'consistent'} M · ${methodLabel()}`;
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

function syncProblemChoices(){
  const time=timeSetting(state),family=equationFamily(state),current=experiments.find(e=>e.id===state.preset);
  $('problemTime').value=time;$('equationFamily').value=family;
  const matches=examplesFor(time,family),options=matches.map(e=>new Option(e.label,e.id));
  if(!matches.some(e=>e.id===state.preset))options.unshift(new Option(`Custom settings · ${current?.label??state.preset}`,state.preset));
  $('quickExperiment').replaceChildren(...options);$('quickExperiment').value=state.preset;
  $('exampleHelp').textContent=current?.description??'Adjust the coefficients and data for this equation.';
  $('experimentContext').textContent=`${time==='steady'?'STEADY':'TRANSIENT'} / ${family} / ${families[family]}`;
  $('experimentTitle').textContent=current?.label??'Custom experiment';$('experimentDescription').textContent=current?.description??'';
  const note=time==='steady'&&family==='C'?'Pure steady convection usually takes inflow data. This example uses compatible traces on every edge. Arbitrary boundary data or some meshes can give a singular Galerkin system; use CD with small ε to study unresolved layers.':time==='steady'&&family==='R'?'For steady pure reaction, c must be nonzero to determine u. Boundary traces in this example match the equilibrium.':'';
  $('problemNote').textContent=note;$('problemNote').classList.toggle('hidden',!note);
  $('stabilization').replaceChildren(...compatibleMethods(state).map(m=>{const option=new Option(m.label,m.id);option.disabled=!m.enabled;return option;}));
}
function syncControls(){
  syncProblemChoices();
  for(const key of ['epsilon','reaction','betaX','betaY','source','amplitude','pattern','diffusion','stabilization','supgScale','afcMu','afcMaxIterations','afcTolerance','afcRelaxation'])$(key).value=state[key];
  for(const side of ['left','right','bottom','top'])$(`bc-${side}`).value=state.boundaries[side];
  $('meshSize').value=state.n;$('meshSizeLabel').textContent=`${state.n} × ${state.n}`;
  $('diffusionSlider').max=Math.max(.15,state.diffusion);$('diffusionSlider').value=state.diffusion;$('diffusionMax').textContent=short(+$('diffusionSlider').max);
  $('afcControls').classList.toggle('hidden',!isAfc(state));$('afcSteadyControls').classList.toggle('hidden',isTransient(state));
  $('afcMu').disabled=state.stabilization==='afc-low';
  document.querySelector('#afcControls > .small-note').textContent=isTransient(state)?'Zalesak-style FCT with lumped mass and SSPRK2 CFL substeps. c ≥ 0; clear A/f overrides.':'BJK-type steady AFC: nonlinear limiter iteration with a reused linear solver (LU through n = 32; sparse ILU-preconditioned Krylov above). c ≥ 0; clear A/f overrides.';
  $('massType').disabled=isAfc(state);
  $('massHelp').textContent=isAfc(state)?'AFC comparison: lumped mass, SSPRK2, automatic CFL substeps. Baseline: consistent mass and the selected time integrator. Reaction stays in the assembled operator.':'Original: consistent M, δ = 0. Comparison: selected mass and stabilization. Reaction cM stays consistent in both. Both use the same time integrator. Up to 240 computed steps; fine meshes save fewer playback frames.';
  $('supgControls').classList.toggle('hidden',state.stabilization!=='supg');$('isotropicControls').classList.toggle('hidden',state.stabilization!=='isotropic');
  const temporal=isTransient(state),family=equationFamily(state),convective=family.includes('C');
  $('boundarySectionNumber').textContent=temporal?'05':'04';$('stabilizationSectionNumber').textContent=temporal?'06':'05';
  $('fieldBadge').textContent=`${family} · ${temporal?'TRANSIENT':'STEADY'}`;
  $('rotationHelp').textContent='Rotation: β(x,y) = (0.5 − y, x − 0.5). Counterclockwise; period 2π. Initial profile: rotation(x,y). Zero traces on all edges.'+(physicalDiffusion(state)>0?' Physical diffusion spreads the shapes; no exact reference is supplied for this bounded-domain problem.':' The exact rotating support stays inside the square.');
  $('rotationHelp').classList.toggle('hidden',state.flow!=='rotation');$('velocityFields').classList.toggle('hidden',!convective||state.flow==='rotation');$('epsilon').disabled=!family.includes('D');$('reaction').disabled=!family.includes('R');
  for(const id of ['timeControls','timeline','timeDiagnostics'])$(id).classList.toggle('hidden',!temporal);
  $('matrixPickerWrap').classList.remove('hidden');
  for(const option of $('matrixPicker').options)option.disabled=(!temporal&&['ML','step','comparisonMass'].includes(option.value))||(isAfc(state)&&temporal&&option.value==='step')||(!isAfc(state)&&option.value==='afcD');
  if($('matrixPicker').selectedOptions[0]?.disabled)$('matrixPicker').value='A';
  for(const key of ['initial','dt','duration','scheme','massType'])$(key).value=state[key];
  $('pdeLabel').textContent=[temporal?'uₜ':'',family.includes('D')?'− εΔu':'',convective?'+ β·∇u':'',family.includes('R')?'+ cu':''].filter(Boolean).join(' ').replace(/^\+ /,'')+' = f';
  $('discreteLabel').textContent=temporal?'Mḃ + Ab = f':'Ab = f';$('inspectorEquation').textContent=temporal?'Mḃ + Ab = f':'Ab = f';
  $('originalFieldName').textContent=$('originalLegendName').textContent=temporal?'Consistent baseline':'Original';
  $('comparisonFieldName').textContent=$('comparisonLegendName').textContent=temporal?'Comparison':isAfc(state)?'With AFC':state.stabilization==='supg'?'With SUPG':state.stabilization==='none'?'Comparison':'With diffusion';
}
function updateSummary(){
  const r=result,afc=r.afc;
  $('afcStatus').textContent=!afc?(trajectory&&frameIndex===0?'Initial state · 0 AFC substeps completed.':'Solve to see AFC cost.'):afc.kind==='steady'?`${afc.converged?'Converged':'ITERATION LIMIT — last iterate'} · ${afc.iterations} iteration${afc.iterations===1?'':'s'} · residual ${afc.nonlinearResidual.toExponential(2)} · ${(100*afc.limitedFraction).toFixed(0)}% limited edges · ${afc.ms.toFixed(0)} ms AFC`:`SSPRK2 · ${afc.substeps} substep${afc.substeps===1?'':'s'} this frame · ${afc.totalSubsteps} so far · ${afc.limiterPasses} limiter passes this frame · ${(100*afc.limitedFraction).toFixed(0)}% limited edges · ${afc.ms.toFixed(0)} ms AFC${afc.dataNonnegative?'':' · signed data: positivity not asserted'}`;
  const affine=afc?.affine,probe=afc?.affineLimiter;
  $('affineStatus').classList.toggle('hidden',!affine||!isAfc(state));
  $('affineStatus').textContent=affine?`Affine test u = 1+x+y · max nodal error ${affine.maxNodalError.toExponential(2)} (Galerkin ${affine.baselineMaxNodalError.toExponential(2)}) · max retained affine edge flux ${probe.maxLimitedFlux.toExponential(2)} · retained-diffusion residual at exact field ${probe.retainedResidual.toExponential(2)} · μ ≥ ${probe.requiredMu.toFixed(4)} removes limiting on this field only.${afc.converged?'':' Solve unconverged: nodal error includes iteration error.'}`:'';
  $('tauSummary').textContent=`Element τ: ${r.tauRange.map(v=>short(v)).join(' – ')}${Math.hypot(...velocity(state))===0?' · no flow, SUPG is zero':''}`;$('nodeCount').textContent=r.mesh.nodes.length;$('elementCount').textContent=r.mesh.triangles.length;$('dofCount').textContent=r.N;
  $('systemSize').textContent=`${r.N} unknowns`;$('matrixShape').textContent=`${r.N} × ${r.N}`;
  const count=Object.keys(state.matrixEdits).length+Object.keys(state.rhsEdits).length+Object.keys(state.nodeEdits).length;$('editCount').textContent=`${count} override${count===1?'':'s'}`;
  $('timing').textContent=`${r.mesh.triangles.length} elements · ${r.ms.toFixed(0)} ms · ${isAfc(state)?(trajectory?`${trajectory.frames.length} saved frames · ${trajectory.afc?.totalSubsteps??'pending'} AFC substeps`:`${afc?.iterations??0} AFC iterations${afc?.converged===false?' · NOT CONVERGED':''}`):trajectory?`${trajectory.computedSteps} time steps · ${state.n>32?'sparse Krylov':'reused LU'}`:state.n>32?'sparse Krylov solve':'partial-pivot solve'}`;
  document.querySelector('[data-field="reference"]').disabled=!r.reference;
  $('exactLegend').classList.toggle('hidden',!r.reference);
  $('fieldBadge').textContent=`${equationFamily(state)} · ${isTransient(state)?'TRANSIENT':'STEADY'}`;
  $('fieldSubtitle').textContent=field==='reference'?'Exact solution sampled at mesh nodes':field==='difference'?'Comparison minus original · shared coordinates':trajectory?`t = ${format(r.time)} · ${field==='stable'&&state.massType==='lumped'?'lumped mass':'consistent mass'} · ${field==='stable'&&isAfc(state)?'SSPRK2 + CFL substeps':state.scheme==='be'?'backward Euler':'Crank–Nicolson'}${field==='stable'?` · ${methodLabel()}`:' · δ = 0'}`:'Piecewise-linear finite element solution';
  let note;
  if(count)note=`${count} algebraic or boundary override${count===1?' is':'s are'} active. The original PDE reference no longer applies. Changes to mesh or triangulation clear these edits.`;
  else if(trajectory&&r.initialMismatch)note=`Initial and boundary data disagree at ${r.initialMismatch} boundary nodes. Prescribed boundary values take precedence; the two runs share the same initial interior field.`;
  else if(state.flow==='rotation'&&physicalDiffusion(state)>0)note='Rotating and diffusing: ε is physical diffusion and acts in both solutions. Compare ε = 0, 0.0001, and 0.001 over T = 2π, then compare FCT with low-order graph diffusion. Positive ε spreads the shapes, so they do not return to their initial profiles after a revolution. No exact reference or L² error is supplied for this bounded-domain case. Zero boundary traces can remove mass as the field spreads.';
  else if(state.flow==='rotation')note='Solid-body rotation: one revolution takes 2π. Scrub the timeline and compare the three shapes with the exact reference. The default n = 32 triangular mesh resolves the slot only coarsely; try 64 or 128. FCT uses the existing lumped SSPRK2 method, not the paper’s CN-FCT algorithms. Zero traces are imposed on all edges; monitor integral drift as well as extrema.';
  else if(affine)note=`Convection-dominated affine test: −0.005Δu + (1,0.3)·∇u = 1.3, exact u = 1+x+y (default coefficients). On the original centrally symmetric meshes, μ = 1 is sufficient for linearity preservation. Try μ = 0.25 and 1, then the optional shifted-node mesh with μ = 0.5 and 1. Check retained edge fluxes as well as nodal error: cancellation can reproduce an affine field despite active limiting. The reported threshold concerns this field only.${afc?.converged===false?' The nonlinear solve is not converged.':''}`;
  else if(isAfc(state)&&state.stabilization==='afc-low')note='The low-order comparison retains all graph diffusion. It needs no nonlinear limiter iteration and is deliberately diffusive. Switch to AFC limited flux correction to recover sharper layers.';
  else if(isAfc(state))note=trajectory?'AFC limits antidiffusive edge fluxes without clipping nodal values. Compare low-order graph diffusion to see what flux correction restores. The comparison uses SSPRK2 substeps; the baseline uses the selected time integrator. Positivity and accuracy are separate checks.':`AFC uses a nonlinear BJK-type limiter and reuses one low-order LU factorization. ${afc?.converged===false?'The iteration limit was reached; this is the last unconverged iterate.':'Compare the sharpened layer with the low-order graph-diffusion solution.'}`;
  else if(state.stabilization==='supg')note=`SUPG adds residual stabilization along the flow (τ = ${r.tauRange.map(v=>short(v)).join(' – ')}). Compare extrema and L² error: oscillations can remain. Switch to isotropic diffusion to compare spreading across the flow; refine the mesh to resolve the layer.`;
  else if(trajectory){
    if(r.initialMismatch)note=`Initial and boundary data disagree at ${r.initialMismatch} boundary nodes. Prescribed boundary data take precedence. Playback compares the same initial interior field.`;
    else if(state.reaction!==0)note=`Reaction c = ${state.reaction}: ${state.reaction>0?'decay':'growth'}. Compare the exact profile at the same time. Added diffusion spreads the field; reaction changes its amplitude. Mass lumping changes only the time derivative, while cM remains consistent.`;
    else if(state.mode==='heat')note='Step forward from the narrow initial pulse. Lumped mass with backward Euler preserves nonnegativity here; compare negative nodes and L² error. Positivity and accuracy are different questions.';
    else if(state.mode==='advection')note='Watch the front move with the flow. Lumping alone can still produce oscillations; compare with added diffusion and switch between Crank–Nicolson and backward Euler to see temporal damping.';
    else note='Physical diffusion spreads the pulse as convection carries it. Compare the effect of extra δK and mass lumping against the exact profile at the same time.';
  }
  else if(state.preset==='steady-c')note='Pure steady convection with compatible affine traces. The Galerkin solution reproduces the exact plane on the starter mesh. Some other meshes give a singular baseline; use CD with positive ε for boundary-layer experiments.';
  else if(['steady-r','steady-cr'].includes(state.preset))note='The exact affine equilibrium balances convection/reaction with its manufactured source. Physical terms absent from the selected family are zero. Edit the source or traces to explore a different balance.';
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
function inspectedMatrix(){if(!result)return null;const key=$('matrixPicker').value,matrix=key==='A'?result.A:key==='step'?result.stableTimeMatrix:result[key];if(matrix)return matrix;$('matrixPicker').value='A';return result.A;}
function setField(next){field=next;document.querySelectorAll('[data-field]').forEach(b=>{b.classList.toggle('active',b.dataset.field===field);b.setAttribute('aria-pressed',b.dataset.field===field);});if(result){updateSummary();render();}}
function selectEditor(next){editor=next;document.querySelectorAll('[data-editor]').forEach(b=>{b.classList.toggle('active',b.dataset.editor===editor);b.setAttribute('aria-pressed',b.dataset.editor===editor);});for(const key of ['matrix','rhs','node'])$(`${key}Editor`).classList.toggle('hidden',key!==editor);$('symmetricLabel').classList.toggle('hidden',editor!=='matrix');updateEditor();}
function validIndex(input,max){const v=Number(input.value);return input.value.trim()!==''&&Number.isInteger(v)&&v>=0&&v<max?v:null;}
function updateEditor(){if(!result)return;
  const i=validIndex($('matrixRow'),result.N),j=validIndex($('matrixCol'),result.N),k=validIndex($('rhsRow'),result.N),node=validIndex($('nodeId'),result.mesh.nodes.length);
  const matrix=inspectedMatrix(),key=$('matrixPicker').value,readOnly=key!=='A';
  $('matrixValue').readOnly=readOnly;$('matrixValueLabel').textContent=key==='A'?'Aᵢⱼ':key==='ML'?'Mᴸᵢⱼ':key==='step'?'Bᵢⱼ':`${key}ᵢⱼ`;
  $('applyEdit').disabled=editor==='matrix'&&readOnly;$('symmetricLabel').classList.toggle('hidden',editor!=='matrix'||readOnly);
  if(i!==null&&j!==null)$('matrixValue').value=Number(matrixGet(matrix,i,j).toPrecision(10));
  if(k!==null)$('rhsValue').value=Number(result.F[k].toPrecision(10));
  if(node!==null){const[x,y]=result.mesh.nodes[node],dof=result.mesh.dof[node],boundary=dof<0;
    $('nodeValue').value=Number(result.raw[node].toPrecision(10));$('nodeValue').disabled=!boundary;
    $('nodeInfo').textContent=`(${short(x)}, ${short(y)}) · ${boundary?'Dirichlet boundary; value is editable':`interior DOF ${dof}; edit its row in A or f`}`;
    if(editor==='node')$('matrixDetail').textContent=`Node ${node}: u = ${format(result.raw[node])}, uδ = ${format(result.stable[node])}`;
  }
  if(editor==='matrix')$('matrixDetail').textContent=i===null||j===null?`Enter row and column in 0…${result.N-1}.`:readOnly?`${key==='step'?'Comparison B':key}[${i}, ${j}] = ${matrixGet(matrix,i,j).toPrecision(6)} · read only`:`A[${i}, ${j}] = ${short(matrixGet(result.A,i,j))} · δK = ${short(addedDiffusion(state)*matrixGet(result.K,i,j))} · SUPG S = ${short(matrixGet(result.S,i,j))}${result.afcD?` · AFC ${trajectory?'low-order D':'retained D'} = ${short(matrixGet(result.effectiveA,i,j)-matrixGet(result.A,i,j))}`:''} · total = ${short(matrixGet(result.effectiveA,i,j))}`;
  if(editor==='rhs')$('matrixDetail').textContent=k===null?`Enter a DOF in 0…${result.N-1}.`:`f[${k}] = ${short(result.F[k])} · stabilized load = ${short(result.effectiveF[k])} · b[${k}] = ${short(result.rawB[k])}`;
  matrixPlot($('matrixPreview'),{...result,A:matrix},i||0,j||0);
}
function clearIndexed(){state.matrixEdits={};state.rhsEdits={};state.nodeEdits={};selected=null;for(const id of ['matrixRow','matrixCol','rhsRow','nodeId'])$(id).value=0;}
function loadPreset(name){const url=new URL(location.href);url.searchParams.set('example',name);history.replaceState(null,'',url);stopPlayback();state=preset(name);selectedTime=0;trajectory=null;$('matrixPicker').value='A';clearIndexed();$('symmetric').checked=!Math.hypot(...velocity(state));field='raw';syncControls();setField('raw');schedule();}

$('problemTime').addEventListener('change',()=>loadPreset(defaultExample($('problemTime').value,$('equationFamily').value)));
$('equationFamily').addEventListener('change',()=>loadPreset(defaultExample($('problemTime').value,$('equationFamily').value)));
$('quickExperiment').addEventListener('change',()=>{if($('quickExperiment').value)loadPreset($('quickExperiment').value);});
$('reset').addEventListener('click',()=>{loadPreset('poisson');camera={yaw:-.66,elevation:.67};toast('Reset to the smooth Poisson experiment.');});
$('meshSize').addEventListener('input',()=>{state.n=+$('meshSize').value;clearIndexed();syncControls();schedule();});
$('pattern').addEventListener('change',()=>{state.pattern=$('pattern').value;clearIndexed();schedule();});
for(const key of ['epsilon','reaction','betaX','betaY','source','amplitude'])$(key).addEventListener('input',()=>{
  const v=key==='source'?$(key).value:$(key).value.trim()===''?NaN:Number($(key).value);state[key]=v;schedule();
});
for(const key of ['initial','dt','duration','scheme','massType'])$(key).addEventListener(['scheme','massType'].includes(key)?'change':'input',()=>{
  state[key]=['dt','duration'].includes(key)?$(key).value.trim()===''?NaN:Number($(key).value):$(key).value;
  if(key==='massType')setField('stable');schedule();
});
$('matrixPicker').addEventListener('change',()=>{selectEditor('matrix');render();});
for(const side of ['left','right','bottom','top'])$(`bc-${side}`).addEventListener('input',()=>{state.boundaries[side]=$(`bc-${side}`).value;schedule();});
$('stabilization').addEventListener('change',()=>{state.stabilization=$('stabilization').value;if(isAfc(state))state.massType='lumped';syncControls();setField('stable');schedule();});
for(const key of ['afcMu','afcMaxIterations','afcTolerance','afcRelaxation'])$(key).addEventListener('input',()=>{state[key]=$(key).value.trim()===''?NaN:Number($(key).value);schedule();});
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
  const sparse=matrixEntries;
  const data={format:'elliptic-lab-v4',description:'Zero-based indices. A and f are the reduced edited base system. Time-dependent runs compare consistent-mass baseline with the selected mass and added diffusion. SUPG uses the physical PDE residual; S and H are its spatial and time corrections. Time-step loads include changing-boundary mass lifting including H. Pe excludes directional SUPG and graph viscosity. AFC transient comparisons use lumped SSPRK2 FCT with CFL substeps, not an implicit step matrix; no linear backward error is reported for them. For transient AFC, stabilizedA and stabilizedF describe the low-order backbone; limited fluxes are applied stage by stage. For steady AFC, stabilizedA includes the retained graph diffusion at the returned iterate, and metadata report convergence.',state:solvedState,experiment:{timeSetting:timeSetting(solvedState),equationTerms:equationFamily(solvedState),example:experiments.find(e=>e.id===solvedState.preset)?.label,method:compatibleMethods(solvedState).find(m=>m.id===solvedState.stabilization)?.label},
    nodes:r.mesh.nodes,triangles:r.mesh.triangles,interiorNodeIds:r.mesh.interior,solutions:{original:Array.from(r.raw),stabilized:Array.from(r.stable),exactAtNodes:r.reference?Array.from(r.reference):null},
    system:{size:r.N,A:sparse(r.A),K:sparse(r.K),M:sparse(r.M),supgSpatial:sparse(r.S),supgTime:sparse(r.H),tauRange:r.tauRange,afc:r.afc,afcD:r.afcD?sparse(r.afcD):null,lumpedMassDiagonal:Array.from(r.lumped),f:Array.from(r.F),stabilizedA:sparse(r.effectiveA),stabilizedF:Array.from(r.effectiveF)},metrics:r.metrics};
  if(trajectory){data.system.comparisonMass=sparse(r.comparisonMass);data.time=r.time;data.computedSteps=trajectory.computedSteps;data.savedFrameStride=trajectory.saveStride;data.trajectory=trajectory.frames.map(f=>({time:f.time,dt:f.dt,original:Array.from(f.raw),comparison:Array.from(f.stable),metrics:f.metrics,afc:f.afc}));data.system.timeStepMatrix=sparse(r.timeMatrix);data.system.comparisonTimeStepMatrix=r.stableTimeMatrix?sparse(r.stableTimeMatrix):null;data.system.timeStepRhs=r.rawRhs?Array.from(r.rawRhs):null;data.system.comparisonTimeStepRhs=r.stableRhs?Array.from(r.stableRhs):null;}
  const url=URL.createObjectURL(new Blob([JSON.stringify(data,null,2)],{type:'application/json'})),a=document.createElement('a');a.href=url;a.download=`cdr-fem-fct-lab-${state.preset}-n${state.n}.json`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);toast('Exported experiment, mesh, matrices, and solutions as JSON.');
});
new ResizeObserver(()=>requestAnimationFrame(render)).observe(document.querySelector('main'));
syncControls();schedule();
