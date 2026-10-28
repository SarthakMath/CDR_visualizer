import {matrix,set,combine,rows as storedRows} from './matrix.js';
import {isAfc,afcTimeStepper,afcReducedLoad} from './afc.js';
import {assemble,expression,boundaryEvaluator,sourceLoad,diagnostics,exactSolution,physicalDiffusion,velocity,velocityAt,isTransient,addedDiffusion,supgLoad} from './fem.js';

export const matrixRows=storedRows;
export function multiplyRows(rows,b){return Float64Array.from(rows,row=>row.reduce((sum,[j,v])=>sum+v*b[j],0));}

import {factorDense,solveFactored} from './linear.js';
export {factorDense,solveFactored} from './linear.js';

function reducedForcing(system,state,load,g,supgF){
  const F=load.slice(),lift=new Float64Array(system.N);
  for(const[row,id,a,k]of system.boundaryTerms){F[row]-=a*g[id];lift[row]-=k*g[id];}
  for(const[i,v]of Object.entries(state.rhsEdits))F[+i]=v;
  const effectiveF=F.map((v,i)=>v+addedDiffusion(state)*lift[i]+supgF[i]);
  for(const[row,id,,,,s]of system.boundaryTerms)effectiveF[row]-=s*g[id];
  return{F,effectiveF};
}
export function computeTrajectory(state,progress=()=>{}){
  const started=performance.now();
  if(!isTransient(state))throw Error('Choose a time-dependent equation.');
  if(!Number.isFinite(state.dt)||state.dt<=0||!Number.isFinite(state.duration)||state.duration<=0)throw Error('Time step and end time must be positive.');
  if(!['be','cn'].includes(state.scheme)||!['consistent','lumped'].includes(state.massType))throw Error('Choose a supported time integrator and mass matrix.');
  const steps=Math.max(1,Math.ceil(state.duration/state.dt-1e-12));
  if(steps>240)throw Error('This interactive run supports at most 240 steps. Increase Δt or reduce the end time.');
  const system=assemble(state),{N,A,K,M,lumped,mesh}=system,theta=state.scheme==='be'?1:.5;
  const ML=matrix(N);for(let i=0;i<N;i++)set(ML,i,i,lumped[i]);
  const afcStepper=isAfc(state)?afcTimeStepper(state,system):null;let afcStats=null;
  const comparisonMass=combine([afcStepper||state.massType==='lumped'?ML:M,1],[system.H,1]),effectiveA=afcStepper?afcStepper.graph.lowOperator:combine([A,1],[K,addedDiffusion(state)],[system.S,1]);
  const identical=addedDiffusion(state)===0&&state.massType==='consistent'&&state.stabilization!=='supg'&&!afcStepper,factorCache=new Map();
  function stepSystems(dt){
    const key=dt.toPrecision(12);if(factorCache.has(key))return factorCache.get(key);
    const make=(mass,operator)=>{
      const B=combine([mass,1],[operator,theta*dt]),R=combine([mass,1],[operator,-(1-theta)*dt]);
      return{B,rows:matrixRows(B,N),rightRows:matrixRows(R,N),factor:factorDense(B,N)};
    };
    const raw=make(M,A),stable=identical||afcStepper?raw:make(comparisonMass,effectiveA),pair={raw,stable};factorCache.set(key,pair);return pair;
  }
  let g=boundaryEvaluator(state,mesh)(0);const evaluateBoundary=boundaryEvaluator(state,mesh),source=expression(state.source),initial=expression(state.initial),eps=physicalDiffusion(state),beta=velocity(state);
  let rawB=new Float64Array(N),stableB;
  for(let i=0;i<N;i++){const[x,y]=mesh.nodes[mesh.interior[i]];rawB[i]=initial(x,y,eps,0,...velocityAt(state,x,y),state.reaction??0);}stableB=rawB.slice();
  let initialMismatch=0;for(const id of mesh.boundary){const[x,y]=mesh.nodes[id];if(Math.abs(initial(x,y,eps,0,...velocityAt(state,x,y),state.reaction??0)-g[id])>1e-8)initialMismatch++;}
  const sourceVaries=/\bt\b/.test(state.source),exact=exactSolution(state);
  let forces=reducedForcing(system,state,system.load,g,system.supgF),lastPair=null;
  if(afcStepper)forces.effectiveF=afcReducedLoad(afcStepper.graph,system.load,g);
  const saveStride=state.n>64?Math.ceil(steps/80):state.n>32?Math.ceil(steps/120):1;
  const frames=[];let lo=0,hi=0,differenceMax=0;
  function snapshot(time,dt,rawRhs=null,stableRhs=null,pair=null){
    const raw=g.slice(),stable=g.slice();mesh.interior.forEach((id,i)=>{raw[id]=rawB[i];stable[id]=stableB[i];});
    const reference=exact?Float64Array.from(mesh.nodes,([x,y])=>exact(x,y,time)):null,exactNow=exact?(x,y)=>exact(x,y,time):null;
    const rawMetric=diagnostics(pair?.raw.B||M,rawRhs||multiplyRows(matrixRows(M,N),rawB),rawB,raw,mesh,exactNow,pair?.raw.rows);
    const stableMetric=identical?{...rawMetric}:diagnostics(pair?.stable.B||comparisonMass,stableRhs||multiplyRows(matrixRows(comparisonMass,N),stableB),stableB,stable,mesh,exactNow,pair?.stable.rows);
    if(afcStepper){stableMetric.residual=null;stableMetric.backwardError=null;}
    if(time===0){rawMetric.residual=stableMetric.residual=null;rawMetric.backwardError=stableMetric.backwardError=null;}
    for(let i=0;i<raw.length;i++){lo=Math.min(lo,raw[i],stable[i],reference?.[i]??Infinity);hi=Math.max(hi,raw[i],stable[i],reference?.[i]??-Infinity);differenceMax=Math.max(differenceMax,Math.abs(stable[i]-raw[i]));}
    frames.push({time,dt,afc:afcStats,raw,stable,reference,F:forces.F,effectiveF:forces.effectiveF,rawRhs,stableRhs,metrics:{raw:rawMetric,stable:stableMetric}});
  }
  snapshot(0,Math.min(state.dt,state.duration));
  for(let step=1;step<=steps;step++){
    const t=Math.min(state.duration,step*state.dt),prevT=(step-1)*state.dt,dt=t-prevT,pair=stepSystems(dt);
    const nextG=evaluateBoundary(t),load=sourceVaries?sourceLoad(mesh,source,state,t):system.load,nextForces=reducedForcing(system,state,load,nextG,sourceVaries?supgLoad(mesh,source,state,t,system.supgElements):system.supgF);
    const rawRhs=multiplyRows(pair.raw.rightRows,rawB),stableRhs=identical?rawRhs: multiplyRows(pair.stable.rightRows,stableB);
    for(let i=0;i<N;i++){
      rawRhs[i]+=dt*(theta*nextForces.F[i]+(1-theta)*forces.F[i]);
      if(!identical&&!afcStepper)stableRhs[i]+=dt*(theta*nextForces.effectiveF[i]+(1-theta)*forces.effectiveF[i]);
    }
    // Consistent mass has an interior-boundary block; full row-sum lumping does not.
    for(const[row,id,,,m,,h]of system.boundaryTerms){const change=nextG[id]-g[id];rawRhs[row]-=m*change;if(!identical&&!afcStepper)stableRhs[row]-=((state.massType==='consistent'?m:0)+h)*change;}
    rawB=solveFactored(pair.raw.factor,rawRhs);
    if(afcStepper){const advanced=afcStepper.advance(stableB,g,prevT,dt);stableB=advanced.b;afcStats=advanced.stats;nextForces.effectiveF=afcReducedLoad(afcStepper.graph,load,nextG);}else stableB=identical?rawB:solveFactored(pair.stable.factor,stableRhs);
    g=nextG;forces=nextForces;lastPair=pair;if(step%saveStride===0||step===steps)snapshot(t,dt,rawRhs,afcStepper?null:stableRhs,pair);
    if(step%20===0||step===steps)progress(step,steps);
  }
  const firstPair=stepSystems(Math.min(state.dt,state.duration)),speed=Math.hypot(...beta);
  return{...system,afc:afcStats,afcD:afcStepper?.graph.D,ML,comparisonMass,effectiveA,frames,steps:frames.length-1,computedSteps:steps,saveStride,theta,initialMismatch,displayRange:[lo,hi],differenceMax,
    timeMatrix:firstPair.raw.B,stableTimeMatrix:afcStepper?null:firstPair.stable.B,lastTimeMatrix:lastPair.raw.B,lastStableTimeMatrix:afcStepper?null:lastPair.stable.B,
    pe:eps>0?speed/(2*state.n*eps):speed>0?Infinity:0,
    effectivePe:eps+addedDiffusion(state)>0?speed/(2*state.n*(eps+addedDiffusion(state))):speed>0?Infinity:0,
    courant:speed*state.dt*state.n,diffusionNumber:eps*state.dt*state.n**2,ms:performance.now()-started};
}
