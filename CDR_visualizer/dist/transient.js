import {assemble,expression,boundaryEvaluator,sourceLoad,diagnostics,exactSolution,physicalDiffusion,velocity,isTransient} from './fem.js';

export function matrixRows(matrix,n){return Array.from({length:n},(_,i)=>{const row=[];for(let j=0;j<n;j++)if(matrix[i*n+j]!==0)row.push([j,matrix[i*n+j]]);return row;});}
export function multiplyRows(rows,b){return Float64Array.from(rows,row=>row.reduce((sum,[j,v])=>sum+v*b[j],0));}

// Factor once per distinct step size, reuse throughout the time evolution.
export function factorDense(matrix,n){
  const a=matrix.slice(),pivots=new Int32Array(n);let norm=0;
  for(const v of a)norm=Math.max(norm,Math.abs(v));
  const threshold=norm*Number.EPSILON*n*8;
  for(let k=0;k<n;k++){
    let p=k;for(let i=k+1;i<n;i++)if(Math.abs(a[i*n+k])>Math.abs(a[p*n+k]))p=i;
    if(!Number.isFinite(a[p*n+k])||Math.abs(a[p*n+k])<=threshold)throw Error('The time-step matrix is singular or numerically ill-conditioned. Restore matrix entries or change the time step.');
    pivots[k]=p;if(p!==k)for(let j=0;j<n;j++){const v=a[k*n+j];a[k*n+j]=a[p*n+j];a[p*n+j]=v;}
    for(let i=k+1;i<n;i++){const f=a[i*n+k]/a[k*n+k];a[i*n+k]=f;if(f!==0)for(let j=k+1;j<n;j++)a[i*n+j]-=f*a[k*n+j];}
  }
  const lower=[],upper=[],diagonal=new Float64Array(n);
  for(let i=0;i<n;i++){lower[i]=[];upper[i]=[];diagonal[i]=a[i*n+i];for(let j=0;j<n;j++)if(a[i*n+j]!==0){if(j<i)lower[i].push([j,a[i*n+j]]);else if(j>i)upper[i].push([j,a[i*n+j]]);}}
  return{lower,upper,diagonal,pivots};
}
export function solveFactored(factor,rhs){
  const{lower,upper,diagonal,pivots}=factor,b=rhs.slice(),n=b.length;
  for(let k=0;k<n;k++){const p=pivots[k],v=b[k];b[k]=b[p];b[p]=v;}
  for(let i=0;i<n;i++)for(const[j,v]of lower[i])b[i]-=v*b[j];
  for(let i=n-1;i>=0;i--){for(const[j,v]of upper[i])b[i]-=v*b[j];b[i]/=diagonal[i];if(!Number.isFinite(b[i]))throw Error('The time evolution produced a non-finite solution. Check data and matrix edits.');}
  return b;
}

function reducedForcing(system,state,load,g){
  const F=load.slice(),lift=new Float64Array(system.N);
  for(const[row,id,a,k]of system.boundaryTerms){F[row]-=a*g[id];lift[row]-=k*g[id];}
  for(const[i,v]of Object.entries(state.rhsEdits))F[+i]=v;
  return{F,effectiveF:F.map((v,i)=>v+state.diffusion*lift[i])};
}
export function computeTrajectory(state,progress=()=>{}){
  const started=performance.now();
  if(!isTransient(state))throw Error('Choose a time-dependent equation.');
  if(!Number.isFinite(state.dt)||state.dt<=0||!Number.isFinite(state.duration)||state.duration<=0)throw Error('Time step and end time must be positive.');
  if(!['be','cn'].includes(state.scheme)||!['consistent','lumped'].includes(state.massType))throw Error('Choose a supported time integrator and mass matrix.');
  const steps=Math.max(1,Math.ceil(state.duration/state.dt-1e-12));
  if(steps>240)throw Error('This interactive run supports at most 240 steps. Increase Δt or reduce the end time.');
  const system=assemble(state),{N,A,K,M,lumped,mesh}=system,theta=state.scheme==='be'?1:.5;
  const ML=new Float64Array(N*N);for(let i=0;i<N;i++)ML[i*N+i]=lumped[i];
  const comparisonMass=state.massType==='lumped'?ML:M,effectiveA=A.map((v,k)=>v+state.diffusion*K[k]);
  const identical=state.diffusion===0&&state.massType==='consistent',factorCache=new Map();
  function stepSystems(dt){
    const key=dt.toPrecision(12);if(factorCache.has(key))return factorCache.get(key);
    const make=(mass,operator)=>{
      const B=mass.map((v,k)=>v+theta*dt*operator[k]),R=mass.map((v,k)=>v-(1-theta)*dt*operator[k]);
      return{B,rows:matrixRows(B,N),rightRows:matrixRows(R,N),factor:factorDense(B,N)};
    };
    const raw=make(M,A),stable=identical?raw:make(comparisonMass,effectiveA),pair={raw,stable};factorCache.set(key,pair);return pair;
  }
  let g=boundaryEvaluator(state,mesh)(0);const evaluateBoundary=boundaryEvaluator(state,mesh),source=expression(state.source),initial=expression(state.initial),eps=physicalDiffusion(state),beta=velocity(state);
  let rawB=new Float64Array(N),stableB;
  for(let i=0;i<N;i++){const[x,y]=mesh.nodes[mesh.interior[i]];rawB[i]=initial(x,y,eps,0,...beta);}stableB=rawB.slice();
  let initialMismatch=0;for(const id of mesh.boundary){const[x,y]=mesh.nodes[id];if(Math.abs(initial(x,y,eps,0,...beta)-g[id])>1e-8)initialMismatch++;}
  const sourceVaries=/\bt\b/.test(state.source),exact=exactSolution(state);
  let forces=reducedForcing(system,state,system.load,g),lastPair=null;
  const frames=[];let lo=0,hi=0,differenceMax=0;
  function snapshot(time,dt,rawRhs=null,stableRhs=null,pair=null){
    const raw=g.slice(),stable=g.slice();mesh.interior.forEach((id,i)=>{raw[id]=rawB[i];stable[id]=stableB[i];});
    const reference=exact?Float64Array.from(mesh.nodes,([x,y])=>exact(x,y,time)):null,exactNow=exact?(x,y)=>exact(x,y,time):null;
    const rawMetric=diagnostics(pair?.raw.B||M,rawRhs||multiplyRows(matrixRows(M,N),rawB),rawB,raw,mesh,exactNow,pair?.raw.rows);
    const stableMetric=identical?{...rawMetric}:diagnostics(pair?.stable.B||comparisonMass,stableRhs||multiplyRows(matrixRows(comparisonMass,N),stableB),stableB,stable,mesh,exactNow,pair?.stable.rows);
    if(time===0){rawMetric.residual=stableMetric.residual=null;rawMetric.backwardError=stableMetric.backwardError=null;}
    for(let i=0;i<raw.length;i++){lo=Math.min(lo,raw[i],stable[i],reference?.[i]??Infinity);hi=Math.max(hi,raw[i],stable[i],reference?.[i]??-Infinity);differenceMax=Math.max(differenceMax,Math.abs(stable[i]-raw[i]));}
    frames.push({time,dt,raw,stable,reference,F:forces.F,effectiveF:forces.effectiveF,rawRhs,stableRhs,metrics:{raw:rawMetric,stable:stableMetric}});
  }
  snapshot(0,Math.min(state.dt,state.duration));
  for(let step=1;step<=steps;step++){
    const t=Math.min(state.duration,step*state.dt),prevT=frames[step-1].time,dt=t-prevT,pair=stepSystems(dt);
    const nextG=evaluateBoundary(t),load=sourceVaries?sourceLoad(mesh,source,state,t):system.load,nextForces=reducedForcing(system,state,load,nextG);
    const rawRhs=multiplyRows(pair.raw.rightRows,rawB),stableRhs=identical?rawRhs: multiplyRows(pair.stable.rightRows,stableB);
    for(let i=0;i<N;i++){
      rawRhs[i]+=dt*(theta*nextForces.F[i]+(1-theta)*forces.F[i]);
      if(!identical)stableRhs[i]+=dt*(theta*nextForces.effectiveF[i]+(1-theta)*forces.effectiveF[i]);
    }
    // Consistent mass has an interior-boundary block; full row-sum lumping does not.
    for(const[row,id,,,m]of system.boundaryTerms){const correction=m*(nextG[id]-g[id]);rawRhs[row]-=correction;if(!identical&&state.massType==='consistent')stableRhs[row]-=correction;}
    rawB=solveFactored(pair.raw.factor,rawRhs);stableB=identical?rawB:solveFactored(pair.stable.factor,stableRhs);
    g=nextG;forces=nextForces;lastPair=pair;snapshot(t,dt,rawRhs,stableRhs,pair);
    if(step%20===0||step===steps)progress(step,steps);
  }
  const firstPair=stepSystems(Math.min(state.dt,state.duration)),speed=Math.hypot(...beta);
  return{...system,ML,effectiveA,frames,steps,theta,initialMismatch,displayRange:[lo,hi],differenceMax,
    timeMatrix:firstPair.raw.B,stableTimeMatrix:firstPair.stable.B,lastTimeMatrix:lastPair.raw.B,lastStableTimeMatrix:lastPair.stable.B,
    pe:eps>0?speed/(2*state.n*eps):speed>0?Infinity:0,
    effectivePe:eps+state.diffusion>0?speed/(2*state.n*(eps+state.diffusion)):speed>0?Infinity:0,
    courant:speed*state.dt*state.n,diffusionNumber:eps*state.dt*state.n**2,ms:performance.now()-started};
}
