import {physicalDiffusion,velocity,expression,sourceLoad,boundaryEvaluator} from './fem.js';
import {factorDense,solveFactored} from './linear.js';
export const isAfc = state => ['afc','afc-low'].includes(state.stabilization);

// Assemble before Dirichlet elimination: graph viscosity needs both a_ij and a_ji,
// including boundary edges. No inferred transpose of a reduced boundary block.
export function buildAfc(state,system){
  if(Object.keys(state.matrixEdits||{}).length||Object.keys(state.rhsEdits||{}).length)throw Error('AFC needs the assembled PDE operator. Clear matrix and load overrides first.');
  if((state.reaction??0)<0)throw Error('This AFC implementation supports c ≥ 0. Choose another method for reaction growth.');
  const {mesh,N}=system,size=mesh.nodes.length,rows=Array.from({length:size},()=>new Map()),mass=new Float64Array(size),beta=velocity(state),eps=physicalDiffusion(state);
  const add=(i,j,v)=>rows[i].set(j,(rows[i].get(j)||0)+v);
  for(const tri of mesh.triangles){
    const p=tri.map(id=>mesh.nodes[id]),twice=(p[1][0]-p[0][0])*(p[2][1]-p[0][1])-(p[2][0]-p[0][0])*(p[1][1]-p[0][1]),area=Math.abs(twice)/2;
    const grad=p.map((_,i)=>[(p[(i+1)%3][1]-p[(i+2)%3][1])/twice,(p[(i+2)%3][0]-p[(i+1)%3][0])/twice]);
    for(let i=0;i<3;i++){mass[tri[i]]+=area/3;for(let j=0;j<3;j++)add(tri[i],tri[j],eps*area*(grad[i][0]*grad[j][0]+grad[i][1]*grad[j][1])+area/3*(beta[0]*grad[j][0]+beta[1]*grad[j][1])+(state.reaction??0)*area*(i===j?1/6:1/12));}
  }
  const edges=[],neighbors=rows.map((r,i)=>[...r.keys()].filter(j=>j!==i)),lowRows=rows.map(r=>new Map(r));
  for(let i=0;i<size;i++)for(const[j,a]of rows[i])if(j>i){const viscosity=Math.max(0,a,rows[j].get(i)||0);if(viscosity>0){edges.push({i,j,viscosity});lowRows[i].set(j,a-viscosity);lowRows[j].set(i,(rows[j].get(i)||0)-viscosity);lowRows[i].set(i,(lowRows[i].get(i)||0)+viscosity);lowRows[j].set(j,(lowRows[j].get(j)||0)+viscosity);}}
  const D=new Float64Array(N*N),lowOperator=new Float64Array(N*N);let maxStep=Infinity;
  for(const id of mesh.interior){const i=mesh.dof[id],diag=lowRows[id].get(id)||0;if(diag>0)maxStep=Math.min(maxStep,.8*mass[id]/diag);for(const[j,v]of lowRows[id]){const col=mesh.dof[j];if(col>=0){lowOperator[i*N+col]=v;D[i*N+col]=v-system.A[i*N+col];}}}
  return {mesh,N,size,rows,lowRows,edges,neighbors,mass,D,lowOperator,maxStep};
}
export function afcReducedLoad(graph,load,g){
  const rhs=load.slice();for(const id of graph.mesh.interior){const i=graph.mesh.dof[id];for(const[j,v]of graph.lowRows[id])if(graph.mesh.dof[j]<0)rhs[i]-=v*g[j];}return rhs;
}
function fullField(graph,b,g){const u=g.slice();graph.mesh.interior.forEach((id,i)=>u[id]=b[i]);return u;}

// Symmetric BJK-type local-extremum limiter, mu=1. Low-order matrix is fixed.
export function steadyFluxes(graph,u){
  const {size,edges,neighbors,mesh}=graph,pplus=new Float64Array(size),pminus=new Float64Array(size),q=new Float64Array(size),rplus=new Float64Array(size).fill(1),rminus=new Float64Array(size).fill(1),flux=new Float64Array(edges.length);
  edges.forEach(({i,j,viscosity:g},k)=>{const f=g*(u[i]-u[j]);flux[k]=f;pplus[i]+=Math.max(0,f);pminus[i]+=Math.min(0,f);pplus[j]+=Math.max(0,-f);pminus[j]+=Math.min(0,-f);q[i]+=g;q[j]+=g;});
  for(const i of mesh.interior){let lo=u[i],hi=u[i];for(const j of neighbors[i]){lo=Math.min(lo,u[j]);hi=Math.max(hi,u[j]);}if(pplus[i]>0)rplus[i]=Math.min(1,q[i]*(hi-u[i])/pplus[i]);if(pminus[i]<0)rminus[i]=Math.min(1,q[i]*(lo-u[i])/pminus[i]);}
  const correction=new Float64Array(graph.N),alpha=new Float64Array(edges.length);let limited=0;
  edges.forEach(({i,j},k)=>{const f=flux[k],a=f>=0?Math.min(rplus[i],rminus[j]):Math.min(rminus[i],rplus[j]);alpha[k]=a;if(a<1-1e-12)limited++;const row=mesh.dof[i],col=mesh.dof[j];if(row>=0)correction[row]+=a*f;if(col>=0)correction[col]-=a*f;});
  return {correction,alpha,limitedFraction:limited/(edges.length||1)};
}
export function solveAfcSteady(state,system){
  const started=performance.now(),graph=buildAfc(state,system),{N,g}=system,lowRhs=afcReducedLoad(graph,system.load,g),factor=factorDense(graph.lowOperator,N);
  const maxIterations=state.afcMaxIterations??500,tolerance=state.afcTolerance??1e-8,relax=state.afcRelaxation??.5;
  if(!Number.isInteger(maxIterations)||maxIterations<1||maxIterations>2000)throw Error('AFC iteration limit must be an integer from 1 to 2000.');
  if(!Number.isFinite(tolerance)||tolerance<1e-12||tolerance>1e-3)throw Error('AFC tolerance must be between 1e-12 and 1e-3.');
  if(!Number.isFinite(relax)||relax<=0||relax>1)throw Error('AFC relaxation must be in (0,1].');
  let b=solveFactored(factor,lowRhs),iterations=0,converged=false,nonlinearResidual=Infinity,limited;
  for(;;){
    limited=state.stabilization==='afc-low'?{correction:new Float64Array(N),alpha:new Float64Array(graph.edges.length),limitedFraction:1}:steadyFluxes(graph,fullField(graph,b,g));
    const rhs=lowRhs.map((v,i)=>v+limited.correction[i]);let norm=0,scale=0;
    for(let i=0;i<N;i++){let v=-rhs[i],row=0;for(let j=0;j<N;j++){const a=graph.lowOperator[i*N+j];v+=a*b[j];row+=Math.abs(a*b[j]);}norm=Math.max(norm,Math.abs(v));scale=Math.max(scale,row+Math.abs(rhs[i]));}
    nonlinearResidual=norm/(scale||1);if(nonlinearResidual<=tolerance){converged=true;break;}if(iterations>=maxIterations)break;
    const next=solveFactored(factor,rhs);b=b.map((v,i)=>(1-relax)*v+relax*next[i]);iterations++;
  }
  // Report the actual final nonlinear operator, not the frozen low-order solve.
  const effectiveA=system.A.slice(),effectiveF=system.F.slice();
  graph.edges.forEach(({i,j,viscosity},k)=>{const v=(1-limited.alpha[k])*viscosity,ri=graph.mesh.dof[i],rj=graph.mesh.dof[j];if(ri>=0){effectiveA[ri*N+ri]+=v;if(rj>=0)effectiveA[ri*N+rj]-=v;else effectiveF[ri]+=v*g[j];}if(rj>=0){effectiveA[rj*N+rj]+=v;if(ri>=0)effectiveA[rj*N+ri]-=v;else effectiveF[rj]+=v*g[i];}});
  return {stableB:b,effectiveA,effectiveF,afcD:graph.D,afc:{iterations,converged,nonlinearResidual,limitedFraction:limited.limitedFraction,ms:performance.now()-started,kind:'steady',edgeCount:graph.edges.length}};
}

// Explicit Euler FCT stage: lumped target mass, graph-viscosity predictor,
// conservative pair fluxes, and Zalesak bounds. No solution-value clipping.
export function fctStage(graph,u,load,dt,nextG,correct=true){
  const {mesh,size,mass,lowRows,edges,neighbors}=graph,low=u.slice();
  for(const id of mesh.interior){let rate=load[mesh.dof[id]];for(const[j,v]of lowRows[id])rate-=v*u[j];low[id]+=dt*rate/mass[id];}
  mesh.boundary.forEach(id=>low[id]=nextG[id]);
  if(!correct)return {u:low,limitedFraction:1,correctionMass:0};
  const plus=new Float64Array(size),minus=new Float64Array(size),rp=new Float64Array(size).fill(1),rm=new Float64Array(size).fill(1),flux=new Float64Array(edges.length);
  edges.forEach(({i,j,viscosity},k)=>{let f=dt*viscosity*(u[i]-u[j]);if(f*(low[i]-low[j])<0)f=0;flux[k]=f;plus[i]+=Math.max(0,f);minus[i]+=Math.min(0,f);plus[j]+=Math.max(0,-f);minus[j]+=Math.min(0,-f);});
  for(const i of mesh.interior){let lo=Math.min(u[i],low[i]),hi=Math.max(u[i],low[i]);for(const j of neighbors[i]){lo=Math.min(lo,u[j],low[j]);hi=Math.max(hi,u[j],low[j]);}if(plus[i]>0)rp[i]=Math.min(1,mass[i]*(hi-low[i])/plus[i]);if(minus[i]<0)rm[i]=Math.min(1,mass[i]*(lo-low[i])/minus[i]);}
  const out=low.slice();let limited=0,correctionMass=0;
  edges.forEach(({i,j},k)=>{const f=flux[k],a=f>=0?Math.min(rp[i],rm[j]):Math.min(rm[i],rp[j]);if(a<1-1e-12)limited++;if(mesh.dof[i]>=0){out[i]+=a*f/mass[i];correctionMass+=a*f;}if(mesh.dof[j]>=0){out[j]-=a*f/mass[j];correctionMass-=a*f;}});
  return {u:out,limitedFraction:limited/(edges.length||1),correctionMass};
}
export function afcTimeStepper(state,system){
  const graph=buildAfc(state,system),source=expression(state.source),boundary=boundaryEvaluator(state,system.mesh),varies=/\bt\b/.test(state.source);let totalSubsteps=0,totalMs=0;
  return {graph,advance(b,g,t,dt){
    const started=performance.now(),substeps=Math.max(1,Math.ceil(dt/graph.maxStep-1e-12));if(!Number.isFinite(substeps)||totalSubsteps+substeps>50000)throw Error('AFC needs too many internal substeps. Reduce diffusion/reaction or mesh size.');
    let u=fullField(graph,b,g),limitedFraction=0,dataNonnegative=u.every(v=>v>=-1e-12);
    for(let k=0;k<substeps;k++){
      const time=t+k*dt/substeps,h=dt/substeps,nextG=boundary(time+h),load=varies?sourceLoad(system.mesh,source,state,time):system.load,nextLoad=varies?sourceLoad(system.mesh,source,state,time+h):system.load;
      dataNonnegative&&=load.every(v=>v>=-1e-12)&&nextLoad.every(v=>v>=-1e-12)&&nextG.every(v=>v>=-1e-12);
      const one=fctStage(graph,u,load,h,nextG,state.stabilization==='afc'),two=fctStage(graph,one.u,nextLoad,h,nextG,state.stabilization==='afc');
      for(const id of graph.mesh.interior)u[id]=.5*(u[id]+two.u[id]);graph.mesh.boundary.forEach(id=>u[id]=nextG[id]);limitedFraction+=.5*(one.limitedFraction+two.limitedFraction)/substeps;
    }
    totalSubsteps+=substeps;totalMs+=performance.now()-started;
    return {b:Float64Array.from(graph.mesh.interior,id=>u[id]),stats:{kind:'transient',substeps,totalSubsteps,limiterPasses:state.stabilization==='afc'?2*substeps:0,limitedFraction,dataNonnegative,iterations:0,converged:true,ms:totalMs,maxInternalDt:graph.maxStep,edgeCount:graph.edges.length}};
  }};
}
