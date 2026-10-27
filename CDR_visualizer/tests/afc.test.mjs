import test from 'node:test';
import assert from 'node:assert/strict';
import {preset,assemble,compute} from '../dist/fem.js';
import {computeTrajectory} from '../dist/transient.js';
import {buildAfc,fctStage,steadyFluxes} from '../dist/afc.js';
const near=(a,b,tol=1e-10)=>assert.ok(Math.abs(a-b)<tol,`${a} differs from ${b}`);
test('full graph viscosity is symmetric, has zero row sums, and yields nonpositive low-order off diagonals',()=>{
 const s=preset('afc-front');s.n=4;const sys=assemble(s),g=buildAfc(s,sys);
 for(const e of g.edges){near(g.rows[e.i].get(e.j)-g.lowRows[e.i].get(e.j),e.viscosity);near(g.rows[e.j].get(e.i)-g.lowRows[e.j].get(e.i),e.viscosity);}
 for(let i=0;i<g.size;i++){let before=0,after=0;for(const v of g.rows[i].values())before+=v;for(const[j,v]of g.lowRows[i]){after+=v;if(i!==j)assert.ok(v<1e-14);}near(before,after);}
 // This verifies boundary graph rows too, not just reduced I-I symmetry.
 assert.ok(g.edges.some(e=>g.mesh.dof[e.i]<0&&g.mesh.dof[e.j]>=0));
});
test('FCT pair corrections preserve weighted mass on a closed two-node graph and stay within local bounds',()=>{
 const mesh={interior:[0,1],boundary:[],dof:Int32Array.from([0,1])},g={mesh,size:2,N:2,mass:Float64Array.from([1,2]),neighbors:[[1],[0]],edges:[{i:0,j:1,viscosity:1}],lowRows:[new Map([[0,1],[1,-1]]),new Map([[0,-1],[1,1]])]};
 const u=Float64Array.from([.9,.1]),load=new Float64Array(2),r=fctStage(g,u,load,.1,new Float64Array(2));
 near(r.u[0]+2*r.u[1],1.1);near(r.correctionMass,0);assert.ok(r.u.every(v=>v>=.1-1e-12&&v<=.9+1e-12));
 const low=fctStage(g,u,load,.1,new Float64Array(2),false);assert.ok(r.u[0]>low.u[0]);
});
test('nonlinear steady limiter has symmetric bounded factors and pairwise conservative internal corrections',()=>{
 const s=preset('afc-layer');s.n=4;const g=buildAfc(s,assemble(s)),u=Float64Array.from(g.mesh.nodes,([x,y])=>x*(1-x)*y*(1-y));
 const result=steadyFluxes(g,u);assert.ok(result.alpha.every(v=>v>=0&&v<=1));
 // Boundary exchange excluded by using zero flux on edges touching boundary.
 const closed={...g,edges:g.edges.filter(e=>g.mesh.dof[e.i]>=0&&g.mesh.dof[e.j]>=0)};
 const c=steadyFluxes(closed,u).correction;near(c.reduce((a,b)=>a+b,0),0);
});
test('steady AFC converges to its nonlinear equation, improves layer overshoot, and flags iteration limits',()=>{
 const s=preset('afc-layer'),r=compute(s);assert.ok(r.afc.converged);assert.ok(r.afc.iterations>1);assert.ok(r.afc.nonlinearResidual<=s.afcTolerance);
 assert.ok(r.metrics.stable.max<1.01);assert.ok(r.metrics.stable.l2<r.metrics.raw.l2);assert.ok(r.metrics.stable.min>=-1e-10);
 s.afcMaxIterations=1;const incomplete=compute(s);assert.equal(incomplete.afc.converged,false);assert.ok(incomplete.afc.nonlinearResidual>s.afcTolerance);
});
test('FCT moving front remains nonnegative, bounded, and sharper than low order; baseline does not change',()=>{
 const s=preset('afc-front'),r=computeTrajectory(s);s.stabilization='afc-low';const low=computeTrajectory(s);
 for(let k=0;k<r.frames.length;k++){const f=r.frames[k];assert.ok(f.metrics.stable.min>=-1e-12);assert.ok(f.metrics.stable.max<=1+1e-10);assert.deepEqual(f.raw,low.frames[k].raw);assert.equal(f.metrics.stable.backwardError,null);}
 assert.ok(r.frames.at(-1).metrics.raw.min<-.2);assert.ok(r.frames.at(-1).metrics.stable.max>low.frames.at(-1).metrics.stable.max);
 assert.ok(r.frames.at(-1).metrics.stable.l2<low.frames.at(-1).metrics.stable.l2);assert.equal(r.stableTimeMatrix,null);
});
test('automatic substeps preserve constants with time-varying boundaries and forcing in every pattern',()=>{
 for(const pattern of ['forward','backward','alternating']){
 const s=preset('afc-front');Object.assign(s,{n:4,pattern,mode:'adr',epsilon:.1,reaction:2,initial:'1',source:'1+c*(1+t)',dt:.2,duration:.4});
 s.boundaries=Object.fromEntries(['left','right','bottom','top'].map(k=>[k,'1+t']));
 const r=computeTrajectory(s);assert.ok(r.afc.totalSubsteps>r.steps);
 // SSPRK2 has O(h^3) local error for reaction with affine source unless stages are exact; this affine solution is reproduced exactly.
 for(const f of r.frames)for(const v of f.stable)near(v,1+f.time);
 }
});
test('signed source data are reported rather than clipped, and unsupported reaction/edits are explicit',()=>{
 const s=preset('afc-front');Object.assign(s,{n:4,initial:'0',source:'-1',duration:.01});s.boundaries={left:'0',right:'0',bottom:'0',top:'0'};
 const r=computeTrajectory(s);assert.ok(r.frames.at(-1).metrics.stable.min<0);assert.equal(r.afc.dataNonnegative,false);
 s.reaction=-1;assert.throws(()=>computeTrajectory(s),/c ≥ 0/);s.reaction=0;s.matrixEdits={'0,0':1};assert.throws(()=>computeTrajectory(s),/overrides/);
 const t=preset('afc-layer');t.afcTolerance=0;assert.throws(()=>compute(t),/tolerance/);
});

test('steady μ controls flux budgets and reproduces the shifted affine patch above its field-specific threshold',()=>{
 const s={...preset('afc-affine'),pattern:'shifted'},half=compute({...s,afcMu:.5}),one=compute(s),zero=compute({...s,afcMu:0});
 for(const r of [half,one,zero]){assert.ok(r.afc.converged);assert.ok(r.afc.affine.baselineMaxNodalError<1e-12);}
 assert.ok(zero.afc.affine.maxNodalError>.05);assert.ok(half.afc.affine.maxNodalError>.02);assert.ok(one.afc.affine.maxNodalError<1e-8);
 assert.ok(half.afc.affineLimiter.retainedResidual>.005);near(one.afc.affineLimiter.retainedResidual,0,1e-14);
 const threshold=one.afc.affineLimiter.requiredMu;assert.ok(threshold>.73&&threshold<.74);
 const at=compute({...s,afcMu:threshold+1e-10});assert.ok(at.afc.affine.maxNodalError<1e-8);
 // At a fixed field, increasing μ never reduces an edge factor.
 const graph=buildAfc(s,assemble(s)),u=Float64Array.from(graph.mesh.nodes,([x,y])=>1+x+y),a=steadyFluxes(graph,u,.5),b=steadyFluxes(graph,u,1);
 a.alpha.forEach((v,i)=>assert.ok(b.alpha[i]>=v-1e-14));
});
test('uniform mesh symmetry can preserve the affine patch despite active limiting',()=>{
 const r=compute({...preset('afc-affine'),pattern:'forward',afcMu:0});
 assert.ok(r.afc.affine.maxNodalError<1e-12);assert.ok(r.afc.limitedFraction>.5);near(r.afc.affineLimiter.retainedResidual,0,1e-14);
});
test('shifted mesh quadrature and cross-section interpolation preserve affine fields',async()=>{
 const {makeMesh,interpolate,sourceLoad,diagnostics,expression}=await import('../dist/fem.js');
 for(const n of [2,3,8,16]){
 const mesh=makeMesh(n,'shifted'),u=Float64Array.from(mesh.nodes,([x,y])=>1+x+y);
 for(let iy=0;iy<=20;iy++)for(let ix=0;ix<=20;ix++)near(interpolate(mesh,u,ix/20,iy/20),1+ix/20+iy/20,1e-12);
 const r=compute({...preset('afc-affine'),n,pattern:'shifted'});assert.ok(r.afc.affine.baselineMaxNodalError<1e-12);
 const metric=diagnostics(new Float64Array(1),new Float64Array(1),new Float64Array(1),u,mesh,(x,y)=>1+x+y);near(metric.integral,2);assert.ok(metric.l2<1e-14);
 // The unit load is the geometric lumped interior mass, not a uniform-area approximation.
 const state={...preset('afc-affine'),n,pattern:'shifted'},sys=assemble(state),load=sourceLoad(mesh,expression('1'),state);
 load.forEach((v,i)=>near(v,sys.lumped[i],1e-14));
 }
});
test('steady μ rejects invalid values and is absent from transient limiter behavior',()=>{
 for(const afcMu of [-1,NaN,Infinity,10.1])assert.throws(()=>compute({...preset('afc-affine'),afcMu}),/μ/);
 const s={...preset('afc-front'),n:4,duration:.01};
 assert.deepEqual(computeTrajectory({...s,afcMu:0}).frames.at(-1).stable,computeTrajectory({...s,afcMu:10}).frames.at(-1).stable);
});

test('μ = 1 leaves affine fluxes unlimited on every original centrally symmetric mesh pattern',()=>{
 for(const pattern of ['forward','backward','alternating'])for(const n of [2,5,8]){
 const s={...preset('afc-affine'),pattern,n,reaction:2},graph=buildAfc(s,assemble(s));
 for(const i of graph.mesh.interior){const [x,y]=graph.mesh.nodes[i];for(const j of graph.neighbors[i]){
   const [a,b]=graph.mesh.nodes[j];assert.ok(graph.neighbors[i].some(k=>Math.abs(graph.mesh.nodes[k][0]-(2*x-a))<1e-12&&Math.abs(graph.mesh.nodes[k][1]-(2*y-b))<1e-12));
 }}
 for(const [a,b,c] of [[1,0,0],[0,1,0],[1,-1,2],[2,3,-1],[-3,1,4],[0,0,1]]){
 const u=Float64Array.from(graph.mesh.nodes,([x,y])=>a*x+b*y+c),limited=steadyFluxes(graph,u,1);
 graph.edges.forEach((e,k)=>{if(graph.mesh.dof[e.i]>=0||graph.mesh.dof[e.j]>=0)near(limited.alpha[k],1,1e-12);});
 }
 }
});
