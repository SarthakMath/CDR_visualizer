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
