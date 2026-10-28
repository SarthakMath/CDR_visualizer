import test from 'node:test';
import assert from 'node:assert/strict';
import {matrix,set,get,add,combine,entries,isSparse} from '../dist/matrix.js';
import {factorDense,solveFactored} from '../dist/linear.js';
import {preset,assemble,compute,makeMesh} from '../dist/fem.js';
import {computeTrajectory} from '../dist/transient.js';
function sparseCopy(a,n){const s=matrix(n,true);for(let i=0;i<n;i++)for(let j=0;j<n;j++)set(s,i,j,a[i*n+j]);return s;}
test('sparse Krylov solutions agree with pivoted LU for diffusion, transport, and SUPG time systems',()=>{
 for(const name of ['poisson','steady-adr','afc-layer','solid-rotation','supg-front']){
  const state={...preset(name),n:8},a=assemble(state),N=a.N;
  const A=['solid-rotation','supg-front'].includes(name)?combine([a.M,1],[a.H,1],[a.A,.02],[a.S,.02]):a.A;
  const rhs=Float64Array.from({length:N},(_,i)=>Math.sin(i+1)+.1);
  const dense=solveFactored(factorDense(A,N),rhs),factor=factorDense(sparseCopy(A,N),N),sparse=solveFactored(factor,rhs);
  assert.ok(Math.max(...sparse.map((v,i)=>Math.abs(v-dense[i])))<1e-7,name);
  const rhs2=rhs.map(v=>v*1.01),again=solveFactored(factor,rhs2);
  assert.ok(Math.max(...again.map((v,i)=>Math.abs(v-1.01*dense[i])))<1e-7,name);
 }
});
test('sparse storage preserves nonlocal edits, cloning, and export through structured clone',()=>{
 const a=matrix(1000);set(a,1,999,3);add(a,1,999,2);set(a,4,4,7);
 const b=structuredClone(combine([a,2]));assert.ok(isSparse(b));assert.equal(get(b,1,999),10);assert.equal(get(b,0,999),0);assert.equal(entries(b).length,2);assert.equal(get(a,1,999),5);
});
test('128 mesh uses sparse matrices and maintains second-order Poisson convergence',()=>{
 const a=compute({...preset('poisson'),n:64}),b=compute({...preset('poisson'),n:128});
 assert.equal(b.N,16129);assert.equal(b.mesh.triangles.length,32768);assert.ok(isSparse(b.A));
 assert.ok(entries(b.A).length<7*b.N);assert.ok(a.metrics.raw.l2/b.metrics.raw.l2>3.9);assert.ok(b.metrics.raw.backwardError<1e-10);
 assert.throws(()=>makeMesh(129));
});
test('large transient runs preserve time steps while thinning saved frames and retaining the endpoint',()=>{
 const s={...preset('rotation-diffusion'),n:65,duration:.081,dt:.001};
 const r=computeTrajectory(s);assert.equal(r.computedSteps,81);assert.equal(r.saveStride,2);assert.equal(r.frames.length,42);assert.equal(r.steps,41);assert.equal(r.frames.at(-1).time,s.duration);
 assert.ok(isSparse(r.M));assert.ok(r.frames.at(-1).stable.every(Number.isFinite));assert.ok(Math.min(...r.frames.at(-1).stable)>-1e-10);assert.equal(r.frames.at(-1).reference,null);
});
test('sparse solves reject an inconsistent singular system instead of showing an unconverged result',()=>{
 const a=matrix(2,true);set(a,0,0,1);set(a,0,1,1);set(a,1,0,1);set(a,1,1,1);
 assert.throws(()=>solveFactored(factorDense(a,2),new Float64Array([1,2])),/did not converge/);
});
