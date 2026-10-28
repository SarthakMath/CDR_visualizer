import test from 'node:test';
import assert from 'node:assert/strict';
import {preset,rotationProfile,velocityAt,exactSolution,assemble,expression} from '../dist/fem.js';
import {buildAfc} from '../dist/afc.js';
import {computeTrajectory} from '../dist/transient.js';
test('paper profiles, inverse rotation, and one-period return',()=>{
 const s=preset('solid-rotation'),u=exactSolution(s);
 assert.equal(rotationProfile(.5,.75),0);assert.equal(rotationProfile(.5,.87),1);
 assert.equal(rotationProfile(.5,.25),1);assert.equal(rotationProfile(.25,.5),.5);
 assert.equal(rotationProfile(0,0),0);assert.equal(expression('rotation(x,y)')(.5,.25),1);
 assert.deepEqual(velocityAt(s,.5,.25),[.25,0]);
 assert.ok(Math.abs(u(.75,.5,Math.PI/2)-1)<1e-12);
 for(const [x,y] of [[.45,.76],[.5,.25],[.25,.5],[.1,.1]])assert.ok(Math.abs(u(x,y,2*Math.PI)-u(x,y,0))<1e-12);
 s.initial='0';assert.equal(exactSolution(s),null);
});
test('affine divergence-free rotation has skew interior convection and matches full AFC assembly',()=>{
 const s={...preset('solid-rotation'),n:8},a=assemble(s),g=buildAfc(s,a),{N,A}=a;
 for(let i=0;i<N;i++)for(let j=0;j<N;j++){
  assert.ok(Math.abs(A[i*N+j]+A[j*N+i])<1e-14);
  assert.ok(Math.abs(A[i*N+j]-(g.rows[a.mesh.interior[i]].get(a.mesh.interior[j])||0))<1e-14);
 }
 // A constant transported field has zero interior residual, including boundary couplings.
 for(const id of a.mesh.interior)assert.ok(Math.abs([...g.rows[id].values()].reduce((x,y)=>x+y,0))<1e-14);
});
test('full revolution FCT remains finite and within initial extrema',()=>{
 const s={...preset('solid-rotation'),n:12,dt:Math.PI/20},r=computeTrajectory(s),hi=Math.max(...r.frames[0].stable);
 for(const f of r.frames){assert.ok(f.stable.every(Number.isFinite));assert.ok(Math.min(...f.stable)>-1e-12);assert.ok(Math.max(...f.stable)<=hi+1e-12);}
 assert.ok(Math.abs(r.frames.at(-1).time-2*Math.PI)<1e-12);
 assert.ok(r.frames.at(-1).reference.every((v,i)=>Math.abs(v-r.frames[0].reference[i])<1e-12));
});
test('rotation SUPG has finite local stabilization and follows the same reference',()=>{
 const s={...preset('solid-rotation'),n:8,stabilization:'supg',duration:.1,dt:.05};
 const r=computeTrajectory(s);assert.ok(r.tauRange[1]>0);assert.ok(r.frames.at(-1).stable.every(Number.isFinite));assert.ok(r.frames.at(-1).reference);
});

test('rotation with diffusion adds epsilon K and keeps the same velocity and initial field',()=>{
 const c={...preset('solid-rotation'),n:8},cd={...preset('rotation-diffusion'),n:8};
 const a=assemble(c),b=assemble(cd);
 assert.equal(cd.mode,'parabolic');assert.equal(cd.initial,c.initial);
 for(let k=0;k<a.A.length;k++)assert.ok(Math.abs(b.A[k]-a.A[k]-cd.epsilon*a.K[k])<1e-14);
 assert.equal(exactSolution(cd),null);
 cd.epsilon=0;assert.ok(exactSolution(cd));
});
test('diffusive rotation stays bounded with FCT and omits a misleading inviscid reference',()=>{
 const s={...preset('rotation-diffusion'),n:12,epsilon:.001,dt:Math.PI/20};
 const r=computeTrajectory(s),hi=Math.max(...r.frames[0].stable);
 for(const f of r.frames){assert.equal(f.reference,null);assert.equal(f.metrics.stable.l2,null);assert.ok(f.stable.every(Number.isFinite));assert.ok(Math.min(...f.stable)>=-1e-12);assert.ok(Math.max(...f.stable)<=hi+1e-12);}
 assert.ok(Math.max(...r.frames.at(-1).stable)<hi);
});
