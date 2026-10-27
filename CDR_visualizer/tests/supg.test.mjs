import test from 'node:test';
import assert from 'node:assert/strict';
import {preset,assemble,compute,exactSolution,supgTau} from '../dist/fem.js';
import {computeTrajectory} from '../dist/transient.js';
const near=(a,b,tol=1e-10)=>assert.ok(Math.abs(a-b)<tol,`${a} differs from ${b}`);
test('SUPG vanishes without flow or at zero strength; disabled δ is ignored',()=>{
 const s=preset();s.stabilization='supg';s.supgScale=1;const r=compute(s);
 assert.ok(r.S.every(v=>v===0)&&r.H.every(v=>v===0));assert.deepEqual(r.raw,r.stable);
 const t=preset('supg-front');Object.assign(t,{n:4,duration:.02,supgScale:0,diffusion:100});const run=computeTrajectory(t);
 for(const f of run.frames)f.raw.forEach((v,i)=>near(v,f.stable[i]));
 t.stabilization='none';const off=computeTrajectory(t);for(const f of off.frames)assert.deepEqual(f.raw,f.stable);
});
test('tau has stable diffusion and zero-diffusion limits and linear strength scaling',()=>{
 const s=preset('supg-layer'),grad=[[-4,0],[4,-4],[0,4]];
 s.epsilon=0;near(supgTau(s,grad),.125);s.supgScale=2;near(supgTau(s,grad),.25);
 s.supgScale=1;s.epsilon=1e8;near(supgTau(s,grad),.25**2/(12*s.epsilon),1e-18);
 s.supgScale=-1;assert.throws(()=>assemble(s),/SUPG strength/);
});
test('steady residual SUPG reproduces an affine ADR patch with nonzero boundary and source',()=>{
 for(const pattern of ['forward','backward','alternating']){
 const s=preset('steady-adr');Object.assign(s,{n:4,pattern,stabilization:'supg',betaX:.7,betaY:.3,reaction:2,source:'vx+2*vy+c*(1+x+2*y)'});
 s.boundaries=Object.fromEntries(['left','right','bottom','top'].map(k=>[k,'1+x+2*y']));
 const r=compute(s);r.mesh.nodes.forEach(([x,y],i)=>near(r.stable[i],1+x+2*y));
 assert.ok(r.S.some(v=>v!==0));assert.ok(r.H.some(v=>v!==0));
 }
});
test('transient SUPG includes H, reaction residual, time-varying load, and boundary H lift',()=>{
 for(const pattern of ['forward','backward','alternating'])for(const scheme of ['be','cn'])for(const massType of ['consistent','lumped']){
 const s=preset('adr-pulse');Object.assign(s,{n:4,pattern,scheme,massType,stabilization:'supg',betaX:.7,betaY:.3,reaction:2,initial:'1+x+2*y',source:'1+vx+2*vy+c*(1+t+x+2*y)',dt:.03,duration:.1});
 s.boundaries=Object.fromEntries(['left','right','bottom','top'].map(k=>[k,'1+t+x+2*y']));
 const r=computeTrajectory(s);for(const f of r.frames)r.mesh.nodes.forEach(([x,y],i)=>near(f.stable[i],1+f.time+x+2*y));
 }
});
test('SUPG adds no streamline stiffness to a crosswind-only profile',()=>{
 const s=preset('supg-layer');s.n=4;const r=assemble(s);
 const y=Float64Array.from(r.mesh.interior,id=>r.mesh.nodes[id][1]**2);
 for(let i=0;i<r.N;i++){let sum=0;for(let j=0;j<r.N;j++)sum+=r.S[i*r.N+j]*y[j];for(const[row,id,,,,v]of r.boundaryTerms)if(row===i)sum+=v*r.mesh.nodes[id][1]**2;near(sum,0);}
});
test('steady and developing outflow layers reduce overshoot and error with SUPG',()=>{
 const s=preset('supg-layer'),r=compute(s);
 assert.ok(r.metrics.raw.max>1.8);assert.ok(r.metrics.stable.max<1);assert.ok(r.metrics.stable.l2<r.metrics.raw.l2);
 const t=preset('forming-layer'),run=computeTrajectory(t),last=run.frames.at(-1),ref=exactSolution(t);
 near(ref(.5,.5,1),(1-Math.exp(-1))*(.5-(Math.exp(-.5/t.epsilon)-Math.exp(-1/t.epsilon))/(1-Math.exp(-1/t.epsilon))));
 assert.ok(last.metrics.raw.max>1);assert.ok(last.metrics.stable.max<.64);assert.ok(last.metrics.stable.l2<last.metrics.raw.l2);
 t.betaX=.5;assert.equal(exactSolution(t),null);
});
test('moving front SUPG reduces extrema and error but retains undershoot',()=>{
 const s=preset('supg-front'),r=computeTrajectory(s),last=r.frames.at(-1);
 assert.ok(last.metrics.stable.min>last.metrics.raw.min);assert.ok(last.metrics.stable.max<last.metrics.raw.max);
 assert.ok(last.metrics.stable.min<-.01);assert.ok(last.metrics.stable.l2<last.metrics.raw.l2);
 const before=r.frames[0];s.stabilization='none';const off=computeTrajectory(s);
 assert.deepEqual(before.raw,off.frames[0].raw);assert.deepEqual(last.raw,off.frames.at(-1).raw);
});
