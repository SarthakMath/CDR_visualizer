import test from 'node:test';
import assert from 'node:assert/strict';
import {preset,assemble,compute,exactSolution,expression,boundaryEvaluator,physicalDiffusion,velocity} from '../dist/fem.js';
import {computeTrajectory} from '../dist/transient.js';
const near=(a,b,tol=1e-10)=>assert.ok(Math.abs(a-b)<tol,`${a} differs from ${b}`);
test('reaction adds cM and the corresponding nonzero boundary lift',()=>{
 const s=preset();s.n=4;s.boundaries={left:'1',right:'1',bottom:'1',top:'1'};
 const base=assemble(s);s.reaction=3;const r=assemble(s);
 for(let i=0;i<r.A.length;i++)near(r.A[i]-base.A[i],3*r.M[i]);
 for(let i=0;i<r.N;i++)near(r.F[i]-base.F[i],-3*r.boundaryTerms.filter(t=>t[0]===i).reduce((sum,t)=>sum+t[4],0));
});
test('steady ADR reproduces an affine manufactured solution on all patterns',()=>{
 for(const pattern of ['forward','backward','alternating']){
 const s=preset('steady-adr');Object.assign(s,{n:4,pattern,reaction:3,source:'vx+2*vy+c*(1+x+2*y)',diffusion:.1});
 s.boundaries=Object.fromEntries(['left','right','bottom','top'].map(k=>[k,'1+x+2*y']));
 const r=compute(s);r.mesh.nodes.forEach(([x,y],i)=>{near(r.raw[i],1+x+2*y);near(r.stable[i],1+x+2*y);});
 }
});
test('transient ADR patch verifies reaction, changing boundary, both integrators and mass choices',()=>{
 for(const scheme of ['be','cn'])for(const massType of ['consistent','lumped']){
 const s=preset('adr-pulse');Object.assign(s,{n:4,scheme,massType,reaction:3,initial:'1+x+2*y',source:'1+vx+2*vy+c*(1+t+x+2*y)',dt:.03,duration:.1,diffusion:.05});
 s.boundaries=Object.fromEntries(['left','right','bottom','top'].map(k=>[k,'1+t+x+2*y']));
 const r=computeTrajectory(s);for(const f of r.frames)r.mesh.nodes.forEach(([x,y],i)=>{near(f.raw[i],1+f.time+x+2*y);near(f.stable[i],1+f.time+x+2*y);});
 }
});
test('pure reaction has the analytic BE/CN amplification for decay, zero and growth',()=>{
 for(const c of [2,0,-1])for(const scheme of ['be','cn']){
 const s=preset('reaction-decay');Object.assign(s,{n:4,reaction:c,scheme,dt:.02,duration:.1});s.boundaries={left:'0',right:'0',bottom:'0',top:'0'};
 const r=computeTrajectory(s),theta=scheme==='be'?1:.5,q=(1-(1-theta)*s.dt*c)/(1+theta*s.dt*c);
 for(const [k,f] of r.frames.entries())for(const id of r.mesh.interior)near(f.raw[id],r.frames[0].raw[id]*q**k);
 }
});
test('reaction preset traces and references satisfy full PDE, including edited coefficients',()=>{
 for(const name of ['reaction-decay','diffusion-reaction','advection-reaction','adr-pulse']){
 const s=preset(name);s.reaction=-.7;s.epsilon=.01;s.betaX=.3;
 const u=exactSolution(s),h=1e-4,x=.38,y=.45,t=.2;
 const ut=(u(x,y,t+h)-u(x,y,t-h))/(2*h),ux=(u(x+h,y,t)-u(x-h,y,t))/(2*h),uy=(u(x,y+h,t)-u(x,y-h,t))/(2*h);
 const lap=(u(x+h,y,t)+u(x-h,y,t)+u(x,y+h,t)+u(x,y-h,t)-4*u(x,y,t))/h**2;
 const [vx,vy]=velocity(s);near(ut-physicalDiffusion(s)*lap+vx*ux+vy*uy+s.reaction*u(x,y,t),0,2e-5);
 const r=assemble({...s,n:4}),g=boundaryEvaluator(s,r.mesh)(t);for(const id of r.mesh.boundary){const[x,y]=r.mesh.nodes[id];near(g[id],u(x,y,t));}
 }
});
test('old references are hidden after adding reaction and invalid coefficients fail',()=>{
 const s=preset();s.reaction=2;assert.equal(exactSolution(s),null);
 s.reaction=NaN;assert.throws(()=>assemble(s),/Reaction/);
 near(expression('exp(-c*t)')(.2,.3,0,.5,0,0,2),Math.exp(-1));
});
test('steady ADR sine reference converges under mesh refinement',()=>{
 const s=preset('steady-adr');s.n=8;const coarse=compute(s).metrics.raw.l2;s.n=16;const fine=compute(s).metrics.raw.l2;
 assert.ok(fine<coarse/2);
});
