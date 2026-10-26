import test from 'node:test';
import assert from 'node:assert/strict';
import {preset,assemble,exactSolution,expression,solveDense,physicalDiffusion} from '../dist/fem.js';
import {computeTrajectory,factorDense,solveFactored} from '../dist/transient.js';
const near=(a,b,tol=1e-10)=>assert.ok(Math.abs(a-b)<tol,`${a} differs from ${b}`);

test('consistent and lumped P1 masses match the hand-computed single-DOF case',()=>{
  const s=preset('heat-pulse');s.n=2;const r=assemble(s);
  near(r.M[0],.125);near(r.lumped[0],.25);
  near(r.M[0]+r.boundaryTerms.reduce((v,t)=>v+t[4],0),r.lumped[0]);
});
test('reusable partial-pivot factors reproduce the independent direct solver on multiple RHS',()=>{
  const A=new Float64Array([0,2,3,4,5,6,7,8,10]),f=factorDense(A,3);
  for(const rhs of [new Float64Array([1,2,3]),new Float64Array([-4,8,1])]){
    const expected=solveDense(A,rhs),actual=solveFactored(f,rhs);actual.forEach((v,i)=>near(v,expected[i]));
  }
});
test('time-varying affine patch test includes consistent-mass boundary lifting, convection, and δK',()=>{
  for(const pattern of ['forward','backward','alternating'])for(const scheme of ['be','cn']){
    const s=preset('moving-pulse');Object.assign(s,{n:4,pattern,scheme,epsilon:.013,betaX:.3,betaY:.2,source:'1.7',initial:'x+2*y',dt:.03,duration:.1,diffusion:.08});
    s.boundaries=Object.fromEntries(['left','right','bottom','top'].map(k=>[k,'t+x+2*y']));
    const r=computeTrajectory(s);
    near(r.frames.at(-1).time,.1);near(r.frames.at(-1).dt,.01);
    for(const f of r.frames)r.mesh.nodes.forEach(([x,y],i)=>{near(f.raw[i],f.time+x+2*y);near(f.stable[i],f.time+x+2*y);});
  }
});
test('Crank–Nicolson integrates a quadratic-in-time, spatially constant manufactured solution',()=>{
  for(const massType of ['consistent','lumped']){
    const s=preset('heat-pulse');Object.assign(s,{n:4,scheme:'cn',massType,source:'2*t',initial:'0',dt:.02,duration:.1});
    s.boundaries=Object.fromEntries(['left','right','bottom','top'].map(k=>[k,'t^2']));
    const r=computeTrajectory(s);for(const f of r.frames)for(const v of [...f.raw,...f.stable])near(v,f.time**2);
  }
});
test('lumping before boundary elimination preserves constant states for all three transient modes',()=>{
  for(const mode of ['heat','advection','parabolic']){
    const s=preset('moving-pulse');Object.assign(s,{mode,n:4,initial:'1',source:'0',dt:.01,duration:.03});
    s.boundaries={left:'1',right:'1',bottom:'1',top:'1'};
    const r=computeTrajectory(s);for(const f of r.frames)for(const v of [...f.raw,...f.stable])near(v,1);
  }
});
test('heat-pulse example demonstrates negativity in consistent mass and positivity with lumped mass + backward Euler',()=>{
  const s=preset('heat-pulse');s.duration=.02;const r=computeTrajectory(s);
  assert.ok(Math.min(...r.frames.map(f=>f.metrics.raw.min)) < -.01);
  assert.ok(Math.min(...r.frames.map(f=>f.metrics.stable.min)) >= -1e-12);
  for(const f of r.frames.slice(1)){assert.ok(f.metrics.raw.backwardError<1e-12);assert.ok(f.metrics.stable.backwardError<1e-12);}
});
test('mass selection affects only the comparison; consistent mass and δ = 0 reproduce baseline exactly',()=>{
  const s=preset('heat-pulse');Object.assign(s,{n:6,duration:.01});const lumped=computeTrajectory(s);s.massType='consistent';const consistent=computeTrajectory(s);
  for(let k=0;k<lumped.frames.length;k++){
    assert.deepEqual(lumped.frames[k].raw,consistent.frames[k].raw);
    assert.deepEqual(consistent.frames[k].raw,consistent.frames[k].stable);
  }
});
test('pure convection reference translates with velocity and ε = 0 needs no steady spatial inversion',()=>{
  const s=preset('moving-front');Object.assign(s,{n:8,betaX:.7,betaY:.1,duration:.05,dt:.01});
  const exact=exactSolution(s),initial=expression(s.initial);near(exact(.35,.51,.1),initial(.28,.5));near(physicalDiffusion(s),0);
  const r=computeTrajectory(s);assert.equal(r.pe,Infinity);assert.ok(r.frames.at(-1).metrics.raw.backwardError<1e-12);
  // Row-sum lumping alone is not a positivity guarantee for central Galerkin transport.
  assert.ok(r.frames.at(-1).metrics.stable.min<0);
});
test('Gaussian exact reference broadens with physical diffusion, and invalid initial data hide it',()=>{
  const s=preset('moving-pulse');s.epsilon=.02;s.betaX=.6;const fn=exactSolution(s),t=.2,w=.0049+4*s.epsilon*t;
  near(fn(.25+s.betaX*t,.4+s.betaY*t,t),.0049/w);
  s.initial='1';assert.equal(exactSolution(s),null);
});
test('time grid, parser variables, and resource limits are explicit',()=>{
  near(expression('t+vx*x+vy*y+eps')(.2,.3,.01,.4,.5,.6),.69);
  const s=preset('heat-pulse');s.dt=0;assert.throws(()=>computeTrajectory(s),/positive/);
  s.dt=1e-6;assert.throws(()=>computeTrajectory(s),/240/);
  s.dt=.1;s.duration=.03;const r=computeTrajectory(s);assert.equal(r.steps,1);near(r.frames[1].time,.03);near(r.frames[1].dt,.03);
  s.n=2;s.duration=1e-15;assert.equal(computeTrajectory(s).frames.length,2);
});
