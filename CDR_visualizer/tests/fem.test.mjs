import test from 'node:test';
import assert from 'node:assert/strict';
import {makeMesh,preset,expression,assemble,compute,solveDense,interpolate,exactSolution} from '../dist/fem.js';
const near=(a,b,tol=1e-10)=>assert.ok(Math.abs(a-b)<tol,`${a} differs from ${b}`);

test('mesh has the requested slope +1 diagonals, counts, area, and interior ordering',()=>{
  const m=makeMesh(4);assert.equal(m.nodes.length,25);assert.equal(m.triangles.length,32);assert.equal(m.interior.length,9);
  assert.deepEqual(m.interior,[6,7,8,11,12,13,16,17,18]);assert.deepEqual(m.triangles.slice(0,2),[[0,1,6],[0,6,5]]);
  let total=0;for(const t of m.triangles){const[a,b,c]=t.map(i=>m.nodes[i]);const twice=(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);assert.ok(twice>0);total+=twice/2;}near(total,1);
});
test('arithmetic parser has conventional precedence and rejects JavaScript and invalid functions',()=>{
  near(expression('-2^2 + 2^-2 + 2^3^2')(0,0),508.25);
  near(expression('sin(pi*x) + exp(-y/eps)')(0.5,0,1),2);
  for(const s of ['window.location','x;alert(1)','sin x','2x','unknown(1)',''])assert.throws(()=>expression(s)(0,0));
  assert.throws(()=>expression('1/(x-x)')(0,0));
});
test('one interior node reproduces the hand-computed Poisson stiffness and constant load',()=>{
  const s=preset();s.n=2;s.source='1';const r=compute(s);
  near(r.A[0],4);near(r.K[0],4);near(r.F[0],.25);near(r.rawB[0],.0625);
});
test('Poisson stencil is symmetric with diagonal 4 and nearest-neighbor -1',()=>{
  const s=preset();s.n=4;const{A,N}=assemble(s);
  near(A[4*N+4],4);for(const j of[1,3,5,7])near(A[4*N+j],-1);
  for(let i=0;i<N;i++)for(let j=0;j<N;j++)near(A[i*N+j],A[j*N+i]);
});
test('affine patch test: convection, nonzero boundaries, and artificial-diffusion lifting',()=>{
  for(const pattern of ['forward','backward','alternating']){
    const s=preset();Object.assign(s,{n:6,pattern,mode:'transport',epsilon:.013,betaX:1.2,betaY:-.7,source:'0.3',diffusion:.07});
    s.boundaries={left:'2*x+3*y-1',right:'2*x+3*y-1',bottom:'2*x+3*y-1',top:'2*x+3*y-1'};
    const r=compute(s);r.mesh.nodes.forEach(([x,y],id)=>{near(r.raw[id],2*x+3*y-1);near(r.stable[id],2*x+3*y-1);});
    near(interpolate(r.mesh,r.raw,.37,.61),2*.37+3*.61-1);
    assert.ok(r.metrics.stable.backwardError<1e-13);
  }
});
test('smooth Poisson has second-order L2 convergence over n = 4, 8, 16, 32',()=>{
  let prev=Infinity;for(const n of [4,8,16,32]){const s=preset();s.n=n;const r=compute(s),err=r.metrics.raw.l2;assert.ok(err<prev/3.5);assert.ok(r.metrics.raw.backwardError<1e-13);prev=err;}assert.ok(prev<.0014);
});
test('sharp manufactured Poisson layer improves under refinement',()=>{
  const coarse=preset('sharp');coarse.n=8;const fine={...coarse,n:32};
  assert.ok(compute(fine).metrics.raw.l2<compute(coarse).metrics.raw.l2/5);
});
test('the documented transport example oscillates, then diffusion reduces its overshoot',()=>{
  const s=preset('layer');s.diffusion=.02625;const r=compute(s);
  near(r.pe,6.25);near(r.effectivePe,1);assert.ok(r.metrics.raw.max>1.8);assert.ok(r.metrics.stable.max<1);assert.ok(r.metrics.stable.min>-1e-10);
  assert.ok(r.metrics.stable.l2<r.metrics.raw.l2);
  near(exactSolution(s)(0,.5),0);near(exactSolution(s)(1,.5),0);near(exactSolution(s)(.5,.5),.5);
});
test('arbitrary matrix and RHS overrides are applied and invalidate the reference',()=>{
  const s=preset();s.n=2;s.source='1';s.matrixEdits={'0,0':8};s.rhsEdits={'0':2};const r=compute(s);
  near(r.rawB[0],.25);assert.equal(r.reference,null);
});
test('individual boundary values override edge data and conflicting corners are averaged',()=>{
  const s=preset();s.n=2;s.boundaries.left='2';s.nodeEdits={'0':7};const r=compute(s);
  near(r.raw[0],7);near(r.raw[3],2);near(r.raw[6],1);assert.equal(r.reference,null);
});
test('pivoting handles a zero leading pivot; singular and invalid inputs fail explicitly',()=>{
  const x=solveDense(new Float64Array([0,2,1,3]),new Float64Array([4,7]));near(x[0],1);near(x[1],2);
  assert.throws(()=>solveDense(new Float64Array([1,2,2,4]),new Float64Array([1,2])),/singular/);
  for(const n of [1,33,2.5,NaN])assert.throws(()=>makeMesh(n));
  const s=preset();s.epsilon=0;assert.throws(()=>compute(s));s.epsilon=1;s.diffusion=-1;assert.throws(()=>compute(s));
});
test('reference responds to physical diffusion and disappears for changed PDE data',()=>{
  const s=preset();s.epsilon=2;near(exactSolution(s)(.5,.5),.5);s.source='1';assert.equal(exactSolution(s),null);
});
