import test from 'node:test';
import assert from 'node:assert/strict';
import {families,experiments,examplesFor,defaultExample,equationFamily,timeSetting,compatibleMethods} from '../dist/experiments.js';
import {preset,compute,physicalDiffusion,velocity} from '../dist/fem.js';
import {computeTrajectory} from '../dist/transient.js';

test('all fourteen steady/transient and C/D/R combinations have matching examples and defaults',()=>{
 assert.equal(new Set(experiments.map(e=>e.id)).size,experiments.length);
 for(const time of ['steady','transient'])for(const family of Object.keys(families)){
  const list=examplesFor(time,family);assert.ok(list.length>0,`${time} ${family}`);
  assert.ok(list.some(e=>e.id===defaultExample(time,family)));
  for(const e of list){const s=preset(e.id);assert.equal(equationFamily(s),family);assert.equal(timeSetting(s),time);}
 }
});
test('every default equation combination assembles and solves with only its selected physical terms',()=>{
 for(const time of ['steady','transient'])for(const family of Object.keys(families)){
  const s=preset(defaultExample(time,family));if(time==='transient')Object.assign(s,{n:4,dt:.01,duration:.02});
  assert.equal(physicalDiffusion(s)>0,family.includes('D'));
  assert.equal(Math.hypot(...velocity(s))>0,family.includes('C'));
  assert.equal(s.reaction!==0,family.includes('R'));
  const r=time==='steady'?compute(s):computeTrajectory(s),u=r.frames?r.frames.at(-1).stable:r.stable;
  assert.ok(u.every(Number.isFinite),`${time} ${family}`);
 }
});
test('steady convection, reaction, and convection-reaction examples reproduce their compatible affine equilibrium',()=>{
 for(const id of ['steady-c','steady-r','steady-cr']){
  const r=compute(preset(id));assert.ok(r.metrics.raw.l2<1e-12,id);assert.ok(r.metrics.raw.backwardError<1e-12);
 }
 const dr=preset('steady-dr'),coarse=compute({...dr,n:4}),fine=compute({...dr,n:8});assert.ok(fine.metrics.raw.l2<coarse.metrics.raw.l2/3);
});
test('pure steady convection reports singular Galerkin meshes rather than silently adding diffusion',()=>{
 const s=preset('steady-c');s.n=16;assert.throws(()=>compute(s),/pure-convection Galerkin matrix is singular/);
 assert.equal(physicalDiffusion(s),0);
});
test('method choices identify BJK versus Zalesak and disable streamline stabilization without convection',()=>{
 for(const e of experiments){const s=preset(e.id),methods=compatibleMethods(s),supg=methods.find(m=>m.id==='supg'),afc=methods.find(m=>m.id==='afc');
  assert.equal(supg.enabled,e.family.includes('C'));
  assert.match(afc.label,e.time==='steady'?/BJK/:/Zalesak/);
  assert.ok(methods.find(m=>m.id===s.stabilization)?.enabled,e.id);
 }
});
