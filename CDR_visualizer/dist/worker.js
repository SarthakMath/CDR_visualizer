import {compute,isTransient} from './fem.js?v=257aed793d6ea050';
import {computeTrajectory} from './transient.js?v=257aed793d6ea050';
self.onmessage=({data:{id,state}})=>{
  try { self.postMessage({id,result:isTransient(state)?computeTrajectory(state,(step,total)=>self.postMessage({id,progress:step/total})):compute(state)}); }
  catch(error){self.postMessage({id,error:error.message});}
};
