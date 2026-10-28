import {preset,isTransient} from './fem.js';

export const families={C:'Convection / advection',D:'Diffusion',R:'Reaction',CD:'Convection + diffusion',CR:'Convection + reaction',DR:'Diffusion + reaction',CDR:'Convection + diffusion + reaction'};
const convective=['transport','advection','parabolic','advection-reaction','adr','steady-adr','steady-convection','steady-convection-reaction'];
const diffusive=['poisson','transport','heat','parabolic','diffusion-reaction','adr','steady-adr','steady-diffusion-reaction'];
const reacting=['reaction','diffusion-reaction','advection-reaction','adr','steady-adr','steady-reaction','steady-convection-reaction','steady-diffusion-reaction'];
export function equationFamily(state){return (convective.includes(state.mode)?'C':'')+(diffusive.includes(state.mode)?'D':'')+(reacting.includes(state.mode)||(state.reaction??0)!==0?'R':'');}
export const timeSetting=state=>isTransient(state)?'transient':'steady';
const entries=[
 ['poisson','Smooth sine field','Refine the mesh and measure diffusion error.'],
 ['sharp','Sharp diffusion layer','A manufactured layer narrower than a coarse mesh.'],
 ['layer','Outflow layer · isotropic diffusion','Compare oscillations with layer smearing.'],
 ['supg-layer','Outflow layer · SUPG','Streamline stabilization on an unresolved boundary layer.'],
 ['afc-layer','Outflow layer · BJK AFC','Compare Galerkin with limited graph diffusion.'],
 ['afc-affine','Affine preservation · μ','Convection-dominated affine patch; check limiter activity and nodal error.'],
 ['steady-c','Compatible affine transport','Pure steady convection with compatible boundary traces.'],
 ['steady-r','Reaction equilibrium','The affine equilibrium cu = f with matching boundary values.'],
 ['steady-cr','Affine transport + reaction','Check the balance of convection, reaction, and forcing.'],
 ['steady-dr','Diffusion–reaction sine field','Manufactured smooth equilibrium with a known reference.'],
 ['steady-adr','CDR sine field','Manufactured convection–diffusion–reaction equilibrium.'],
 ['moving-front','Moving front · mass comparison','Watch a steep pulse move; lumping alone can leave oscillations.'],
 ['supg-front','Moving front · SUPG','Compare directional stabilization with the translating reference.'],
 ['rotation-diffusion','Rotating and diffusing · three shapes','Rotate the cylinder, cone, and hump while physical diffusion smooths their edges.'],
 ['solid-rotation','Solid-body rotation · three shapes','Rotate a slotted cylinder, cone, and smooth hump; compare shape preservation after one revolution (T = 2π).'],
 ['afc-front','Moving front · FCT','Zalesak-style limiting with lumped SSPRK2 substeps.'],
 ['heat-pulse','Diffusing pulse','Compare consistent and lumped mass near a narrow positive pulse.'],
 ['moving-pulse','Moving and diffusing pulse','Separate physical diffusion from convection and added diffusion.'],
 ['reaction-decay','Pure decay / growth','Change the reaction rate and time integrator.'],
 ['diffusion-reaction','Spreading and decay','Diffusion broadens a pulse while reaction changes amplitude.'],
 ['advection-reaction','Moving and decaying pulse','Follow convection and reaction without physical diffusion.'],
 ['adr-pulse','Moving, spreading, and decay','All three physical terms acting on an analytic Gaussian.'],
 ['forming-layer','Developing outflow layer','A manufactured transient CDR boundary layer.']
];
export const experiments=entries.map(([id,label,description])=>{const s=preset(id);return{id,label,description,time:timeSetting(s),family:equationFamily(s)};});
export function examplesFor(time,family){return experiments.filter(e=>e.time===time&&e.family===family);}
const defaults={steady:{C:'steady-c',D:'poisson',R:'steady-r',CD:'afc-layer',CR:'steady-cr',DR:'steady-dr',CDR:'steady-adr'},transient:{C:'afc-front',D:'heat-pulse',R:'reaction-decay',CD:'moving-pulse',CR:'advection-reaction',DR:'diffusion-reaction',CDR:'adr-pulse'}};
export function defaultExample(time,family){return defaults[time]?.[family]??'afc-layer';}
export function compatibleMethods(state){const steady=!isTransient(state),flow=equationFamily(state).includes('C');return [
 {id:'none',label:'Galerkin · no stabilization',enabled:true},
 {id:'isotropic',label:'Isotropic artificial diffusion · δK',enabled:true},
 {id:'supg',label:'SUPG · full streamline residual',enabled:flow},
 {id:'afc-low',label:'Graph diffusion · low order',enabled:true},
 {id:'afc',label:steady?'AFC · BJK-type limiter':'FCT · Zalesak-style limiter (SSPRK2)',enabled:true}
];}
