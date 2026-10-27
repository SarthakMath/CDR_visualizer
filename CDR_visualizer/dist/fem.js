import {isAfc,solveAfcSteady} from './afc.js';
// Linear (P1) triangular Galerkin finite elements. No numerical dependencies.
export const layerTrace = 'x-(exp((x-1)/eps)-exp(-1/eps))/(1-exp(-1/eps))';
export const sharpQ = '(x-(exp((x-1)/0.025)-exp(-40))/(1-exp(-40)))';
export const isTransient = state => ['heat','advection','parabolic','reaction','diffusion-reaction','advection-reaction','adr'].includes(state.mode);
export const velocity = state => ['transport','advection','parabolic','advection-reaction','adr','steady-adr'].includes(state.mode)?[state.betaX,state.betaY]:[0,0];
export const physicalDiffusion = state => ['advection','advection-reaction','reaction'].includes(state.mode)?0:state.epsilon;
export const addedDiffusion = state => ['none','supg','afc','afc-low'].includes(state.stabilization)?0:state.diffusion;
export const frontInitial = '0.5*(tanh((x-0.15)/0.015)-tanh((x-0.35)/0.015))*exp(-((y-0.5)/0.18)^2)';
const movingFront = frontInitial.replace(/\bx\b/g,'(x-vx*t)').replace(/\by\b/g,'(y-vy*t)');
const gaussian = (cx,cy,s2) => `${s2}/(${s2}+4*eps*t)*exp(-((x-${cx}-vx*t)^2+(y-${cy}-vy*t)^2)/(${s2}+4*eps*t))`;
export function preset(name = 'poisson') {
  const s = { preset:name, mode:'poisson', n:16, pattern:'forward', epsilon:1,
    betaX:0, betaY:0, reaction:0, source:'2*pi^2*sin(pi*x)*sin(pi*y)', amplitude:1,
    boundaries:{left:'0',right:'0',bottom:'0',top:'0'}, diffusion:0,stabilization:'isotropic',supgScale:1,afcMu:1,afcMaxIterations:500,afcTolerance:1e-8,afcRelaxation:.5,
    matrixEdits:{}, rhsEdits:{}, nodeEdits:{},initial:'sin(pi*x)*sin(pi*y)',dt:0.005,duration:0.8,scheme:'be',massType:'consistent' };
  if (name === 'layer') Object.assign(s, {mode:'transport',epsilon:0.005,betaX:1,source:'1',
    boundaries:{left:'0',right:'0',bottom:layerTrace,top:layerTrace}});
  if (name === 'sharp') s.source = `(exp((x-1)/0.025)/(0.025^2*(1-exp(-40)))+pi^2*${sharpQ})*sin(pi*y)`;
  if(name==='moving-front')Object.assign(s,{mode:'advection',n:24,epsilon:0,betaX:.55,betaY:0,source:'0',initial:frontInitial,scheme:'cn',massType:'lumped',boundaries:Object.fromEntries(['left','right','bottom','top'].map(k=>[k,movingFront]))});
  if(name==='heat-pulse')Object.assign(s,{mode:'heat',epsilon:.01,source:'0',initial:'exp(-((x-0.5)^2+(y-0.5)^2)/0.000625)',dt:.001,duration:.1,massType:'lumped',boundaries:Object.fromEntries(['left','right','bottom','top'].map(k=>[k,gaussian(.5,.5,.000625)]))});
  if(name==='moving-pulse')Object.assign(s,{mode:'parabolic',n:24,epsilon:.002,betaX:.5,betaY:.1,source:'0',initial:'exp(-((x-0.25)^2+(y-0.4)^2)/0.0049)',massType:'lumped',boundaries:Object.fromEntries(['left','right','bottom','top'].map(k=>[k,gaussian(.25,.4,.0049)]))});
  if(['reaction-decay','diffusion-reaction','advection-reaction','adr-pulse'].includes(name)){
    const modes={'reaction-decay':'reaction','diffusion-reaction':'diffusion-reaction','advection-reaction':'advection-reaction','adr-pulse':'adr'};
    Object.assign(s,{mode:modes[name],n:24,epsilon:['reaction-decay','advection-reaction'].includes(name)?0:.002,betaX:['reaction-decay','diffusion-reaction'].includes(name)?0:.5,betaY:['reaction-decay','diffusion-reaction'].includes(name)?0:.1,reaction:2,source:'0',initial:'exp(-((x-0.25)^2+(y-0.4)^2)/0.0049)',scheme:'cn',massType:'consistent',boundaries:Object.fromEntries(['left','right','bottom','top'].map(k=>[k,`exp(-c*t)*(${gaussian(.25,.4,.0049)})`]))});
  }
  if(name==='steady-adr')Object.assign(s,{mode:'steady-adr',epsilon:.02,betaX:1,reaction:2,source:'(2*eps*pi^2+c)*sin(pi*x)*sin(pi*y)+pi*vx*cos(pi*x)*sin(pi*y)+pi*vy*sin(pi*x)*cos(pi*y)'});
  if(name==='supg-layer')Object.assign(s,preset('layer'),{preset:name,stabilization:'supg'});
  if(name==='supg-front')Object.assign(s,preset('moving-front'),{preset:name,stabilization:'supg',massType:'consistent',boundaries:Object.fromEntries(['left','right','bottom','top'].map(k=>[k,`exp(-c*t)*(${movingFront})`]))});
  if(name==='forming-layer'){
    const q=`(${layerTrace})`,trace=`(1-exp(-t))*${q}`;
    Object.assign(s,{mode:'adr',epsilon:.005,betaX:1,reaction:.5,source:`exp(-t)*${q}+(1-exp(-t))*(1+c*${q})`,initial:'0',duration:1,dt:.01,scheme:'cn',stabilization:'supg',boundaries:{left:'0',right:'0',bottom:trace,top:trace}});
  }
  if(name==='afc-front')Object.assign(s,preset('supg-front'),{preset:name,stabilization:'afc',massType:'lumped'});
  if(name==='afc-layer')Object.assign(s,preset('supg-layer'),{preset:name,stabilization:'afc'});
  if(name==='afc-affine')Object.assign(s,{mode:'transport',n:8,pattern:'forward',epsilon:.005,betaX:1,betaY:.3,source:'vx+vy+c*(1+x+y)',boundaries:Object.fromEntries(['left','right','bottom','top'].map(k=>[k,'1+x+y'])),stabilization:'afc',afcMu:1,afcTolerance:1e-11,afcMaxIterations:1500});
  return s;
}

// Small arithmetic parser. Expressions never execute JavaScript.
export function expression(source) {
  if (typeof source !== 'string' || source.length > 800) throw Error('Use an expression shorter than 800 characters.');
  const tokens=[]; const re=/\s*(?:(\d+(?:\.\d*)?(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?)|([A-Za-z_]\w*)|([+\-*/^(),]))/gy;
  let offset=0;
  while(offset<source.length) {
    if (!source.slice(offset).trim()) break;
    re.lastIndex=offset; const m=re.exec(source);
    if(!m) throw Error(`Unexpected character near “${source.slice(offset,offset+12)}”.`);
    tokens.push(m[1] ? {number:Number(m[1])} : m[2] || m[3]); offset=re.lastIndex;
  }
  let p=0;
  const funcs={sin:Math.sin,cos:Math.cos,tan:Math.tan,exp:Math.exp,log:Math.log,sqrt:Math.sqrt,abs:Math.abs,tanh:Math.tanh};
  function atom(){
    const t=tokens[p++];
    if(t && typeof t==='object') return ()=>t.number;
    if(t==='(') {const a=sum(); if(tokens[p++]!==')') throw Error('Missing closing parenthesis.'); return a;}
    if(['x','y','eps','t','vx','vy','c'].includes(t)) return (...v)=>v[['x','y','eps','t','vx','vy','c'].indexOf(t)];
    if(t==='pi') return ()=>Math.PI;
    if(t==='e') return ()=>Math.E;
    if(Object.hasOwn(funcs,t)) {if(tokens[p++]!=='(')throw Error(`${t} needs parentheses.`);const a=sum();if(tokens[p++]!==')')throw Error('Missing closing parenthesis.');return (...v)=>funcs[t](a(...v));}
    throw Error(`Unknown term “${t ?? ''}”. Use x, y, t, eps, vx, vy, c, pi and standard math functions.`);
  }
  function power(){let a=atom(); if(tokens[p]==='^'){p++;const b=unary();const prev=a;a=(...v)=>prev(...v)**b(...v);}return a;}
  function unary(){if(tokens[p]==='-'){p++;const a=unary();return (...v)=>-a(...v);}if(tokens[p]==='+'){p++;return unary();}return power();}
  function product(){let a=unary();while(tokens[p]==='*'||tokens[p]==='/'){const op=tokens[p++],b=unary(),prev=a;a=(...v)=>op==='*'?prev(...v)*b(...v):prev(...v)/b(...v);}return a;}
  function sum(){let a=product();while(tokens[p]==='+'||tokens[p]==='-'){const op=tokens[p++],b=product(),prev=a;a=(...v)=>op==='+'?prev(...v)+b(...v):prev(...v)-b(...v);}return a;}
  const fn=sum();if(p!==tokens.length)throw Error('Unexpected term. Use * for multiplication.');
  return (x,y,eps=1,t=0,vx=0,vy=0,c=0)=>{const v=fn(x,y,eps,t,vx,vy,c);if(!Number.isFinite(v))throw Error('An expression produced a non-finite value. Check denominators and function domains.');return v;};
}

export function makeMesh(n,pattern='forward') {
  if(!Number.isInteger(n)||n<2||n>32)throw Error('Mesh subdivisions must be an integer from 2 to 32.');
  if(!['forward','backward','alternating','shifted'].includes(pattern))throw Error('Unknown triangulation pattern.');
  const nodes=[],triangles=[],interior=[],boundary=[],dof=new Int32Array((n+1)**2).fill(-1);
  for(let j=0;j<=n;j++)for(let i=0;i<=n;i++){
    const id=j*(n+1)+i;nodes.push([i/n,j/n]);
    if(i===0||j===0||i===n||j===n)boundary.push(id);else{dof[id]=interior.length;interior.push(id);}
  }
  for(let j=0;j<n;j++)for(let i=0;i<n;i++){
    const a=j*(n+1)+i,b=a+1,c=a+n+1,d=c+1;
    if(['forward','shifted'].includes(pattern)||(pattern==='alternating'&&(i+j)%2===0))triangles.push([a,b,d],[a,d,c]);
    else triangles.push([a,b,c],[b,d,c]);
  }
  if(pattern==='shifted'){const k=Math.floor(n/2),id=k*(n+1)+k;nodes[id][0]+=.35/n;nodes[id][1]+=.2/n;}
  return {n,pattern,nodes,triangles,interior,boundary,dof};
}

// Degree-five, seven-point Dunavant rule; weights sum to one.
export const quadrature=[
  [[1/3,1/3,1/3],0.225],
  ...[[0.05971587178977,0.470142064105115,0.470142064105115],[0.470142064105115,0.05971587178977,0.470142064105115],[0.470142064105115,0.470142064105115,0.05971587178977]].map(q=>[q,0.132394152788506]),
  ...[[0.797426985353087,0.101286507323456,0.101286507323456],[0.101286507323456,0.797426985353087,0.101286507323456],[0.101286507323456,0.101286507323456,0.797426985353087]].map(q=>[q,0.125939180544827])
];

export function boundaryEvaluator(state,mesh){
  const epsilon=physicalDiffusion(state),beta=velocity(state);
  const boundaries=Object.fromEntries(Object.entries(state.boundaries).map(([k,v])=>[k,expression(v)]));
  return time=>{
  const g=new Float64Array(mesh.nodes.length);
  for(const id of mesh.boundary){
    const [x,y]=mesh.nodes[id],values=[];
    if(x===0)values.push(boundaries.left(x,y,epsilon,time,...beta,state.reaction??0));
    if(x===1)values.push(boundaries.right(x,y,epsilon,time,...beta,state.reaction??0));
    if(y===0)values.push(boundaries.bottom(x,y,epsilon,time,...beta,state.reaction??0));
    if(y===1)values.push(boundaries.top(x,y,epsilon,time,...beta,state.reaction??0));
    g[id]=values.reduce((a,b)=>a+b,0)/values.length;
  }
  for(const [id,value] of Object.entries(state.nodeEdits||{})) {
    if(!Number.isInteger(+id)||+id<0||+id>=g.length||mesh.dof[+id]>=0)throw Error('Boundary node edit refers to an invalid node.');
    if(!Number.isFinite(value))throw Error('Boundary values must be finite.');g[+id]=value;
  }return g;
  };
}

function triangleArea(ps){return Math.abs((ps[1][0]-ps[0][0])*(ps[2][1]-ps[0][1])-(ps[2][0]-ps[0][0])*(ps[1][1]-ps[0][1]))/2;}

export function sourceLoad(mesh,source,state,time=0){
  const F=new Float64Array(mesh.interior.length),epsilon=physicalDiffusion(state),beta=velocity(state);
  for(const tri of mesh.triangles){
    const ps=tri.map(id=>mesh.nodes[id]),area=triangleArea(ps);
    for(const[q,w]of quadrature){
      const x=q.reduce((s,v,i)=>s+v*ps[i][0],0),y=q.reduce((s,v,i)=>s+v*ps[i][1],0),v=source(x,y,epsilon,time,...beta,state.reaction??0)*state.amplitude*area*w;
      for(let i=0;i<3;i++){const row=mesh.dof[tri[i]];if(row>=0)F[row]+=v*q[i];}
    }
  }return F;
}

// Semidiscrete residual SUPG. The element Laplacian vanishes for affine P1
// basis functions and constant diffusivity. H multiplies the time derivative.
export function supgTau(state,grad){
  if(state.stabilization!=='supg')return 0;
  const beta=velocity(state),speed=Math.hypot(...beta);if(speed===0)return 0;
  const h=2*speed/grad.reduce((sum,g)=>sum+Math.abs(beta[0]*g[0]+beta[1]*g[1]),0);
  const eps=physicalDiffusion(state),pe=eps>0?speed*h/(2*eps):Infinity;
  const factor=pe<1e-3?pe/3-pe**3/45:pe===Infinity?1:1/Math.tanh(pe)-1/pe;
  const tau=h/(2*speed)*factor;
  return (state.supgScale??1)*tau/Math.hypot(1,(state.reaction??0)*tau);
}
export function supgLoad(mesh,source,state,time,elements){
  const load=new Float64Array(mesh.interior.length),eps=physicalDiffusion(state),beta=velocity(state);
  for(const {tri,area,tau,stream} of elements){
    if(tau===0)continue;
    let integral=0;
    for(const[q,w]of quadrature){const x=q.reduce((a,v,i)=>a+v*mesh.nodes[tri[i]][0],0),y=q.reduce((a,v,i)=>a+v*mesh.nodes[tri[i]][1],0);integral+=area*w*state.amplitude*source(x,y,eps,time,...beta,state.reaction??0);}
    tri.forEach((id,i)=>{const row=mesh.dof[id];if(row>=0)load[row]+=tau*stream[i]*integral;});
  }return load;
}

export function assemble(state) {
  const n=state.n,epsilon=physicalDiffusion(state);
  if(!Number.isFinite(epsilon)||epsilon<0||(!isTransient(state)&&epsilon===0))throw Error('Physical diffusion ε must be positive for steady problems and nonnegative for time-dependent problems.');
  for(const key of ['betaX','betaY','amplitude','diffusion'])if(!Number.isFinite(state[key]))throw Error(`${key} must be finite.`);
  if(!Number.isFinite(state.reaction??0))throw Error('Reaction coefficient c must be finite.');
  if(!['none','isotropic','supg','afc','afc-low'].includes(state.stabilization??'isotropic'))throw Error('Unknown stabilization method.');
  if(!Number.isFinite(state.supgScale??1)||(state.supgScale??1)<0)throw Error('SUPG strength must be finite and nonnegative.');
  if(state.diffusion<0)throw Error('Artificial diffusion must be nonnegative.');
  const mesh=makeMesh(n,state.pattern),N=mesh.interior.length;
  const A=new Float64Array(N*N),K=new Float64Array(N*N),F=new Float64Array(N),lift=new Float64Array(N);
  const S=new Float64Array(N*N),H=new Float64Array(N*N),supgElements=[];
  const M=new Float64Array(N*N),lumped=new Float64Array(N),boundaryTerms=[];
  const g=boundaryEvaluator(state,mesh)(0),source=expression(state.source),load=sourceLoad(mesh,source,state);
  F.set(load);
  const beta=velocity(state);
  for(const tri of mesh.triangles){
    const ps=tri.map(id=>mesh.nodes[id]);
    const twice=(ps[1][0]-ps[0][0])*(ps[2][1]-ps[0][1])-(ps[2][0]-ps[0][0])*(ps[1][1]-ps[0][1]);
    const area=Math.abs(twice)/2;
    const grad=ps.map((_,i)=>[(ps[(i+1)%3][1]-ps[(i+2)%3][1])/twice,(ps[(i+2)%3][0]-ps[(i+1)%3][0])/twice]);
    const tau=supgTau(state,grad),stream=grad.map(g=>beta[0]*g[0]+beta[1]*g[1]);
    supgElements.push({tri,area,tau,stream});
    for(let i=0;i<3;i++){
      const row=mesh.dof[tri[i]];if(row<0)continue;lumped[row]+=area/3;
      for(let j=0;j<3;j++){
        const col=mesh.dof[tri[j]],k=area*(grad[i][0]*grad[j][0]+grad[i][1]*grad[j][1]);
        const m=area*(i===j?1/6:1/12);
        const a=epsilon*k+area/3*(beta[0]*grad[j][0]+beta[1]*grad[j][1])+(state.reaction??0)*m;
        const h=tau*stream[i]*area/3,s=tau*stream[i]*area*stream[j]+(state.reaction??0)*h;
        if(col<0){F[row]-=a*g[tri[j]];lift[row]-=k*g[tri[j]];boundaryTerms.push([row,tri[j],a,k,m,s,h]);}
        else{A[row*N+col]+=a;K[row*N+col]+=k;M[row*N+col]+=m;S[row*N+col]+=s;H[row*N+col]+=h;}
      }
    }
  }
  for(const [key,value] of Object.entries(state.matrixEdits||{})){
    const [i,j]=key.split(',').map(Number);
    if(!Number.isInteger(i)||!Number.isInteger(j)||i<0||j<0||i>=N||j>=N||!Number.isFinite(value))throw Error('Invalid matrix edit.');
    A[i*N+j]=value;
  }
  for(const [key,value] of Object.entries(state.rhsEdits||{})){
    const i=+key;if(!Number.isInteger(i)||i<0||i>=N||!Number.isFinite(value))throw Error('Invalid right-hand-side edit.');F[i]=value;
  }
  const supgF=supgLoad(mesh,source,state,0,supgElements);
  return {mesh,N,A,K,M,S,H,supgElements,supgF,tauRange:[Math.min(...supgElements.map(e=>e.tau)),Math.max(...supgElements.map(e=>e.tau))],lumped,boundaryTerms,load,F,lift,g};
}

export function solveDense(A,rhs) {
  const n=rhs.length,a=A.slice(),b=rhs.slice();
  let norm=0;for(const v of a)norm=Math.max(norm,Math.abs(v));
  if(!Number.isFinite(norm)||norm===0)throw Error('The edited matrix is singular. Restore entries or reset the experiment.');
  const threshold=norm*Number.EPSILON*n*8;
  // Partial pivoting also supports nonsymmetric and manually edited systems.
  for(let k=0;k<n;k++){
    let pivot=k;for(let i=k+1;i<n;i++)if(Math.abs(a[i*n+k])>Math.abs(a[pivot*n+k]))pivot=i;
    if(Math.abs(a[pivot*n+k])<=threshold)throw Error('The edited matrix is singular or numerically ill-conditioned. Restore entries or reset the experiment.');
    if(pivot!==k){for(let j=k;j<n;j++){const t=a[k*n+j];a[k*n+j]=a[pivot*n+j];a[pivot*n+j]=t;}const t=b[k];b[k]=b[pivot];b[pivot]=t;}
    for(let i=k+1;i<n;i++){
      const factor=a[i*n+k]/a[k*n+k];if(factor===0)continue;a[i*n+k]=0;
      for(let j=k+1;j<n;j++)a[i*n+j]-=factor*a[k*n+j];b[i]-=factor*b[k];
    }
  }
  for(let i=n-1;i>=0;i--){for(let j=i+1;j<n;j++)b[i]-=a[i*n+j]*b[j];b[i]/=a[i*n+i];if(!Number.isFinite(b[i]))throw Error('The edited system did not produce a finite solution.');}
  return b;
}

export function exactSolution(state){
  const clean=preset(state.preset);
  if(state.mode!==clean.mode||state.source!==clean.source||state.amplitude!==1||(!isTransient(state)&&(state.betaX!==clean.betaX||state.betaY!==clean.betaY))||Object.keys(state.matrixEdits).length||Object.keys(state.rhsEdits).length||Object.keys(state.nodeEdits).length||Object.keys(clean.boundaries).some(k=>state.boundaries[k]!==clean.boundaries[k]))return null;
  if(state.preset==='afc-affine')return(x,y)=>1+x+y;
  if(state.preset==='steady-adr')return(x,y)=>Math.sin(Math.PI*x)*Math.sin(Math.PI*y);
  if(state.preset==='forming-layer'){if(state.initial!==clean.initial||state.betaX!==1||state.betaY!==0)return null;const q=expression(layerTrace);return(x,y,t=0)=>(1-Math.exp(-t))*q(x,y,state.epsilon);}
  if(['supg-front','afc-front'].includes(state.preset)){if(state.initial!==clean.initial)return null;const fn=expression(frontInitial),[vx,vy]=velocity(state);return(x,y,t=0)=>Math.exp(-(state.reaction??0)*t)*fn(x-vx*t,y-vy*t);}
  if((state.reaction??0)!==0&&!['reaction-decay','diffusion-reaction','advection-reaction','adr-pulse'].includes(state.preset))return null;
  if(isTransient(state)){
    if(state.initial!==clean.initial)return null;
    const eps=physicalDiffusion(state),[vx,vy]=velocity(state);
    if(state.preset==='moving-front'){const fn=expression(frontInitial);return(x,y,t=0)=>fn(x-vx*t,y-vy*t,eps);}
    if(!['heat-pulse','moving-pulse','reaction-decay','diffusion-reaction','advection-reaction','adr-pulse'].includes(state.preset))return null;
    const[cx,cy,s2]=state.preset==='heat-pulse'?[.5,.5,.000625]:[.25,.4,.0049];
    return(x,y,t=0)=>{const w=s2+4*eps*t;return Math.exp(-(state.reaction??0)*t)*s2/w*Math.exp(-((x-cx-vx*t)**2+(y-cy-vy*t)**2)/w);};
  }
  if(['layer','supg-layer','afc-layer'].includes(state.preset))return (x)=>x-(Math.exp((x-1)/state.epsilon)-Math.exp(-1/state.epsilon))/(-Math.expm1(-1/state.epsilon));
  if(state.preset==='sharp'){const q=expression(sharpQ);return(x,y)=>q(x,y)*Math.sin(Math.PI*y)/state.epsilon;}
  return(x,y)=>Math.sin(Math.PI*x)*Math.sin(Math.PI*y)/state.epsilon;
}

export function interpolate(mesh,u,x,y){
  const n=mesh.n,i=Math.min(n-1,Math.max(0,Math.floor(x*n))),j=Math.min(n-1,Math.max(0,Math.floor(y*n)));
  const candidates=[];
  const radius=mesh.pattern==='shifted'?1:0;
  for(let cy=Math.max(0,j-radius);cy<=Math.min(n-1,j+radius);cy++)for(let cx=Math.max(0,i-radius);cx<=Math.min(n-1,i+radius);cx++)candidates.push(...mesh.triangles.slice(2*(cy*n+cx),2*(cy*n+cx)+2));
  for(const ids of candidates){
    const [a,b,c]=ids.map(id=>mesh.nodes[id]),det=(b[1]-c[1])*(a[0]-c[0])+(c[0]-b[0])*(a[1]-c[1]);
    const l0=((b[1]-c[1])*(x-c[0])+(c[0]-b[0])*(y-c[1]))/det,l1=((c[1]-a[1])*(x-c[0])+(a[0]-c[0])*(y-c[1]))/det,l2=1-l0-l1;
    if(Math.min(l0,l1,l2)>-1e-9)return l0*u[ids[0]]+l1*u[ids[1]]+l2*u[ids[2]];
  }return NaN;
}

export function diagnostics(A,F,b,u,mesh,exact,rows=null){
  let residual=0,normF=0,normA=0,normB=0;
  for(let i=0;i<b.length;i++){
    let sum=-F[i],rowNorm=0;
    if(rows){for(const[j,v]of rows[i]){sum+=v*b[j];rowNorm+=Math.abs(v);}}
    else for(let j=0;j<b.length;j++){sum+=A[i*b.length+j]*b[j];rowNorm+=Math.abs(A[i*b.length+j]);}
    residual=Math.max(residual,Math.abs(sum));normF=Math.max(normF,Math.abs(F[i]));normA=Math.max(normA,rowNorm);normB=Math.max(normB,Math.abs(b[i]));
  }
  let error2=0,integral=0;
  for(const ids of mesh.triangles){
    const area=triangleArea(ids.map(id=>mesh.nodes[id]));
    integral+=area*(u[ids[0]]+u[ids[1]]+u[ids[2]])/3;
    if(exact)for(const [q,w] of quadrature){
    const x=q.reduce((s,v,i)=>s+v*mesh.nodes[ids[i]][0],0),y=q.reduce((s,v,i)=>s+v*mesh.nodes[ids[i]][1],0),uh=q.reduce((s,v,i)=>s+v*u[ids[i]],0);
    error2+=w*(uh-exact(x,y))**2*area;
  }
  }
  return {min:Math.min(...u),max:Math.max(...u),residual,backwardError:residual/(normA*normB+normF||1),l2:exact?Math.sqrt(error2):null,integral,negativeNodes:Array.from(u).filter(v=>v < -1e-10).length};
}

export function compute(state){
  const started=performance.now(),system=assemble(state),{A,K,F,lift,N,mesh,g}=system,delta=addedDiffusion(state);
  let effectiveA=A.map((a,k)=>a+delta*K[k]+system.S[k]),effectiveF=F.map((f,i)=>f+delta*lift[i]+system.supgF[i]);
  for(const[row,id,,,,s]of system.boundaryTerms)effectiveF[row]-=s*g[id];
  const rawB=solveDense(A,F);let afcResult=null,stableB=delta>0||state.stabilization==='supg'?solveDense(effectiveA,effectiveF):rawB.slice();
  if(isAfc(state)){afcResult=solveAfcSteady(state,system);({stableB,effectiveA,effectiveF}=afcResult);}
  const raw=g.slice(),stable=g.slice();mesh.interior.forEach((id,i)=>{raw[id]=rawB[i];stable[id]=stableB[i];});
  const exact=exactSolution(state),reference=exact?mesh.nodes.map(([x,y])=>exact(x,y)):null;
  const metrics={raw:diagnostics(A,F,rawB,raw,mesh,exact),stable:diagnostics(effectiveA,effectiveF,stableB,stable,mesh,exact)};
  if(state.preset==='afc-affine'&&reference&&afcResult){
    afcResult.afc.affine={maxNodalError:Math.max(...stable.map((v,i)=>Math.abs(v-reference[i]))),baselineMaxNodalError:Math.max(...raw.map((v,i)=>Math.abs(v-reference[i])))};
  }
  return {...system,...(afcResult?{afc:afcResult.afc,afcD:afcResult.afcD}:{}),raw,stable,rawB,stableB,reference,effectiveA,effectiveF,metrics,ms:performance.now()-started,
    pe:Math.hypot(...velocity(state))>0?Math.hypot(state.betaX,state.betaY)/(2*state.n*state.epsilon):0,
    effectivePe:Math.hypot(...velocity(state))>0?Math.hypot(state.betaX,state.betaY)/(2*state.n*(state.epsilon+delta)):0};
}
