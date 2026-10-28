// ILU(0)-preconditioned Krylov solves. Only preconditioning pivots are guarded;
// convergence is always checked against the unmodified matrix and right side.
import {rows} from './matrix.js';
const dot=(a,b)=>{let s=0;for(let i=0;i<a.length;i++)s+=a[i]*b[i];return s;};
const norm=a=>Math.sqrt(dot(a,a));
function multiply(f,x){const y=new Float64Array(f.n);for(let i=0;i<f.n;i++)for(let k=f.ptr[i];k<f.ptr[i+1];k++)y[i]+=f.val[k]*x[f.col[k]];return y;}
function precondition(f,b){const x=b.slice();for(let i=0;i<f.n;i++)for(let k=f.ptr[i];k<f.diag[i];k++)x[i]-=f.lu[k]*x[f.col[k]];for(let i=f.n-1;i>=0;i--){for(let k=f.diag[i]+1;k<f.ptr[i+1];k++)x[i]-=f.lu[k]*x[f.col[k]];x[i]/=f.lu[f.diag[i]];}return x;}
export function factorSparse(a){
 const n=a.n,rr=rows(a).map((r,i)=>{if(!r.some(([j])=>j===i))r.push([i,0]);return r.sort((a,b)=>a[0]-b[0]);}),ptr=new Int32Array(n+1),diag=new Int32Array(n);
 for(let i=0;i<n;i++)ptr[i+1]=ptr[i]+rr[i].length;
 const col=new Int32Array(ptr[n]),val=new Float64Array(ptr[n]);
 for(let i=0;i<n;i++)rr[i].forEach(([j,v],k)=>{col[ptr[i]+k]=j;val[ptr[i]+k]=v;if(i===j)diag[i]=ptr[i]+k;});
 const lu=val.slice(),position=new Int32Array(n).fill(-1);
 for(let i=0;i<n;i++){
  let scale=0;for(let k=ptr[i];k<ptr[i+1];k++){position[col[k]]=k;scale=Math.max(scale,Math.abs(val[k]));}
  if(scale===0)throw Error('The matrix has an empty row and is singular. Restore entries or choose a well-posed problem.');
  for(let k=ptr[i];k<diag[i];k++){const j=col[k],v=lu[k]/lu[diag[j]];lu[k]=v;for(let q=diag[j]+1;q<ptr[j+1];q++){const at=position[col[q]];if(at>=0)lu[at]-=v*lu[q];}}
  if(Math.abs(lu[diag[i]])<scale*1e-10)lu[diag[i]]=(lu[diag[i]]<0?-1:1)*scale*1e-10;
  for(let k=ptr[i];k<ptr[i+1];k++)position[col[k]]=-1;
 }
 return{kind:'sparse-factor',n,ptr,col,val,lu,diag,last:null,iterations:0};
}
function residual(f,b,x){const ax=multiply(f,x);return b.map((v,i)=>v-ax[i]);}
export function solveSparse(f,b){
 const n=f.n,bnorm=norm(b);if(bnorm===0){f.last=new Float64Array(n);return f.last.slice();}
 const tol=bnorm*1e-11,maxIterations=2000;let x=f.last?.slice()??new Float64Array(n),r=residual(f,b,x),shadow=r.slice(),p=new Float64Array(n),v=new Float64Array(n),rhoOld=1,alpha=1,omega=1;
 const finish=()=>{const err=norm(residual(f,b,x));if(Number.isFinite(err)&&err<=tol*5){f.last=x.slice();return x;}return null;};
 if(norm(r)<=tol)return finish();
 for(let it=0;it<maxIterations;it++){
  const rho=dot(shadow,r);if(!Number.isFinite(rho)||Math.abs(rho)<1e-290||Math.abs(omega)<1e-290)break;
  const beta=(rho/rhoOld)*(alpha/omega);for(let i=0;i<n;i++)p[i]=r[i]+beta*(p[i]-omega*v[i]);
  const ph=precondition(f,p);v=multiply(f,ph);const den=dot(shadow,v);if(!Number.isFinite(den)||Math.abs(den)<1e-290)break;alpha=rho/den;
  const s=r.map((z,i)=>z-alpha*v[i]);
  if(norm(s)<=tol){for(let i=0;i<n;i++)x[i]+=alpha*ph[i];const answer=finish();if(answer){f.iterations=it+1;return answer;}break;}
  const sh=precondition(f,s),t=multiply(f,sh),tt=dot(t,t);if(!Number.isFinite(tt)||tt<1e-290)break;omega=dot(t,s)/tt;
  for(let i=0;i<n;i++){x[i]+=alpha*ph[i]+omega*sh[i];r[i]=s[i]-omega*t[i];}
  if(norm(r)<=tol){const answer=finish();if(answer){f.iterations=it+1;return answer;}break;}
  rhoOld=rho;
 }
 // Restarted right-preconditioned GMRES is a fallback for BiCGSTAB breakdown.
 x=new Float64Array(n);
 for(let cycle=0;cycle<50;cycle++){
  r=residual(f,b,x);const rn=norm(r);if(rn<=tol){const answer=finish();if(answer)return answer;}
  const V=[r.map(v=>v/rn)],Z=[],H=Array.from({length:31},()=>new Float64Array(30)),cs=new Float64Array(30),sn=new Float64Array(30),g=new Float64Array(31);g[0]=rn;let m=0;
  for(let j=0;j<30;j++){
   const z=precondition(f,V[j]);Z.push(z);const w=multiply(f,z);
   // Reorthogonalization limits loss of orthogonality on difficult edited systems.
   for(let pass=0;pass<2;pass++)for(let i=0;i<=j;i++){const d=dot(w,V[i]);H[i][j]+=d;for(let k=0;k<n;k++)w[k]-=d*V[i][k];}
   H[j+1][j]=norm(w);const breakdown=H[j+1][j]<1e-290;if(!breakdown)V.push(w.map(v=>v/H[j+1][j]));
   for(let i=0;i<j;i++){const a=H[i][j],b=H[i+1][j];H[i][j]=cs[i]*a+sn[i]*b;H[i+1][j]=-sn[i]*a+cs[i]*b;}
   const d=Math.hypot(H[j][j],H[j+1][j]);if(!Number.isFinite(d)||d===0)break;
   cs[j]=H[j][j]/d;sn[j]=H[j+1][j]/d;H[j][j]=d;g[j+1]=-sn[j]*g[j];g[j]*=cs[j];m=j+1;
   if(Math.abs(g[j+1])<=tol||breakdown)break;
  }
  if(!m)break;const y=g.slice(0,m);for(let i=m-1;i>=0;i--){for(let j=i+1;j<m;j++)y[i]-=H[i][j]*y[j];y[i]/=H[i][i];}
  for(let j=0;j<m;j++)for(let i=0;i<n;i++)x[i]+=y[j]*Z[j][i];
  const answer=finish();if(answer){f.iterations=maxIterations+(cycle+1)*30;return answer;}
 }
 throw Error('The sparse linear solver did not converge to its residual tolerance. Reduce the mesh, adjust the problem or stabilization, or clear matrix overrides. No unconverged linear solution is displayed.');
}
