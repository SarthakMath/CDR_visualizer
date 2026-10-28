import {isSparse} from './matrix.js';
import {factorSparse,solveSparse} from './sparse-solver.js';
// Factor once per distinct step size, reuse throughout the time evolution.
export function factorDense(matrix,n){
  if(isSparse(matrix))return factorSparse(matrix);
  const a=matrix.slice(),pivots=new Int32Array(n);let norm=0;
  for(const v of a)norm=Math.max(norm,Math.abs(v));
  const threshold=norm*Number.EPSILON*n*8;
  for(let k=0;k<n;k++){
    let p=k;for(let i=k+1;i<n;i++)if(Math.abs(a[i*n+k])>Math.abs(a[p*n+k]))p=i;
    if(!Number.isFinite(a[p*n+k])||Math.abs(a[p*n+k])<=threshold)throw Error('The time-step matrix is singular or numerically ill-conditioned. Restore matrix entries or change the time step.');
    pivots[k]=p;if(p!==k)for(let j=0;j<n;j++){const v=a[k*n+j];a[k*n+j]=a[p*n+j];a[p*n+j]=v;}
    for(let i=k+1;i<n;i++){const f=a[i*n+k]/a[k*n+k];a[i*n+k]=f;if(f!==0)for(let j=k+1;j<n;j++)a[i*n+j]-=f*a[k*n+j];}
  }
  const lower=[],upper=[],diagonal=new Float64Array(n);
  for(let i=0;i<n;i++){lower[i]=[];upper[i]=[];diagonal[i]=a[i*n+i];for(let j=0;j<n;j++)if(a[i*n+j]!==0){if(j<i)lower[i].push([j,a[i*n+j]]);else if(j>i)upper[i].push([j,a[i*n+j]]);}}
  return{lower,upper,diagonal,pivots};
}
export function solveFactored(factor,rhs){
  if(factor.kind==='sparse-factor')return solveSparse(factor,rhs);
  const{lower,upper,diagonal,pivots}=factor,b=rhs.slice(),n=b.length;
  for(let k=0;k<n;k++){const p=pivots[k],v=b[k];b[k]=b[p];b[p]=v;}
  for(let i=0;i<n;i++)for(const[j,v]of lower[i])b[i]-=v*b[j];
  for(let i=n-1;i>=0;i--){for(const[j,v]of upper[i])b[i]-=v*b[j];b[i]/=diagonal[i];if(!Number.isFinite(b[i]))throw Error('The time evolution produced a non-finite solution. Check data and matrix edits.');}
  return b;
}

