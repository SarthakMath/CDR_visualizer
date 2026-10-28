// Plain structured-cloneable sparse rows above the original dense mesh limit.
export const isSparse=a=>a?.kind==='sparse';
export function matrix(n,sparse=n>961){return sparse?{kind:'sparse',n,rows:Array.from({length:n},()=>new Map())}:new Float64Array(n*n);}
export function get(a,i,j){return isSparse(a)?a.rows[i].get(j)||0:a[i*Math.sqrt(a.length)+j];}
export function set(a,i,j,v){if(isSparse(a)){if(v===0)a.rows[i].delete(j);else a.rows[i].set(j,v);}else a[i*Math.sqrt(a.length)+j]=v;}
export function add(a,i,j,v){if(v!==0)set(a,i,j,get(a,i,j)+v);}
export function rows(a,n=isSparse(a)?a.n:Math.sqrt(a.length)){return isSparse(a)?a.rows.map(r=>Array.from(r)):Array.from({length:n},(_,i)=>{const r=[];for(let j=0;j<n;j++)if(a[i*n+j]!==0)r.push([j,a[i*n+j]]);return r;});}
export function clone(a){return isSparse(a)?{kind:'sparse',n:a.n,rows:a.rows.map(r=>new Map(r))}:a.slice();}
export function combine(...terms){const n=isSparse(terms[0][0])?terms[0][0].n:Math.sqrt(terms[0][0].length),out=matrix(n,terms.some(([a])=>isSparse(a)));for(const[a,scale]of terms){if(!scale)continue;if(isSparse(a)){a.rows.forEach((r,i)=>{for(const[j,v]of r)add(out,i,j,scale*v);});}else if(isSparse(out)){for(let i=0;i<n;i++)for(let j=0;j<n;j++)add(out,i,j,scale*a[i*n+j]);}else for(let k=0;k<a.length;k++)out[k]+=scale*a[k];}return out;}
export function entries(a){const out=[];rows(a).forEach((r,i)=>{for(const[j,v]of r)out.push([i,j,v]);});return out;}
