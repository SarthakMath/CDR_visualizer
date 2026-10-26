import {interpolate,exactSolution,isTransient} from './fem.js';
const palette=[[36,74,128],[67,142,179],[127,188,182],[224,214,161],[239,177,107],[220,116,68],[187,67,59]];
export function color(t,alpha=1){
  t=Math.max(0,Math.min(1,t));const k=t*(palette.length-1),i=Math.min(palette.length-2,Math.floor(k)),v=k-i;
  return `rgba(${palette[i].map((c,j)=>Math.round(c*(1-v)+palette[i+1][j]*v)).join(',')},${alpha})`;
}
function setup(canvas){const r=canvas.getBoundingClientRect(),dpr=Math.min(window.devicePixelRatio||1,2);canvas.width=Math.round(r.width*dpr);canvas.height=Math.round(r.height*dpr);const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);ctx.clearRect(0,0,r.width,r.height);return{ctx,w:r.width,h:r.height};}
function line(ctx,pts,stroke,width=1,dash=[]){ctx.beginPath();pts.forEach((p,i)=>i?ctx.lineTo(p[0],p[1]):ctx.moveTo(p[0],p[1]));ctx.strokeStyle=stroke;ctx.lineWidth=width;ctx.setLineDash(dash);ctx.stroke();ctx.setLineDash([]);}
function text(ctx,str,x,y,align='center',fill='#9aa7b3',font='9px "IBM Plex Mono", monospace'){ctx.fillStyle=fill;ctx.font=font;ctx.textAlign=align;ctx.fillText(str,x,y);}
export function short(v){if(Math.abs(v)<1e-12)return'0';if(Math.abs(v)<0.001||Math.abs(v)>=10000)return v.toExponential(2);return Number(v.toFixed(3)).toString();}

export function meshPlot(canvas,mesh,selected){const{ctx,w,h}=setup(canvas),size=Math.min(w-58,h-30),left=(w-size)/2,top=(h-size)/2-3;const points=mesh.nodes.map(([x,y])=>[left+x*size,top+(1-y)*size]);
  ctx.fillStyle='#f0f5f7';ctx.fillRect(left,top,size,size);
  for(const ids of mesh.triangles)line(ctx,[...ids,ids[0]].map(i=>points[i]),'#99b3bf',.45);
  line(ctx,[[left,top],[left+size,top],[left+size,top+size],[left,top+size],[left,top]],'#688999',1);
  text(ctx,'0',left,top+size+13);text(ctx,'1',left+size,top+size+13);text(ctx,'1',left-11,top+3);text(ctx,'x',left+size+12,top+size+3);text(ctx,'y',left-11,top-6);
  if(selected!==null&&points[selected]){ctx.beginPath();ctx.arc(...points[selected],3.5,0,Math.PI*2);ctx.fillStyle='#dd673b';ctx.fill();}
  return points;
}

export function fieldValues(result,field){if(field==='difference')return result.stable.map((v,i)=>v-result.raw[i]);if(field==='reference')return result.reference||result.raw;return result[field];}
export function surfacePlot(canvas,result,opts){const{ctx,w,h}=setup(canvas),{mesh}=result,values=fieldValues(result,opts.field);
  let lo,hi;if(opts.field==='difference'){hi=result.differenceMax??Math.max(1e-12,...values.map(Math.abs));hi=Math.max(hi,1e-12);lo=-hi;}else if(result.displayRange){[lo,hi]=result.displayRange;}else{lo=Math.min(...result.raw,...result.stable,...(result.reference||[]));hi=Math.max(...result.raw,...result.stable,...(result.reference||[]));}
  const range=hi-lo||1;const norm=v=>(v-lo)/range;let points,project;
  if(opts.view==='heatmap'){
    const size=Math.min(w-115,h-58),left=(w-size)/2-8,top=(h-size)/2-8;
    project=(x,y,z)=>[left+x*size,top+(1-y)*size,0];
    points=mesh.nodes.map(([x,y],i)=>project(x,y,values[i]));
    for(let t=0;t<=4;t++){text(ctx,String(t/4),left+t*size/4,top+size+18);text(ctx,String(t/4),left-12,top+(1-t/4)*size+3,'right');}
    text(ctx,'x',left+size+13,top+size+3);text(ctx,'y',left-13,top-8);
  }else{
    const scale=Math.min(w*.56,h*.82),ca=Math.cos(opts.yaw),sa=Math.sin(opts.yaw),ce=Math.cos(opts.elevation),se=Math.sin(opts.elevation),baseY=h*.55,centerX=w*.48;
    project=(x,y,z)=>{x-=.5;y-=.5;const depth=x*sa+y*ca;return[centerX+(x*ca-y*sa)*scale,baseY+(depth*se-(z-lo)/range*.62*ce)*scale,depth*ce+(z-lo)/range*.62*se];};
    points=mesh.nodes.map(([x,y],i)=>project(x,y,values[i]));
    // Domain grid remains at the bottom of the shared vertical range.
    for(let t=0;t<=10;t++){line(ctx,[project(t/10,0,lo),project(t/10,1,lo)],'#e7edf2',.7);line(ctx,[project(0,t/10,lo),project(1,t/10,lo)],'#e7edf2',.7);}
    for(const t of [0,.25,.5,.75,1]){const p=project(t,0,lo);text(ctx,String(t),p[0]+2,p[1]+19);const q=project(0,t,lo);text(ctx,String(t),q[0]-13,q[1]+3,'right');}
    const px=project(1.12,0,lo),py=project(0,1.12,lo);text(ctx,'x',px[0],px[1]+10,'center','#8b9dac','italic 13px Georgia');text(ctx,'y',py[0]-7,py[1],'center','#8b9dac','italic 13px Georgia');
    const corner=[[0,0],[1,0],[1,1],[0,1]].sort((a,b)=>project(b[0],b[1],lo)[0]-project(a[0],a[1],lo)[0])[0];
    const z0=project(...corner,lo),z1=project(...corner,hi);line(ctx,[z0,z1],'#cfd9e2',.8);
    for(let t=0;t<=4;t++){const p=project(...corner,lo+range*t/4);line(ctx,[[p[0]-3,p[1]],[p[0]+3,p[1]]],'#cfd9e2',.8);text(ctx,short(lo+range*t/4),p[0]+8,p[1]+3,'left');}text(ctx,'u',z1[0],z1[1]-13,'center','#8b9dac','italic 13px Georgia');
  }
  // Sort far faces first for orthographic rendering.
  const ordered=mesh.triangles.map(ids=>({ids,depth:ids.reduce((s,i)=>s+points[i][2],0)})).sort((a,b)=>a.depth-b.depth);
  for(const{ids}of ordered){
    const pts=ids.map(i=>points[i]),v=ids.reduce((s,i)=>s+values[i],0)/3;
    let fill=color(norm(v));
    if(opts.view==='heatmap'){
      const [a,b,c]=pts,dx=b[0]-a[0],dy=b[1]-a[1],ex=c[0]-a[0],ey=c[1]-a[1],det=dx*ey-ex*dy;
      const dv=values[ids[1]]-values[ids[0]],ev=values[ids[2]]-values[ids[0]],gx=(dv*ey-ev*dy)/det,gy=(dx*ev-ex*dv)/det,g2=gx*gx+gy*gy;
      if(g2>1e-24){const start=(lo-values[ids[0]])/g2,end=(hi-values[ids[0]])/g2;fill=ctx.createLinearGradient(a[0]+gx*start,a[1]+gy*start,a[0]+gx*end,a[1]+gy*end);for(let k=0;k<palette.length;k++)fill.addColorStop(k/(palette.length-1),color(k/(palette.length-1)));}
    }
    ctx.beginPath();pts.forEach((p,i)=>i?ctx.lineTo(p[0],p[1]):ctx.moveTo(p[0],p[1]));ctx.closePath();ctx.fillStyle=fill;ctx.fill();
    ctx.strokeStyle=opts.mesh?'rgba(31,59,77,.25)':fill;ctx.lineWidth=opts.mesh?.55:.7;ctx.stroke();
  }
  const section=[];for(let k=0;k<=200;k++){const x=k/200,y=opts.slice,uh=interpolate(mesh,values,x,y);section.push(project(x,y,uh));}line(ctx,section,'rgba(255,255,255,.85)',1.4,[4,3]);
  if(opts.selected!==null&&points[opts.selected]){const p=points[opts.selected];ctx.beginPath();ctx.arc(p[0],p[1],4,0,Math.PI*2);ctx.fillStyle='#fff';ctx.fill();ctx.strokeStyle='#344b5b';ctx.lineWidth=1.5;ctx.stroke();}
  return {points,lo,hi};
}

export function crossPlot(canvas,result,state,slice){const{ctx,w,h}=setup(canvas),left=45,right=20,top=20,bottom=34,cw=w-left-right,ch=h-top-bottom,exact=exactSolution(state),count=800;
  const raw=[],stable=[],reference=[];
  for(let i=0;i<=count;i++){const x=i/count;raw.push(interpolate(result.mesh,result.raw,x,slice));stable.push(interpolate(result.mesh,result.stable,x,slice));if(exact)reference.push(exact(x,slice,result.time||0));}
  let lo=Math.min(0,result.displayRange?.[0]??Infinity,...raw,...stable,...reference),hi=Math.max(result.displayRange?.[1]??-Infinity,...raw,...stable,...reference),pad=(hi-lo||1)*.08;lo-=pad;hi+=pad;
  const xy=(x,y)=>[left+x*cw,top+(hi-y)/(hi-lo)*ch];
  for(let i=0;i<=4;i++){const y=lo+(hi-lo)*i/4,p=xy(0,y);line(ctx,[[left,p[1]],[w-right,p[1]]],'#edf1f5',1);text(ctx,short(y),left-9,p[1]+3,'right',undefined,'8px "IBM Plex Mono", monospace');const x=i/4,q=xy(x,lo);text(ctx,String(x),q[0],h-15,'center',undefined,'8px "IBM Plex Mono", monospace');}
  const draw=(vs,stroke,width,dash=[])=>line(ctx,vs.map((v,i)=>xy(i/count,v)),stroke,width,dash);
  if(reference.length)draw(reference,'#8194a7',1.5,[4,4]);
  if(state.diffusion>0||(isTransient(state)&&state.massType==='lumped'))draw(stable,'#208c88',2);
  draw(raw,'#dd673b',1.8);
  if(result.mesh.n<=32)for(let i=0;i<=result.mesh.n;i++){const x=i/result.mesh.n,p=xy(x,interpolate(result.mesh,result.raw,x,slice));ctx.beginPath();ctx.arc(p[0],p[1],2,0,Math.PI*2);ctx.fillStyle='#dd673b';ctx.fill();}
  text(ctx,'x',w-10,h-15,'center',undefined,'italic 11px Georgia');
}

export function matrixPlot(canvas,result,row,col){const{ctx,w,h}=setup(canvas),{A,N}=result,pad=7,cw=w-2*pad,ch=h-2*pad;let max=0;for(const v of A)max=Math.max(max,Math.abs(v));
  for(let i=0;i<N;i++)for(let j=0;j<N;j++){const v=A[i*N+j];if(Math.abs(v)>1e-14){ctx.fillStyle=v>0?'#d89a7c':'#7d9cad';ctx.fillRect(pad+j/N*cw,pad+i/N*ch,Math.max(.6,cw/N),Math.max(.6,ch/N));}}
  ctx.strokeStyle='#d2613a';ctx.lineWidth=.8;ctx.strokeRect(pad+col/N*cw-1,pad+row/N*ch-1,Math.max(2,cw/N)+2,Math.max(2,ch/N)+2);
}
