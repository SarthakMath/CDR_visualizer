// One content version for the entire module graph, including worker imports.
// Run before serving or deployment; --check verifies without modifying files.
import {readFileSync,writeFileSync,readdirSync} from 'node:fs';
import {createHash} from 'node:crypto';
const dist=new URL('../dist/',import.meta.url);
const names=readdirSync(dist).filter(n=>/\.(js|css|html)$/.test(n)).sort();
const local=/(['"])(\.\/[\w.-]+\.(?:js|css))(?:\?v=[\w-]+)?\1/g;
const sources=names.map(name=>[name,readFileSync(new URL(name,dist),'utf8').replace(local,(_,q,path)=>q+path+q)]);
const hash=createHash('sha256');for(const[name,source]of sources)hash.update(name+'\0'+source+'\0');
const version=hash.digest('hex').slice(0,16);let stale=false;
for(const[name,source]of sources){
 const versioned=source.replace(local,(_,q,path)=>`${q}${path}?v=${version}${q}`),file=new URL(name,dist);
 if(readFileSync(file,'utf8')!==versioned){stale=true;if(!process.argv.includes('--check'))writeFileSync(file,versioned);}
}
if(stale&&process.argv.includes('--check')){console.error('Asset versions are stale. Run node scripts/version-assets.mjs from CDR_visualizer.');process.exitCode=1;}
else console.log(`Asset version ${version}${process.argv.includes('--check')?' verified':''}`);
