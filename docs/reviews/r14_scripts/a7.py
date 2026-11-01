import json, numpy as np, pandas as pd, os
from scipy.stats import spearmanr
M=pd.read_csv('data/labels/manifest_log10A.tsv.gz',sep='\t')
rs=[json.loads(l) for l in open(os.path.expanduser('~/rtsc_targeted_work/D1-labels/labels/alexsc.jsonl'))]
d={r['id']:(r.get('dosef'),r.get('debye'),r.get('nsites')) for r in rs}
c=M[M.source=='cerq'].copy()
c['dosef']=[d.get(i,(None,))[0] for i in c.id]; c['debye']=[d.get(i,(None,None))[1] for i in c.id]
c=c.dropna(subset=['dosef']); c['dosef']=c.dosef.astype(float)
c['dpa']=c.dosef/c.nsites
print('n',len(c),'spearman A~dosef',spearmanr(c.A_K,c.dosef)[0],'A~dosef/atom',spearmanr(c.A_K,c.dpa)[0],'lam~dosef/atom',spearmanr(c.lam,c.dpa)[0],'A~debye',spearmanr(c.A_K,c.debye.astype(float))[0])
print('A~x_light*dpa',spearmanr(c.A_K,c.dpa*(1+c.x_light))[0])
