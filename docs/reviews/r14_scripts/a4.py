import json, numpy as np, pandas as pd, os
S=[json.loads(l) for l in open(os.path.expanduser('~/rtsc_targeted_work/R14-review/smear.jsonl'))]
M=pd.read_csv('data/labels/manifest_log10A.tsv.gz',sep='\t')
L=pd.read_csv('data/labels/labels_all.tsv.gz',sep='\t')
tr=set(M[M.source=='cerq'].id)
A=np.array([[np.nan if a is None else a for a in s['A']] for s in S],float)
lam=np.array([[np.nan if a is None else a for a in s['lam']] for s in S],float)
m=np.array([[np.nan if a is None else a for a in s['m']] for s in S],float)
ids=np.array([s['id'] for s in S]); intr=np.isin(ids,list(tr))
print('n',len(S),'in train',intr.sum())
sig=[0.005*k for k in range(1,11)]
for k in range(10):
    a=A[intr,k]; print(f'sigma {sig[k]:.3f} median A {np.nanmedian(a):.1f} p99 {np.nanpercentile(a,99):.0f} max {np.nanmax(a):.0f} n>=300 {(a>=300).sum()} n>=500 {(a>=500).sum()} n>=1000 {(a>=1000).sum()}  max m {np.nanmax(m[intr,k]):.3g}')
r=A[intr,0]/A[intr,5]; print('A(0.005)/A(0.030) quantiles',np.nanquantile(r,[.05,.25,.5,.75,.95,.99]).round(3))
mx=np.nanmax(A[intr],1); print('envelope max_sigma A: p99',np.nanpercentile(mx,99),'max',np.nanmax(mx),'n>=500',(mx>=500).sum(),'n>=1000',(mx>=1000).sum())
o=np.argsort(-A[intr,0])[:15]
F=dict(zip(L.id,L.formula))
for i in o:
    ii=np.where(intr)[0][i]; print(ids[ii],F.get(ids[ii]),np.round(A[ii],0).tolist(), 'lam',np.round(lam[ii,[0,5]],2).tolist())
# hydrides
H=dict(zip(L.id,L.x_H))
xh=np.array([H.get(i,0) for i in ids]); h=intr&(xh>0)
print('hydrides A(0.005)/A(0.030) median',np.nanmedian(A[h,0]/A[h,5]), 'nonH',np.nanmedian(A[intr&(xh==0),0]/A[intr&(xh==0),5]))
