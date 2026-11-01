import numpy as np, pandas as pd, json, os
M=pd.read_csv('data/labels/manifest_log10A.tsv.gz',sep='\t')
L=pd.read_csv('data/labels/labels_all.tsv.gz',sep='\t')
c=M[M.source=='cerq'].groupby(['formula','nsites']).A_K.median(); j=M[M.source=='jarvis'].groupby(['formula','nsites']).A_K.median()
b=pd.concat([c.rename('c'),j.rename('j')],axis=1).dropna(); r=b.j/b.c
print('same formula+nsites in both:',len(b),'A_jar/A_cerq quantiles',r.quantile([.05,.25,.5,.75,.95]).round(2).tolist(),'spearman',b.corr('spearman').iloc[0,1])
kc=L[(L.source=='cerq')&L.id.isin(M.id)]
print('kept cerq wmin in [-30,-5):',((kc.wmin_cm<-5)&(kc.wmin_cm>=-30)).sum())
G=json.load(open(os.path.expanduser('~/rtsc_targeted_work/D1-labels/gen/proxies_comp.json')))
for k,per in G['per_structure'].items():
    kn=np.array([p['knn_log10A'] for p in per if not p['in_train']]); print(k.rsplit('/',1)[-1],'novel n',len(kn),'knn mean K',round(10**kn.mean(),1))
