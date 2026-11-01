import numpy as np, pandas as pd
from scipy.stats import spearmanr
M=pd.read_csv('data/labels/manifest_log10A.tsv.gz',sep='\t')
la=np.log10(M.A_K); lm=np.log10(M.lam_w2_K2); ll=np.log10(M.lam); w=2*M.A_K/M.lam; lw=np.log10(w)
for nm,x,targets in [('A',la,[1299.44,2060]),('m',lm,[2.696e6,6e6])]:
    print(nm,'mean',x.mean().round(3),'sd',x.std().round(3),'max',10**x.max(),[ (t, round((np.log10(t)-x.mean())/x.std(),2)) for t in targets])
print('var logA',la.var(),'var loglam',ll.var(),'var log<w>',lw.var(),'2cov',2*np.cov(ll,lw)[0,1])
print('spearman A vs x_H',spearmanr(M.A_K,M.x_H)[0],'A vs x_light',spearmanr(M.A_K,M.x_light)[0],'A vs <w>',spearmanr(M.A_K,w)[0],'A vs lam',spearmanr(M.A_K,M.lam)[0], 'm vs x_light',spearmanr(M.lam_w2_K2,M.x_light)[0])
# bound Tc
tb=np.minimum(M.A_K/4.33148, 0.18273*np.sqrt(M.lam_w2_K2)); print('bound-Tc max',tb.max(),'p99',np.percentile(tb,99))
print('A-binding fraction', np.mean(M.A_K/4.33148 < 0.18273*np.sqrt(M.lam_w2_K2)))
h=M[M.x_H>0]; print('hydrides n',len(h),'A median',h.A_K.median(),'nonH',M[M.x_H==0].A_K.median())
