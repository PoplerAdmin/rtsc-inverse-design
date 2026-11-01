import pandas as pd, numpy as np, json, os
from pymatgen.core import Composition
from scipy.stats import spearmanr, mannwhitneyu
D='data/labels/'
M=pd.read_csv(D+'manifest_log10A.tsv.gz',sep='\t')
def vec(f):
    c=Composition(f); v=np.zeros(100)
    for e,a in c.items(): v[e.Z]=a/c.num_atoms
    return v
M['vec']=M.formula.map(vec); M['chemsys']=M.formula.map(lambda f:'-'.join(sorted(e.symbol for e in Composition(f).elements)))
M['anon']=M.formula.map(lambda f:Composition(f).anonymized_formula)
tr=M[M.split=='train']; te=M[M.split=='test']
X=np.stack(tr.vec); y=np.log10(tr.A_K.values)
def knn(v,X,y,k=5):
    d=np.abs(X-v).sum(1); i=np.argsort(d)[:k]; return np.median(y[i]), d[i[0]]
P=np.array([knn(v,X,y) for v in te.vec]); yt=np.log10(te.A_K.values)
print('test kNN spearman',spearmanr(P[:,0],yt)[0])
dn=P[:,1]; print('nearest-train L1 dist quantiles',np.quantile(dn,[.1,.25,.5,.75,.9]).round(3))
for lo,hi in [(0,.25),(.25,.5),(.5,2.1)]:
    m=(dn>=lo)&(dn<hi); print(f'dist[{lo},{hi}) n={m.sum()} spearman={spearmanr(P[m,0],yt[m])[0]:.3f} mae={np.abs(P[m,0]-yt[m]).mean():.3f}')
# chemsys overlap
trcs=set(tr.chemsys); print('test rows whose chemsys appears in train', te.chemsys.isin(trcs).mean())
# test whose same anon formula + chemsys differs by one element substitution: count train rows sharing anon formula and >= n-1 elements
trd=tr.groupby('anon')
def onesub(r):
    els=set(r.chemsys.split('-'))
    if r.anon not in trd.groups: return False
    for cs in tr.loc[trd.groups[r.anon],'chemsys']:
        if len(els & set(cs.split('-')))>=len(els)-1: return True
    return False
print('test rows with a 1-element-substitution analog in train', te.apply(onesub,axis=1).mean())
# grouped-by-chemsys split kNN re-validation (train/test from all rows)
rng=np.random.default_rng(0); cs=np.array(sorted(set(M.chemsys))); rng.shuffle(cs); tcs=set(cs[:len(cs)//10])
A=M[~M.chemsys.isin(tcs)]; B=M[M.chemsys.isin(tcs)]
XA=np.stack(A.vec); yA=np.log10(A.A_K.values)
PB=np.array([knn(v,XA,yA)[0] for v in B.vec]); yB=np.log10(B.A_K.values)
print('chemsys-disjoint kNN: n',len(B),'spearman',spearmanr(PB,yB)[0],'mae',np.abs(PB-yB).mean(),'base',np.abs(yB-yA.mean()).mean())
# Leave-element-family: hold out all hydrides rows top? skip
# kNN ceiling: max kNN possible
print('max train log10A', y.max(), 10**y.max())
# generated sets
G=json.load(open(os.path.expanduser('~/rtsc_targeted_work/D1-labels/gen/proxies_comp.json')))
ks=list(G['per_structure'].keys())
metals=set('Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Cs Ba La Ce Pr Nd Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi'.split())
ref=None
for k in ks:
    per=G['per_structure'][k]
    kn=np.array([p['knn_log10A'] for p in per]); it=np.array([p['in_train'] for p in per])
    mf=np.array([not any(e.symbol in metals for e in Composition(p['formula']).elements) for p in per])
    lo,hi=np.percentile([np.mean(np.random.default_rng(i).choice(kn,len(kn))) for i in range(2000)],[2.5,97.5])
    print(k.rsplit('/',1)[-1], 'n',len(per),'knn mean K',round(10**kn.mean(),1),'95%CI of mean (K)',round(10**lo,1),round(10**hi,1),'in_train',it.mean().round(3),'metalfree',mf.sum(), 'knn mean metal-only K', round(10**kn[~mf].mean(),1), 'max knn K',round(10**kn.max(),1))
names=[k.rsplit('/',1)[-1] for k in ks]; print(names)
per={k.rsplit('/',1)[-1]:G['per_structure'][k] for k in ks}
a=np.array([p['knn_log10A'] for p in per['A100']]); b=np.array([p['knn_log10A'] for p in per['A300']]); c=np.array([p['knn_log10A'] for p in per['A1500']])
print('A300 vs A100 knn MWU p', mannwhitneyu(b,a).pvalue, 'A1500 vs A300', mannwhitneyu(c,b).pvalue)
for key in ['x_light','x_H']:
    a=np.array([p[key] for p in per['A100']]); b=np.array([p[key] for p in per['A300']]); print(key,'A300 vs A100 p',mannwhitneyu(b,a).pvalue)
# count tests in table
PH=json.load(open(os.path.expanduser('~/rtsc_targeted_work/D1-labels/gen/proxies_phonon.json')))
for k,v in PH.items(): print(k.rsplit('/',1)[-1], {kk:(round(vv,4) if isinstance(vv,float) else vv) for kk,vv in v.items() if kk!='top_proxy'})
