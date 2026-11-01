import os, gzip, numpy as np, json, sys
R=os.path.expanduser('~/rtsc_targeted_work/D1-labels/raw/mc_3kbt5/extract')
K=11.604518*13605.693
out=open(sys.argv[1],'w')
for b in sorted(os.listdir(R)):
    bd=os.path.join(R,b)
    if not os.path.isdir(bd): continue
    for m in sorted(os.listdir(bd)):
        d=os.path.join(bd,m); rec={'id':m.rsplit('_',1)[-1]}
        A=[];L=[];M2=[];N=[]
        for k in range(1,11):
            fn=os.path.join(d,f'a2F.dos{k}.gz')
            try:
                a=np.loadtxt(gzip.open(fn,'rt'),comments=['#',' lambda','lambda'],usecols=(0,1))
            except Exception: A.append(None);L.append(None);M2.append(None);N.append(None);continue
            w=a[:,0]*K; f=a[:,1]; p=w>0
            A.append(float(np.trapezoid(np.clip(f[p],0,None),w[p])))
            L.append(float(2*np.trapezoid(np.clip(f[p],0,None)/w[p],w[p])))
            M2.append(float(2*np.trapezoid(np.clip(f[p],0,None)*w[p],w[p])))
            q=w<0; N.append(float(np.trapezoid(np.abs(f[q]),w[q])) if q.sum()>1 else 0.0)
        rec.update(A=A,lam=L,m=M2,neg=N)
        try:
            rec['lambda_file']=open(os.path.join(d,'lambda')).read()[:800]
        except Exception: pass
        out.write(json.dumps(rec)+'\n')
print('DONE')
