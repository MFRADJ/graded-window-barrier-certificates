#!/usr/bin/env python3
"""F3 research-check arithmetic (not a theorem prover, not an Arb certificate).
Checks independent formula substitutions and finite examples in the five historically edited proofs.
Uses mpmath 1.3.0. Does not edit input manuscripts or old certificates.
"""
import json, hashlib, math, platform, sys
import mpmath as mp
mp.mp.dps = 75
LOGS=[]

def check(label, cond, details):
    ok=bool(cond)
    LOGS.append(dict(name=label,status='PASS' if ok else 'FAIL', details=str(details)))
    print(f'{"PASS" if ok else "FAIL"} {label}: {details}')


def theta(j,b):
    return mp.jtheta(j,0,mp.exp(-b/2))
def P(b):
    return mp.qp(mp.exp(-b),mp.exp(-b))
def A(b):
    return (theta(3,b)**4+theta(4,b)**4)/(2*P(b)**4)
def B(b):
    return theta(2,b)**4/(2*P(b)**4)
def Psi(b):
    return (theta(3,b)**4+theta(2,b)**4)/(P(b)**4)

# Lemma A.2: identity and its F1 derivative from elementary modular factors.
for d in ['0.4','0.9','2','4','5','8']:
    beta=mp.mpf(d)
    star=4*mp.pi**2/beta
    lgA=mp.log(A(beta))
    rhs=(star-beta)/6-mp.log(2)+mp.log(Psi(star))
    err=abs(lgA-rhs)
    check('A.2 modular logA beta='+d,err < mp.mpf('1e-55'),mp.nstr(err,7))
    def psi(b): return mp.log(Psi(b))
    f1_direct=mp.diff(lambda z:z*mp.log(A(z)),beta)+mp.log(2)
    f1_formula=-beta/3+psi(star)-star*mp.diff(psi,star)
    de=abs(f1_direct-f1_formula)
    check('A.2 F1 derivative beta='+d,de<mp.mpf('1e-54'),mp.nstr(de,7))

# Lemma 6.4: two independent isotropic code weight bounds and thermodynamic envelope.
# For every k and m, the hypothetical saturated sector multiplicity is no larger than
# min(2^m, binom(m,k) 2^k). This tests the inequality used, but not isotropy itself.
for m in (2,5,11):
    for beta_txt in ('1','2','4','5','7'):
        beta=mp.mpf(beta_txt);ell=mp.log(B(beta)/A(beta))
        v=ell+beta*mp.diff(lambda z:mp.log(B(z)/A(z)),beta)
        g=mp.log(A(beta))+beta*mp.diff(lambda z:mp.log(A(z)),beta)
        F1=g+mp.log(2);F2=g+mp.log(1+2*mp.exp(v))
        errs=[]
        for w in range(m+1):
            mult=min(2**m,math.comb(m,w)*2**w)
            deriv=mp.log(mult)/m+g+mp.mpf(w)/m*v
            errs.append(float(deriv-min(F1,F2)))
        check(f'6.4 sector envelope m={m}, beta={beta_txt}',max(errs)<1e-32,f'max excess {max(errs):.3g}')

# Lemma A.4: endpoint box expressions using high precision (NOT rigorous intervals).
def C(b):return B(b)/mp.exp(-b/2)
def deltaA(b):return -mp.diff(A,b)
def U(left,right):
    return mp.log(A(left)) - left*deltaA(right)/A(left)+mp.log(1+2*mp.exp(-left)*C(left)/A(right)*mp.exp(right*deltaA(left)/A(right)))
box=[]
for i in range(10):
    lo=4+mp.mpf(i)/10;hi=lo+mp.mpf(1)/10
    u=U(lo,hi);box.append(u)
check('A.4 ten rational-box evaluations (numeric)',max(box)<0,'largest='+mp.nstr(max(box),14))

# Corollary 6.13: log(2lambda)/lambda decreasing for lambda >=2.
for lam in (2,3,5,10,30,100):
    derivative=(1-mp.log(2*lam))/lam**2
    check('6.13 tolerance monotonic lambda='+str(lam),derivative<0,mp.nstr(derivative,8))

# Lemma 4.4: elementary factors in Chebyshev estimate and uniform size condition.
for d in (mp.mpf('0.1'),mp.mpf('1'),mp.mpf('4')):
    V=mp.mpf('2.7');m=int(mp.ceil(64*V/(9*d*d)))+10
    assert 2*mp.sqrt(m*V)<mp.mpf('3')/4*m*d
    check('4.4 finite-length threshold delta='+mp.nstr(d,3), True,'m='+str(m))

result={'method':'Independent numerical/algebraic spot checks, not complete formal verification','python':sys.version.split()[0], 'mpmath':mp.__version__,'suite_count':len(LOGS),'pass_count':sum(x['status']=='PASS' for x in LOGS),'fail_count':sum(x['status']=='FAIL' for x in LOGS),'tests':LOGS}
from pathlib import Path
path=Path(__file__).with_name('F3_CONTROLES_RESULTATS.json')
path.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print(f"SUMMARY {result['pass_count']}/{result['suite_count']}")
if result['fail_count']:sys.exit(1)
