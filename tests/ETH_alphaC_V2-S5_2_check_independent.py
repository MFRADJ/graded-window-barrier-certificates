#!/usr/bin/env python3
"""Independent small-scale exact checks for the alphaC V2-S5.2 manuscript.

Python standard library only. All Hilbert-series coefficients and dimensions are
exact integers. Logarithms/real constants are floating-point checks, not
interval certificates. This script does not validate the full analytic proof.
"""
import itertools
import math
import json
from pathlib import Path

KAPPA = math.pi / math.log(2) * math.sqrt(8/3)
INV_SQ2 = 1 / math.sqrt(2)


def conv(p, q, L):
    r = [0] * L
    for i, x in enumerate(p[:L]):
        if x:
            for j, y in enumerate(q[:L - i]):
                if y:
                    r[i + j] += x * y
    return r


def power(p, m, L):
    out = [1] + [0] * (L-1)
    while m:
        if m & 1:
            out = conv(out, p, L)
        m //= 2
        if m:
            p = conv(p, p, L)
    return out


def characters(L):
    # x = q^(1/2): theta_3 = sum_{n in Z} x^(n^2)
    theta3, theta4 = [0] * L, [0] * L
    for n in range(math.isqrt(L-1) + 1):
        if n*n >= L:
            break
        theta3[n*n] += 1 if n == 0 else 2
        theta4[n*n] += (1 if n == 0 else 2*(-1)**n)
    n3, n4 = power(theta3, 4, L), power(theta4, 4, L)
    # 1/P(q)^4 with P(q)=prod_(n>=1)(1-q^n), q=x^2
    invp = [1] + [0]*(L-1)
    for size in range(2, L, 2):
        for _ in range(4):
            for j in range(size, L):
                invp[j] += invp[j-size]
    N = [(a+b)//2 for a,b in zip(n3,n4)]
    assert all((a+b)%2==0 for a,b in zip(n3,n4))
    D = [(a-b)//2 for a,b in zip(n3,n4)]
    assert all((a-b)%2==0 for a,b in zip(n3,n4))
    A = conv(N, invp, L)
    B = conv(D, invp, L)
    assert A[:7] == [1,0,28,0,134,0,568] and B[:7] == [0,8,0,64,0,288,0], (A[:7],B[:7])
    return A, B


def sector(A,B,m,w,L):
    return conv(power(A,m-w,L),power(B,w,L),L)


def cost(m, counts, h2, A, B):
    L = len(A)
    d = [0]*L
    min_sector_margin = float('inf')
    for w, count in enumerate(counts):
        if not count:
            continue
        sec = sector(A,B,m,w,L)
        for i in range(L):
            d[i] += count*sec[i]
        e2=h2
        cs = sum(sec[:e2+1]); E=e2/2
        lower = KAPPA*math.sqrt(m*E)-4*m/3-20*(m*E)**.25-1
        min_sector_margin=min(min_sector_margin,math.log2(cs)-lower)
        assert all(sec[j] <= sec[j+2] for j in range(w,len(sec)-2,2)), 'non-monotonic sector'
    if d[h2] < 2:
        return None
    C = list(itertools.accumulate(d))
    h=h2/2
    lh=d[h2].bit_length()-1
    def win(left,right):
        left=max(0,left)
        if right<left:return 0
        return C[right]-(C[left-1] if left else 0)
    best=-1
    best_t=None
    for t2 in range(-2*h2,2*h2+1):
        wd1=win(t2,t2+2*h2)
        wd2=win(-t2,2*h2-t2)
        val=math.log2(max(1,wd1))+math.log2(max(1,wd2))
        if val>best:
            best=val;best_t=t2/2
    rho=best/(4*lh)
    k=(sum(counts)).bit_length()-1
    X=KAPPA*math.sqrt(m*h)
    e1=4*m/3+24*(m*h)**.25+1
    e2=e1+1
    low=INV_SQ2-0.19*math.sqrt(m/h)-1.7*(m*h)**(-.25)
    upper=INV_SQ2+0.4*math.sqrt(m/h)+5*(m*h)**(-.25)
    exact_center=math.log2(C[2*h2])/(2*math.log2(C[h2]))
    assert rho>=exact_center-1e-9
    assert rho>=low-1e-9
    if X>e2:
        assert rho<= INV_SQ2*(X+k/math.sqrt(2))/(X-e2)+1e-9
    if m*h>=4000:
        assert rho<=upper+1e-9
    assert lh>=X-e2-1e-9
    return {"m":m,"k":k,"h":h,"lambda_h":lh,"rho":round(rho,8),"t_star":best_t,
            "min_sector_margin":round(min_sector_margin,5),"margin_lower":round(rho-low,6)}


def main():
    samples=[]
    # Counts of weights of real additive subgroups:
    # vacuum {0}, and the coordinatewise rank-m lagrangian <(1,0)_i>.
    testcases=[]
    for m in (1,2,3):
        hvals=[15*m, 15*m+.5, 18*m]
        zero=[1]+[0]*m
        lagrangian=[math.comb(m,w) for w in range(m+1)]
        for h in hvals:
            for name,counts in [("vacuum",zero),("lagrangian",lagrangian)]:
                if (h%1) and name=="vacuum":
                    continue
                testcases.append((m,name,counts,int(2*h)))
    L=8*max(h2 for _,_,_,h2 in testcases)//2+3 # 4h energy = 8h units in x
    A,B=characters(L)
    for m,name,counts,h2 in testcases:
        s=cost(m,counts,h2,A,B)
        if s: s["code"]=name; samples.append(s)
    # Analogous real upper bound and safe inequalities
    inequalities={
      "20*2^(.25)<24": 20*2**.25 <24,
      "lower sqrt coefficient<0.19": (2/3+1/math.sqrt(2))/KAPPA <.19,
      "lower fourth-root coefficient<1.7":12.5/KAPPA<1.7,
      "b<0.5":((4/3)/math.sqrt(15)+26*4000**(-.25))/KAPPA<.5,
      "upper sqrt coefficient<0.4":(1+4*math.sqrt(2)/3)/KAPPA<.4,
      "upper fourth-root coefficient<5":26*math.sqrt(2)/KAPPA<5
    }
    assert all(inequalities.values())
    result={"status":"PASS (small exact tests; real inequalities evaluated in float)",
            "cases":len(samples),"kappa":KAPPA,"inequalities":inequalities,"samples":samples}
    p=Path(__file__).with_name("ETH_alphaC_V2-S5_2_checks.json")
    p.write_text(json.dumps(result,indent=2))
    print('RESULT:',result['status'],'cases=',len(samples),'kappa=',KAPPA)
    for x in samples:print(x)
    print('WROTE:',p)

if __name__=='__main__': main()
