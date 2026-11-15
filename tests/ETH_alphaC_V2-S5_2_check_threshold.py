#!/usr/bin/env python3
"""Independent integer-series check *inside* mh>=4000 for vacuum code.
The series is computed using divisor-sum and logarithmic-derivative recurrences;
it does not call Arb/python-flint or the V2-S4 review scripts.
"""
import math, json
from pathlib import Path

def sigma1(N):
    a=[0]*(N+1)
    for d in range(1,N+1):
        for k in range(d,N+1,d):a[k]+=d
    return a

def d4_theta(N):
    # theta_{D4}(q)=1+24q+24q^2+96q^3+24q^4+...
    # r_D4(n) = 24*sum_{d|n, d odd}d for n>=1
    a=[0]*(N+1);a[0]=1
    for d in range(1,N+1,2):
        for k in range(d,N+1,d):a[k]+=24*d
    assert a[:5]==[1,24,24,96,24]
    return a

def nth_power(F,m,N):
    # n C_n = sum_{k=1}^n ((m+1)*k-n) F_k C_{n-k}
    c=[0]*(N+1);c[0]=1
    for n in range(1,N+1):
        ss=sum(((m+1)*k-n)*F[k]*c[n-k] for k in range(1,n+1))
        assert ss%n==0
        c[n]=ss//n
    return c

def eta_power(r,N,divisors):
    # P(q)^{-r} = exp(r sum_{n>=1}sigma1(n) q^n/n)
    p=[0]*(N+1);p[0]=1
    for n in range(1,N+1):
        ss=r*sum(divisors[k]*p[n-k] for k in range(1,n+1))
        assert ss%n==0
        p[n]=ss//n
    return p

def main(m=20,h=300):
    assert h>=15*m and m*h>=4000
    N=4*h
    sig=sigma1(N)
    T=d4_theta(N)
    top=nth_power(T,m,N)
    bottom=eta_power(4*m,N,sig)
    d=[0]*(N+1)
    for i,x in enumerate(top):
        if x:
            for j in range(N-i+1):
                d[i+j]+=x*bottom[j]
    assert d[0]==1 and d[1]==28*m
    c=[];z=0
    for x in d:
        z+=x;c.append(z)
    def window(a,b):
        # a,b are endpoint halves, inclusive; step 1/2
        lo=max(0,math.ceil(a/2))
        hi=min(N,math.floor(b/2))
        return 0 if hi<lo else c[hi]-(c[lo-1] if lo else 0)
    lam=d[h].bit_length()-1
    best=(-1,None)
    for t2 in range(-4*h,4*h+1):
        x=window(t2,t2+4*h)
        y=window(-t2,4*h-t2)
        value=math.log2(max(1,x))+math.log2(max(1,y))
        if value>best[0]:best=(value,t2/2)
    rho=best[0]/(4*lam)
    k=0;K=math.pi/math.log(2)*math.sqrt(8/3)
    X=K*math.sqrt(m*h)
    e1=4*m/3+24*(m*h)**.25+1;e2=e1+1
    bounds={
       'lower':1/math.sqrt(2)-.19*math.sqrt(m/h)-1.7*(m*h)**(-.25),
       'upper':1/math.sqrt(2)+.4*math.sqrt(m/h)+5*(m*h)**(-.25),
       'upper_strict':(1/math.sqrt(2))*(X+k/math.sqrt(2))/(X-e2),
       'X_minus_e2':X-e2,
    }
    assert rho>=bounds['lower'] and rho<=bounds['upper'] and rho<=bounds['upper_strict']
    assert lam>=bounds['X_minus_e2']
    out={'m':m,'h':h,'mh':m*h,'code':'vacuum','lambda_h':lam,'rho_win':rho,
         'maximizer_t':best[1],'bounds':bounds,'checks':'PASS','method':'exact integer series + float comparisons'}
    p=Path(__file__).with_name('ETH_alphaC_V2-S5_2_threshold.json')
    p.write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
