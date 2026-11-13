#!/usr/bin/env python3
"""Independent interval cross-check of Appendix A, three regimes (0,4], [4,5], [5,infty).
Uses mpmath.iv at 80 decimal places, not Arb. All sign decisions are on intervals.
Series truncations have explicit analytic upper bounds at q0=e^-1.
This independent checker does not replace the original Arb certificate.
"""
from mpmath import iv, __version__ as mpmath_version
import sys

iv.dps = 80
N = 401

def conv(x, y, n):
    z=[0]*n
    for i, a in enumerate(x[:n]):
        if a:
            for j, b in enumerate(y[:n-i]):
                z[i+j]+=a*b
    return z

def fourth(x):
    t=conv(x,x,N)
    return conv(t,t,N)

th3=[0]*N;th4=[0]*N
for j in range(-30,31):
    if j*j<N:
        th3[j*j]+=1
        th4[j*j]+=1 if j % 2 == 0 else -1
T3,T4=fourth(th3),fourth(th4)
T2=[T3[i]-T4[i] for i in range(N)]
H=(N+1)//2
inv=[1]+[0]*(H-1)
for n in range(1,H):
    for k in range(4):
        for i in range(n,H):
            inv[i]+=inv[i-n]
P4=[0]*N
for i in range(H):
    if 2*i<N:P4[2*i]=inv[i]
Acoeff=conv([(T3[i]+T4[i])//2 for i in range(N)],P4,N)
Bcoeff=conv([T2[i]//2 for i in range(N)],P4,N)
Pcoeff=conv([T3[i]+T2[i] for i in range(N)],P4,N)
a=[Acoeff[2*i] for i in range(H)]
c=[Bcoeff[2*i+1] for i in range(H-1)]
if any(type(x) is not int for x in a+c+Pcoeff):
    raise TypeError("Series coefficients must be exact integers")
if not (a[:6]==[1,28,134,568,1809,5316]):
    raise AssertionError("Exact coefficient invariant failed")
if not (c[:6]==[8,64,288,1024,3152,8704]):
    raise AssertionError("Exact coefficient invariant failed")
if not (all(x>=0 for x in a+c+Pcoeff)):
    raise AssertionError("Exact coefficient invariant failed")

ONE=iv.mpf(1); TWO=iv.mpf(2)
q0=iv.exp(-ONE)
r0=iv.sqrt(q0)
th3_upper=ONE+TWO*r0/(ONE-r0)
th2_upper=TWO*iv.exp(-ONE/8)/(ONE-r0)
invP4_upper=iv.exp(4*q0/(ONE-q0)**2)
GA=th3_upper**4*invP4_upper
GC=GA/r0
GP=(th3_upper**4+th2_upper**4)*invP4_upper

def upper_interval(x):
    return iv.mpf([0,x.b])

def geo(r,n):return r**n/(ONE-r)
def ngeo(r,n):return iv.mpf(n)*r**n/(ONE-r)**2

def adc(beta):
    beta=iv.mpf(beta);q=iv.exp(-beta);r=q/q0
    A=sum((iv.mpf(a[j])*q**j for j in range(len(a))),iv.mpf(0))
    D=sum((iv.mpf(j*a[j])*q**j for j in range(1,len(a))),iv.mpf(0))
    C=sum((iv.mpf(c[j])*q**j for j in range(len(c))),iv.mpf(0))
    return A+upper_interval(GA*geo(r,len(a))),D+upper_interval(GA*ngeo(r,len(a))),C+upper_interval(GC*geo(r,len(c)))

def sprime(y):
    y=iv.mpf(y);s=iv.exp(-y/2);r=s/r0
    total=sum((iv.mpf(j*Pcoeff[j])*s**j/2 for j in range(1,N)), iv.mpf(0))
    return total+upper_interval(GP*ngeo(r,N)/2)

def u_bound(left,right,v_left,v_right):
    a1,d1,c1=v_left
    a2,d2,c2=v_right
    return iv.ln(a1)-left*d2/a1+iv.ln(ONE+TWO*iv.exp(-left)*c1*iv.exp(right*d1/a2)/a2)

def say(name, v, comparator):
    print(name, str(v))
    print('  interval bound certified:',comparator)
    if comparator is not True:
        raise AssertionError(name+' not rigorously established')

print('method: mpmath.iv interval arithmetic, version',mpmath_version,'precision dps',iv.dps)
print('Analytic truncation bounds shared with the Arb checker')
pi=iv.pi
z=pi**2
R=z*(TWO+z)*sprime(z)
lim=4*pi**2/3
say('(a) R(pi^2) < 4*pi^2/3',R,R<lim)
print('  threshold:',lim)

nodes=[iv.mpf(4)+iv.mpf(j)/10 for j in range(11)]
vals=[adc(p) for p in nodes]
ub=[u_bound(nodes[i],nodes[i+1],vals[i],vals[i+1]) for i in range(10)]
worst=max(v.b for v in ub)
say('(b) upper bound on each of 10 boxes of [4,5]',iv.mpf([min(v.a for v in ub),worst]),all(v<0 for v in ub))
print('  largest upper endpoint:',worst)

A5,D5,C5=vals[-1]
e5=A5-ONE
lhs=28*(4-e5)/(ONE+e5)
rhs=TWO*C5*iv.exp(5*D5)
margin=lhs-rhs
say('(c) positive large-beta bracket',margin,margin>0)
print('  lhs:',lhs,'rhs:',rhs)

print('==> MAIN (0,4]/[4,5]/[5,infty) INTERVAL CROSS-CHECK: PASS')
