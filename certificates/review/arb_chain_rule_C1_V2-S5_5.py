"""Certificat Arb (independant) des constantes NON certifiees par le dossier :
lemme C.1 (regle de chaine), log2(1+x)<=4x^{1/4}, seuils derives."""
import sys
from flint import arb, ctx
ctx.prec=256
ok=[]
def chk(n,c): ok.append(bool(c)); print(('PASS ' if c else 'FAIL ')+n)
pi=arb.pi(); c=2*pi**2; y=arb('3e-9'); u=16*y
yp=arb('6e-8'); ypp=arb('1.32e-6'); one=1-arb('1e-7')
chk('24y+4u^2/(1-u) < 1e-7', 24*y+4*u**2/(1-u) < arb('1e-7'))
chk('24+64((1-u)^-2-1) < 25', 24+64*((1-u)**-2-1) < 25)
chk('8+64((1-u)^-2-1) < 9', 8+64*((1-u)**-2-1) < 9)
chk('2048(1-u)^-3 < 2049', 2048*(1-u)**-3 < 2049)
chk('F_y>0: 24-64((1-u)^-2-1)>0', 24-64*((1-u)**-2-1) > 0)
chk('G_y<0: -8+64((1-u)^-2-1)<0', -8+64*((1-u)**-2-1) < 0)
g2=25*ypp/one+(2049/one+625/one**2)*yp**2
r2=34*ypp/one+(4098/one+(625+81)/one**2)*yp**2
chk("|g''| bound < 1e-4  [%s]"%g2.str(5), g2 < arb('1e-4'))
chk("|r''| bound < 1e-4  [%s]"%r2.str(5), r2 < arb('1e-4'))
chk("|g'|<=25 y'/(1-1e-7) < 1e-4", 25*yp/one < arb('1e-4'))
chk("|r'|<=34 y'/(1-1e-7) < 1e-4", 34*yp/one < arb('1e-4'))
chk("|g|,|r| <= -log(1-2e-7) < 1e-4", -(1-arb('2e-7')).log() < arb('1e-4'))
# log2(1+x) <= 4 x^{1/4} : f'(x)>0 <=> x^{-3/4}(1+x) ln2 > 1, min en x=3
chk('4*3^{-3/4} ln2 > 1', 4*arb(3)**(arb(-3)/4)*arb(2).log() > 1)
print("VERDICT %d/%d"%(sum(ok),len(ok))); sys.exit(0 if all(ok) else 1)
