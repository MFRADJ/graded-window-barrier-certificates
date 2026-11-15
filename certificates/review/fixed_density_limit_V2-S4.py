#!/usr/bin/env python3
"""V2-S4 -- limite a densite fixee (m -> inf, h/m -> r, k/m -> theta) pour les codes
a caractere produit Phi^m : Phi = A (S={0}), A+B (S=<(1,0)_i>, k=m), A+3B (S=K^m, k=2m).
Exposant de Chernoff/grandes deviations : phi(r) = (1/ln2) inf_beta [log Phi(e^-beta) + beta r],
rho_inf(r) = phi(2r) / (2 phi(r)). Heuristique standard (non certifiee) servant de repere."""
import mpmath as mp
mp.mp.dps = 30


def A(q):
    t3 = mp.jtheta(3, 0, mp.sqrt(q)); t4 = mp.jtheta(4, 0, mp.sqrt(q))
    return (t3 ** 4 + t4 ** 4) / (2 * mp.qp(q) ** 4)


def B(q):
    return mp.jtheta(2, 0, mp.sqrt(q)) ** 4 / (2 * mp.qp(q) ** 4)


PHI = {'{0} (theta=0)': lambda q: A(q), '<(1,0)> (theta=1)': lambda q: A(q) + B(q),
       'K^m (theta=2)': lambda q: A(q) + 3 * B(q)}
KAPPA = mp.pi / mp.log(2) * mp.sqrt(mp.mpf(8) / 3)


def phi(F, r):
    f = lambda b: mp.log(F(mp.exp(-b))) + b * r
    b0 = mp.sqrt(2 * mp.pi ** 2 / 3 / r)
    b = mp.findroot(lambda b: mp.diff(f, b), b0)
    return f(b) / mp.log(2)


if __name__ == '__main__':
    print('%-20s %6s %10s %10s %12s %14s' % ('code', 'r', 'phi(r)', 'rho_inf', 'rho-1/r2', 'formule appr.'))
    for name, F in PHI.items():
        th = {'{0} (theta=0)': 0, '<(1,0)> (theta=1)': 1, 'K^m (theta=2)': 2}[name]
        for r in (2, 5, 15, 60, 400):
            p1, p2 = phi(F, r), phi(F, 2 * r)
            rho = p2 / (2 * p1)
            appr = (KAPPA * mp.sqrt(2 * r - mp.mpf(1) / 6) - 1 + th) / (2 * (KAPPA * mp.sqrt(r - mp.mpf(1) / 6) - 1 + th))
            print('%-20s %6s %10s %10s %+12s %14s' % (name, r, mp.nstr(p1, 8), mp.nstr(rho, 7),
                  mp.nstr(rho - 1 / mp.sqrt(2), 5), mp.nstr(appr, 7)))
    print('\ncoefficient asymptotique (sqrt2-1)/(2 kappa) =', mp.nstr((mp.sqrt(2) - 1) / (2 * KAPPA), 6))
