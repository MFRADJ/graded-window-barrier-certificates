#!/usr/bin/env python3
"""V2-S4 -- regime h = r*m a densite fixee, m grand (test de l'uniformite en m).
Codes a caractere ferme : S={0} -> A^m ; S=<(1,0)_i> (k=m, mots impairs) -> (A+B)^m ;
sous-code pair de ce dernier (k=m-1) -> ((A+B)^m+(A-B)^m)/2 ;
hors hypothese : S=K^m (k=2m) -> (A+3B)^m.
"""
import math, sys
from flint import fmpz_poly
import importlib.util, os
_here = os.path.dirname(os.path.abspath(__file__))
_src = next(os.path.join(_here, f) for f in ('audit_exact_V2-S4.py', 'audit_exact_v2s3.py')
            if os.path.exists(os.path.join(_here, f)))
_spec = importlib.util.spec_from_file_location('audit_exact', _src)
_mod = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_mod)
series, mul, powtrunc, coeffs, lg, KAPPA, ISQ2 = (_mod.series, _mod.mul, _mod.powtrunc,
                                                   _mod.coeffs, _mod.lg, _mod.KAPPA, _mod.ISQ2)


def rho_at(d, h2):
    import itertools
    C = list(itertools.accumulate(d))
    def win(lo, hi):
        lo = max(lo, 0)
        return 0 if hi < lo else C[hi] - (C[lo - 1] if lo > 0 else 0)
    lam = d[h2].bit_length() - 1          # floor(log2 d) exact
    best = max(lg(max(1, win(t, t + 2 * h2))) + lg(max(1, win(-t, 2 * h2 - t)))
               for t in range(-2 * h2, 2 * h2 + 1))
    return best / (4 * lam), lam, lg(C[h2]), lg(C[2 * h2])


def run(r, ms, hs_half=False):
    print('\n== h = %s m %s' % (r, '(+1/2)' if hs_half else ''))
    print('%4s %-10s %3s %8s %9s %9s %9s %9s %8s' % ('m', 'code', 'k', 'h', 'rho', 'borne(i)', 'borne(ii)', 'rho-1/r2', 'marge(2)'))
    for m in ms:
        h2 = int(2 * r * m) + (1 if hs_half else 0)
        h = h2 / 2
        L = 4 * h2 + 2
        A, B = series(L)
        codes = []
        Am = powtrunc(A, m, L)
        ApB = powtrunc(A + B, m, L)
        AmB = powtrunc(A - B, m, L)
        codes.append(('{0}', 0, coeffs(Am, L)))
        codes.append(('<(1,0)>', m, coeffs(ApB, L)))
        ev = [(x + y) // 2 for x, y in zip(coeffs(ApB, L), coeffs(AmB, L))]
        codes.append(('pair', m - 1, ev))
        codes.append(('K^m(k=2m)', 2 * m, coeffs(powtrunc(A + 3 * B, m, L), L)))
        X = KAPPA * math.sqrt(m * h)
        e1 = 4 * m / 3 + 24 * (m * h) ** 0.25 + 1; e2 = e1 + 1
        # marge de (2) pour le secteur 0 (bits), a l'energie h
        import itertools
        N0 = sum(coeffs(Am, L)[:h2 + 1])
        marg2 = lg(N0) - (X - 4 * m / 3 - 20 * (m * h) ** 0.25 - 1)
        for name, k, d in codes:
            if d[h2] == 0:
                print('%4d %-10s %3d %8.1f   (niveau non admissible)' % (m, name, k, h)); continue
            rho, lam, Dh, D2h = rho_at(d, h2)
            bi = ISQ2 * (X - e1 / math.sqrt(2)) / (X + k)
            bii = ISQ2 * (X + k / math.sqrt(2)) / (X - e2) if X > e2 else float('inf')
            print('%4d %-10s %3d %8.1f %9.5f %9.5f %9.4f %+9.5f %8.1f' % (m, name, k, h, rho, bi, bii, rho - ISQ2, marg2))
        sys.stdout.flush()


if __name__ == '__main__':
    run(15, [1, 2, 4, 8, 16, 24, 32, 48])
    run(15, [1, 3, 9, 27], hs_half=True)
    run(60, [1, 2, 4, 8, 12])
    run(400, [1, 2, 3])
