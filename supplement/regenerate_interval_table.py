#!/usr/bin/env python3
"""Export Appendix A's table from mpmath.iv, using exact endpoint rounding.

The imported checker first verifies all three numerical conditions. This
table is an mpmath.iv cross-check, distinct from the Arb certificate.
No endpoint is converted to binary floating point for display decisions.
"""
import argparse
from contextlib import redirect_stdout
import csv
from fractions import Fraction
from io import StringIO
from pathlib import Path
import runpy


def endpoint(value):
    """Convert an mpmath binary endpoint to its exact rational value."""
    sign, mantissa, exponent, bits = value
    if bits < 0:
        raise ValueError("A finite interval endpoint is required")
    return (-1 if sign else 1) * Fraction(mantissa) * Fraction(2) ** exponent


def fixed(value, digits, mode):
    """Round a rational to a decimal string by integer arithmetic."""
    scaled = value * 10 ** digits
    lower, remainder = divmod(scaled.numerator, scaled.denominator)
    if mode == 'upper':
        number = lower + bool(remainder)
    elif mode == 'lower':
        number = lower
    elif mode == 'nearest':
        twice = 2 * remainder
        number = lower + (twice > scaled.denominator or
                          (twice == scaled.denominator and lower % 2 != 0))
    else:
        raise ValueError(mode)
    sign = '-' if number < 0 else ''
    whole, fraction = divmod(abs(number), 10 ** digits)
    return f'{sign}{whole}.{fraction:0{digits}d}'


def bounds(value):
    return tuple(endpoint(v) for v in value._mpi_)


def approximation(value):
    lo, hi = bounds(value)
    text = fixed(lo, 10, 'nearest')
    if text != fixed(hi, 10, 'nearest'):
        raise AssertionError('Interval does not determine a unique rounding')
    return text


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=root/'results'/'interval_table.csv')
    args = parser.parse_args()
    with redirect_stdout(StringIO()):
        data = runpy.run_path(str(Path(__file__).with_name('certify_mpmath_iv_crosscheck.py')))
    header = ['beta', 'A', 'Delta_A', 'C', 'U_upper']
    header += [f'{name}_{end}' for name in ['A', 'Delta_A', 'C', 'U'] for end in ['lower', 'upper_endpoint']]
    rows = []
    for i, values in enumerate(data['vals']):
        beta = f'{40+i:02d}'
        beta = beta[:-1] + '.' + beta[-1]
        row = [beta, *(approximation(v) for v in values)]
        u = data['ub'][i] if i < 10 else None
        row.append(fixed(bounds(u)[1], 4, 'upper') if u is not None else '')
        for v in (*values, u):
            row.extend([fixed(bounds(v)[0], 30, 'lower'), fixed(bounds(v)[1], 30, 'upper')]
                       if v is not None else ['', ''])
        rows.append(row)
        print('&'.join(row[:5]).rstrip('&') + ('&--' if u is None else '') + r'\\')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)
    print('==> INTERVAL TABLE: PASS (11 rows; exact endpoint rounding)')


if __name__ == '__main__':
    main()
