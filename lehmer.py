"""Mahler measure tools for Lehmer's problem, in pure Python (no dependencies).

The Mahler measure of P(x) = a_n (x - r_1) ... (x - r_n) is M(P) = |a_n| * prod max(1, |r_i|).
Lehmer (1933) asked whether M(P) can be arbitrarily close to 1 without being 1.

Polynomials are lists of integer coefficients, LOWEST degree first:
x^3 - x - 1 is [-1, -1, 0, 1]. The command line takes coefficients HIGHEST degree first.

Commands:
    python lehmer.py measure 1 1 0 -1 -1 -1 -1 -1 0 1 1
                                  report on one polynomial (coefficients highest degree first):
                                  M from the roots, a certified interval for M, the exact
                                  Kronecker test, cyclotomic factors, root positions
    python lehmer.py search D [LIMIT]
                                  every monic reciprocal polynomial of degree 2D with
                                  coefficients in {-1, 0, 1} and 1 < M < LIMIT (default 1.3)
    python lehmer.py smyth N      check Smyth's bound M >= 1.3247... on every non-reciprocal
                                  monic polynomial of degree <= N with coefficients in {-1, 0, 1}
    python lehmer.py pierce N 1 -1 -1
                                  Lehmer's numbers Delta_k = prod (r_i^k - 1), k = 1..N, for a
                                  monic polynomial (coefficients highest degree first)

By MV2 of Munim, Inc. Made with AI agents working in MV2's name. MIT licence.
"""

import cmath
import itertools
import math
import sys
from functools import lru_cache

# Lehmer's polynomial x^10 + x^9 - x^7 - x^6 - x^5 - x^4 - x^3 + x + 1 and its measure
LEHMER = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1]
LEHMER_NUMBER = 1.1762808182599175
# x^3 - x - 1: its real root theta_0 is the smallest Pisot number and Smyth's constant
PLASTIC = [-1, -1, 0, 1]
THETA0 = 1.3247179572447460


# ---------------------------------------------------------------- integer polynomials

def trim(p):
    """Drop zero coefficients at the top."""
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def strip_x(p):
    """Return (k, q) with p = x^k q and q(0) != 0."""
    p = trim(p)
    if not p:
        raise ValueError("zero polynomial")
    k = 0
    while p[k] == 0:
        k += 1
    return k, p[k:]


def poly_mul(a, b):
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def divmod_monic(p, d):
    """Exact division of an integer polynomial by a monic one: (quotient, remainder)."""
    p, d = trim(p), trim(d)
    if not d or d[-1] != 1:
        raise ValueError("divisor must be monic")
    n = len(d) - 1
    r = p[:]
    if len(r) - 1 < n:
        return [], r
    q = [0] * (len(r) - n)
    for i in range(len(r) - 1, n - 1, -1):
        c = r[i]
        if c:
            q[i - n] = c
            for j in range(n + 1):
                r[i - n + j] -= c * d[j]
    return trim(q), trim(r)


def binomial(n, k):
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))


def is_reciprocal(p):
    """True if p(x) = +-x^n p(1/x) once factors of x are removed."""
    _, p = strip_x(p)
    r = p[::-1]
    return r == p or r == [-c for c in p]


def poly_str(p):
    p = trim(p)
    if not p:
        return "0"
    terms = []
    for i in range(len(p) - 1, -1, -1):
        c = p[i]
        if c == 0:
            continue
        a = abs(c)
        if i == 0:
            body = str(a)
        else:
            mono = "x" if i == 1 else "x^%d" % i
            body = mono if a == 1 else "%d%s" % (a, mono)
        terms.append(("-" if c < 0 else "+", body))
    out = ("-" if terms[0][0] == "-" else "") + terms[0][1]
    for sign, body in terms[1:]:
        out += " %s %s" % (sign, body)
    return out


# ---------------------------------------------------------------- cyclotomic polynomials

def totient(m):
    result, n, k = m, m, 2
    while k * k <= n:
        if n % k == 0:
            while n % k == 0:
                n //= k
            result -= result // k
        k += 1
    if n > 1:
        result -= result // n
    return result


@lru_cache(maxsize=None)
def cyclotomic(m):
    """Phi_m as a tuple of coefficients, from x^m - 1 = prod over d | m of Phi_d."""
    if m < 1:
        raise ValueError("m must be positive")
    num = [-1] + [0] * (m - 1) + [1]
    for d in range(1, m):
        if m % d == 0:
            num, r = divmod_monic(num, list(cyclotomic(d)))
            if r:
                raise ArithmeticError("Phi_%d does not divide x^%d - 1" % (d, m))
    return tuple(num)


def cyclotomic_part(p):
    """Split p = x^k * prod Phi_m^e_m * rest, where rest has no cyclotomic factor.

    Returns (k, {m: e_m}, rest). Exact. Only m with phi(m) <= deg p can occur, and
    phi(m) >= sqrt(m / 2), so m <= 2 deg(p)^2 + 2 covers every candidate.
    """
    k, p = strip_x(p)
    exps = {}
    n = len(p) - 1
    for m in range(1, 2 * n * n + 3):
        if totient(m) > len(p) - 1:
            continue
        phi = list(cyclotomic(m))
        while len(p) >= len(phi):
            q, r = divmod_monic(p, phi)
            if r:
                break
            p = q
            exps[m] = exps.get(m, 0) + 1
    return k, exps, p


# ---------------------------------------------------------------- Kronecker's theorem, exactly

def graeffe(p):
    """Root squaring: the integer polynomial whose roots are the squares of the roots of p.

    Writing p(x) = E(x^2) + x O(x^2), the result is (-1)^n (E(y)^2 - y O(y)^2); its leading
    coefficient is a_n^2, so M(graeffe(p)) = M(p)^2.
    """
    p = trim(p)
    n = len(p) - 1
    e2 = poly_mul(p[0::2], p[0::2])
    o2 = poly_mul(p[1::2], p[1::2])
    out = [0] * (n + 1)
    for i, c in enumerate(e2):
        out[i] += c
    for i, c in enumerate(o2):
        out[i + 1] -= c
    if n % 2:
        out = [-c for c in out]
    return out


def is_kronecker(p):
    """True iff M(p) = 1, i.e. p = +-x^k times a product of cyclotomic polynomials. Exact.

    This is Kronecker's proof run as an algorithm. If M(p) = 1, every root-squared
    polynomial p_k has |coefficient of x^j| <= C(n, j); the sequence reaches a fixed point
    (a product of Phi_m with m odd), and a fixed point forces every root onto the unit circle.
    If M(p) > 1 then M(p_k) = M(p)^(2^k) grows, and once it passes 2^n some coefficient must
    break its binomial bound.
    """
    _, p = strip_x(p)
    if abs(p[-1]) != 1 or abs(p[0]) != 1:
        return False
    if p[-1] < 0:
        p = [-c for c in p]
    n = len(p) - 1
    if n == 0:
        return True
    bound = [binomial(n, j) for j in range(n + 1)]
    q = graeffe(p)
    while True:
        if any(abs(c) > b for c, b in zip(q, bound)):
            return False
        nxt = graeffe(q)
        if nxt == q:
            return True
        q = nxt


def mahler_bounds(p, steps=12):
    """A certified interval [lo, hi] for M(p) from exact integer root squaring.

    For q = p after `steps` squarings (degree n), M(q) = M(p)^(2^steps), and
        max_j |q_j| / C(n, j) <= M(q) <= ||q||_2      (coefficient bound; Landau).
    Only the final logarithms and the 2^steps-th root are in floating point.
    """
    _, p = strip_x(p)
    n = len(p) - 1
    if n == 0:
        v = float(abs(p[0]))
        return v, v
    q = p
    for _ in range(steps):
        q = graeffe(q)
    scale = float(2 ** steps)
    lo = max(math.log(abs(c)) - math.log(binomial(n, j)) for j, c in enumerate(q) if c)
    hi = 0.5 * math.log(sum(c * c for c in q))
    return math.exp(lo / scale), math.exp(hi / scale)


# ---------------------------------------------------------------- numerical roots and M

def _horner(a, z):
    """Value and derivative at z of the polynomial with coefficients a (lowest first)."""
    v = 0j
    dv = 0j
    for c in reversed(a):
        dv = dv * z + v
        v = v * z + c
    return v, dv


def roots(p, max_iter=300, tol=1e-14):
    """All complex roots of p by the Aberth-Ehrlich iteration (double precision)."""
    p = trim(p)
    n = len(p) - 1
    if n < 1:
        return []
    a = [complex(c) / p[-1] for c in p]
    radius = max(abs(a[i]) ** (1.0 / (n - i)) for i in range(n))
    if radius == 0:
        return [0j] * n
    z = [radius * cmath.exp(2j * math.pi * k / n + 0.4j) for k in range(n)]
    for _ in range(max_iter):
        biggest = 0.0
        for i in range(n):
            v, dv = _horner(a, z[i])
            if v == 0:
                continue
            s = 0j
            for j in range(n):
                if j != i and z[i] != z[j]:
                    s += 1 / (z[i] - z[j])
            denom = dv - v * s
            if denom == 0:
                continue
            w = v / denom
            z[i] -= w
            biggest = max(biggest, abs(w) / max(1.0, abs(z[i])))
        if biggest < tol:
            break
    return z


def mahler_measure(p):
    """M(p). Cyclotomic factors and x are removed exactly first (they contribute 1); the
    rest is measured from its roots in double precision."""
    _, _, rest = cyclotomic_part(p)
    m = float(abs(rest[-1]))
    for r in roots(rest):
        if abs(r) > 1:
            m *= abs(r)
    return m


def jensen_log_measure(p, points=1024):
    """log M(p) as the average of log |p| over the unit circle (Jensen's formula), by the
    midpoint rule. Accurate when p has no roots on or very near the unit circle."""
    a = trim(p)
    total = 0.0
    for k in range(points):
        z = cmath.exp(2j * math.pi * (k + 0.5) / points)
        total += math.log(abs(_horner(a, z)[0]))
    return total / points


def root_profile(p, tol=1e-6):
    """Numbers of roots (outside, on, inside) the unit circle, numerically, ignoring x^k."""
    _, p = strip_x(p)
    out = on = inside = 0
    for r in roots(p):
        if abs(r) > 1 + tol:
            out += 1
        elif abs(r) < 1 - tol:
            inside += 1
        else:
            on += 1
    return out, on, inside


# ---------------------------------------------------------------- Lehmer's numbers Delta_k

def det_int(m):
    """Exact determinant of an integer matrix (Bareiss fraction-free elimination)."""
    a = [row[:] for row in m]
    n = len(a)
    if n == 0:
        return 1
    sign = 1
    prev = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            for r in range(k + 1, n):
                if a[r][k] != 0:
                    a[k], a[r] = a[r], a[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * a[k][k] - a[i][k] * a[k][j]) // prev
        prev = a[k][k]
    return sign * a[n - 1][n - 1]


def pierce(p, k):
    """Delta_k = prod (r_i^k - 1) over the roots of a monic p, exactly.

    prod over k-th roots of unity z of p(z) is the determinant of multiplication by p on
    Z[x]/(x^k - 1), a circulant matrix; and prod_i (r_i^k - 1) = (-1)^(n k) prod_z p(z).
    """
    p = trim(p)
    if not p or p[-1] != 1:
        raise ValueError("Lehmer's numbers need a monic polynomial")
    if k < 1:
        raise ValueError("k must be positive")
    n = len(p) - 1
    c = [[0] * k for _ in range(k)]
    for j in range(k):
        for i, coef in enumerate(p):
            c[(i + j) % k][j] += coef
    return (-1) ** (n * k) * det_int(c)


# ---------------------------------------------------------------- searches

def reciprocal_family(half_degree, coeffs=(-1, 0, 1)):
    """Monic reciprocal polynomials of degree 2*half_degree, middle coefficients from coeffs."""
    for tail in itertools.product(coeffs, repeat=half_degree):
        lower = [1] + list(tail)
        yield lower + lower[-2::-1]


def nonreciprocal_family(max_degree, coeffs=(-1, 0, 1)):
    """Monic non-reciprocal polynomials of degree <= max_degree with p(0) != 0."""
    for n in range(1, max_degree + 1):
        for low in itertools.product(coeffs, repeat=n):
            if low[0] == 0:
                continue
            p = list(low) + [1]
            if not is_reciprocal(p):
                yield p


def small_measures(polys, limit=1.3):
    """Sorted (M, p) for every p with 1 < M(p) < limit. M = 1 is decided exactly."""
    found = []
    for p in polys:
        if is_kronecker(p):
            continue
        m = mahler_measure(p)
        if m < limit:
            found.append((m, p))
    found.sort()
    return found


def smyth_check(max_degree):
    """(count, smallest M, a polynomial attaining it) over nonreciprocal_family(max_degree)."""
    count = 0
    best = None
    for p in nonreciprocal_family(max_degree):
        count += 1
        m = mahler_measure(p)
        if best is None or m < best[0]:
            best = (m, p)
    return count, best[0], best[1]


# ---------------------------------------------------------------- command line

def report(p):
    k, exps, rest = cyclotomic_part(p)
    lo, hi = mahler_bounds(p)
    out, on, inside = root_profile(p)
    factors = ["x^%d" % k] if k else []
    factors += ["Phi_%d^%d" % (m, e) if e > 1 else "Phi_%d" % m for m, e in sorted(exps.items())]
    print("polynomial:        %s" % poly_str(p))
    print("Mahler measure:    %.12f   (from the roots)" % mahler_measure(p))
    print("certified:         %.9f <= M <= %.9f   (12 exact root squarings)" % (lo, hi))
    print("M = 1 exactly:     %s   (Kronecker test)" % ("yes" if is_kronecker(p) else "no"))
    print("cyclotomic part:   %s" % (" ".join(factors) if factors else "none"))
    print("reciprocal:        %s" % ("yes" if is_reciprocal(p) else "no"))
    print("roots out/on/in:   %d / %d / %d   (unit circle, numerical)" % (out, on, inside))


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in ("-h", "--help", "help"):
        print(__doc__)
        return 0
    cmd, rest = args[0], args[1:]
    try:
        if cmd == "measure":
            if not rest:
                raise ValueError("give the coefficients, highest degree first")
            report([int(t) for t in rest][::-1])
            return 0
        if cmd == "search":
            half = int(rest[0]) if rest else 5
            limit = float(rest[1]) if len(rest) > 1 else 1.3
            found = small_measures(reciprocal_family(half), limit)
            print("%d polynomials of degree %d with 1 < M < %g:" % (len(found), 2 * half, limit))
            for m, p in found:
                print("%.12f  %s" % (m, poly_str(p)))
            return 0
        if cmd == "smyth":
            n = int(rest[0]) if rest else 6
            count, m, p = smyth_check(n)
            print("%d non-reciprocal polynomials checked; smallest M = %.12f for %s"
                  % (count, m, poly_str(p)))
            print("Smyth's bound theta_0 = %.12f: %s" % (THETA0, "holds" if m >= THETA0 - 1e-9 else "FAILS"))
            return 0
        if cmd == "pierce":
            if len(rest) < 2:
                raise ValueError("give N and the coefficients, highest degree first")
            n = int(rest[0])
            p = [int(t) for t in rest[1:]][::-1]
            for k in range(1, n + 1):
                print("%d  %d" % (k, pierce(p, k)))
            return 0
    except ValueError as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2
    print("unknown command %r; run with --help" % cmd, file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
