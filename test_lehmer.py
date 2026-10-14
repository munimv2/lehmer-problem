import itertools
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import lehmer as L  # noqa: E402

GOLDEN = (1 + 5 ** 0.5) / 2


def lucas(n):
    a, b = 2, 1
    for _ in range(n):
        a, b = b, a + b
    return a


class TestBasics(unittest.TestCase):
    def test_lehmer_polynomial(self):
        self.assertEqual(L.poly_str(L.LEHMER), "x^10 + x^9 - x^7 - x^6 - x^5 - x^4 - x^3 + x + 1")
        self.assertEqual(L.poly_str(L.PLASTIC), "x^3 - x - 1")
        self.assertEqual(L.poly_str([3, 0, -2]), "-2x^2 + 3")

    def test_divmod_monic(self):
        p = L.poly_mul([1, 2, 3], [-1, 0, 1])
        q, r = L.divmod_monic(p, [-1, 0, 1])
        self.assertEqual(q, [1, 2, 3])
        self.assertEqual(r, [])
        q, r = L.divmod_monic([5, 0, 1], [1, 1])
        self.assertEqual(L.poly_mul(q, [1, 1])[0] + r[0], 5)
        with self.assertRaises(ValueError):
            L.divmod_monic([1, 2], [1, 2])

    def test_reciprocal(self):
        self.assertTrue(L.is_reciprocal(L.LEHMER))
        self.assertTrue(L.is_reciprocal([-1, 1]))           # x - 1 is anti-reciprocal
        self.assertTrue(L.is_reciprocal([0, 0, 1, 1]))      # x^2 (x + 1)
        self.assertFalse(L.is_reciprocal(L.PLASTIC))


class TestCyclotomic(unittest.TestCase):
    def test_small(self):
        self.assertEqual(L.cyclotomic(1), (-1, 1))
        self.assertEqual(L.cyclotomic(2), (1, 1))
        self.assertEqual(L.cyclotomic(3), (1, 1, 1))
        self.assertEqual(L.cyclotomic(4), (1, 0, 1))
        self.assertEqual(L.cyclotomic(6), (1, -1, 1))
        self.assertEqual(L.cyclotomic(12), (1, 0, -1, 0, 1))

    def test_degree_is_totient(self):
        for m in range(1, 61):
            self.assertEqual(len(L.cyclotomic(m)) - 1, L.totient(m), m)

    def test_phi_105(self):
        # 105 = 3 * 5 * 7 is the first index with a coefficient outside {-1, 0, 1}
        for m in range(1, 105):
            self.assertTrue(set(L.cyclotomic(m)) <= {-1, 0, 1}, m)
        self.assertIn(-2, L.cyclotomic(105))

    def test_cyclotomic_part(self):
        phi3, phi7 = list(L.cyclotomic(3)), list(L.cyclotomic(7))
        p = L.poly_mul([0, 0, 1], L.poly_mul(L.LEHMER, L.poly_mul(phi7, L.poly_mul(phi3, phi3))))
        k, exps, rest = L.cyclotomic_part(p)
        self.assertEqual(k, 2)
        self.assertEqual(exps, {3: 2, 7: 1})
        self.assertEqual(rest, L.LEHMER)

    def test_lehmer_has_no_cyclotomic_factor(self):
        self.assertEqual(L.cyclotomic_part(L.LEHMER), (0, {}, L.LEHMER))


class TestKronecker(unittest.TestCase):
    def test_graeffe(self):
        self.assertEqual(L.graeffe([-2, 1]), [-4, 1])
        self.assertEqual(L.graeffe([1, 0, 1]), [1, 2, 1])            # Phi_4 -> Phi_2^2
        self.assertEqual(L.graeffe([1, 1]), [-1, 1])                 # Phi_2 -> Phi_1
        self.assertEqual(L.graeffe(list(L.cyclotomic(15))), list(L.cyclotomic(15)))
        self.assertEqual(L.graeffe([1, -1, 1]), [1, 1, 1])           # Phi_6 -> Phi_3

    def test_cyclotomic_products(self):
        self.assertTrue(L.is_kronecker([1]))
        self.assertTrue(L.is_kronecker([-1]))
        self.assertTrue(L.is_kronecker([1, -1]))
        self.assertTrue(L.is_kronecker([0, 0, 0] + list(L.cyclotomic(5))))
        p = L.poly_mul(list(L.cyclotomic(8)), L.poly_mul(list(L.cyclotomic(9)), list(L.cyclotomic(9))))
        self.assertTrue(L.is_kronecker(p))
        self.assertTrue(L.is_kronecker(list(L.cyclotomic(105))))
        self.assertTrue(L.is_kronecker(list(L.cyclotomic(64))))

    def test_not_kronecker(self):
        self.assertFalse(L.is_kronecker(L.LEHMER))
        self.assertFalse(L.is_kronecker(L.PLASTIC))
        self.assertFalse(L.is_kronecker([1, 2]))
        self.assertFalse(L.is_kronecker([2, 0, 1]))
        self.assertFalse(L.is_kronecker([-1, -1, 1]))

    def test_agrees_with_factorisation(self):
        # every integer polynomial of degree <= 4 with coefficients in {-1, 0, 1}
        for n in range(1, 5):
            for low in itertools.product((-1, 0, 1), repeat=n):
                for lead in (-1, 1):
                    p = list(low) + [lead]
                    rest = L.cyclotomic_part(p)[2]
                    self.assertEqual(L.is_kronecker(p), rest in ([1], [-1]), p)


class TestMeasure(unittest.TestCase):
    def test_exact_values(self):
        self.assertAlmostEqual(L.mahler_measure([-2, 1]), 2.0, places=12)
        self.assertAlmostEqual(L.mahler_measure([-1, 2]), 2.0, places=12)
        self.assertAlmostEqual(L.mahler_measure([-1, 0, 3]), 3.0, places=12)
        self.assertAlmostEqual(L.mahler_measure([-1, -1, 1]), GOLDEN, places=12)
        self.assertAlmostEqual(L.mahler_measure([1, -3, 1]), GOLDEN ** 2, places=12)
        self.assertAlmostEqual(L.mahler_measure(L.PLASTIC), L.THETA0, places=12)
        self.assertAlmostEqual(L.mahler_measure(L.LEHMER), L.LEHMER_NUMBER, places=12)
        self.assertEqual(L.mahler_measure(list(L.cyclotomic(30))), 1.0)

    def test_theta0_is_root(self):
        t = L.THETA0
        self.assertAlmostEqual(t ** 3 - t - 1, 0.0, places=13)
        x = L.LEHMER_NUMBER
        self.assertAlmostEqual(sum(c * x ** i for i, c in enumerate(L.LEHMER)), 0.0, places=12)

    def test_multiplicative(self):
        p = L.poly_mul(L.LEHMER, L.poly_mul(L.PLASTIC, [0, -2, 1]))
        self.assertAlmostEqual(L.mahler_measure(p), 2 * L.LEHMER_NUMBER * L.THETA0, places=10)

    def test_jensen(self):
        # (2x - 1)(x - 3) = 2x^2 - 7x + 3 has M = 6; no test polynomial has roots near |z| = 1
        for p in ([-1, -1, 1], L.PLASTIC, [-2, 1], [3, -7, 2]):
            self.assertAlmostEqual(L.jensen_log_measure(p), math.log(L.mahler_measure(p)), places=10)
        self.assertAlmostEqual(L.mahler_measure([3, -7, 2]), 6.0, places=12)

    def test_certified_bounds(self):
        eps = 1e-12  # the final logarithm and root are in floating point
        lo, hi = L.mahler_bounds(L.LEHMER, steps=16)
        self.assertLessEqual(lo, L.LEHMER_NUMBER + eps)
        self.assertGreaterEqual(hi, L.LEHMER_NUMBER - eps)
        self.assertLess(hi - lo, 1e-3)
        self.assertGreater(lo, 1.17)
        lo, hi = L.mahler_bounds([-2, 1], steps=10)
        self.assertAlmostEqual(lo, 2.0, places=9)
        self.assertAlmostEqual(hi, 2.0, places=9)
        lo, hi = L.mahler_bounds(L.PLASTIC, steps=14)
        self.assertLessEqual(lo, L.THETA0 + eps)
        self.assertGreaterEqual(hi, L.THETA0 - eps)
        self.assertGreater(lo, 1.32)

    def test_bounds_of_cyclotomic_product(self):
        lo, hi = L.mahler_bounds(L.poly_mul([0, 1], list(L.cyclotomic(12))), steps=12)
        self.assertGreaterEqual(lo, 1.0)
        self.assertLess(hi, 1.01)

    def test_salem_profile(self):
        # Lehmer's number is a Salem number: one root outside, its inverse inside, 8 on the circle
        self.assertEqual(L.root_profile(L.LEHMER), (1, 8, 1))
        self.assertEqual(L.root_profile(L.PLASTIC), (1, 0, 2))


class TestPierce(unittest.TestCase):
    def test_mersenne(self):
        for k in range(1, 13):
            self.assertEqual(L.pierce([-2, 1], k), 2 ** k - 1)

    def test_lucas(self):
        # x^2 - x - 1: Delta_k = (alpha beta)^k - L_k + 1 = (-1)^k - L_k + 1
        for k in range(1, 16):
            self.assertEqual(L.pierce([-1, -1, 1], k), (-1) ** k - lucas(k) + 1)

    def test_growth_rate(self):
        d = L.pierce([-1, -1, 1], 40)
        self.assertAlmostEqual(math.log(abs(d)) / 40, math.log(GOLDEN), places=6)

    def test_cyclotomic_vanishes(self):
        self.assertEqual(L.pierce([1, 1, 1], 3), 0)

    def test_needs_monic(self):
        with self.assertRaises(ValueError):
            L.pierce([1, 2], 3)


class TestSearches(unittest.TestCase):
    def test_lehmer_is_smallest_in_degree_10_family(self):
        found = L.small_measures(L.reciprocal_family(5), limit=1.3)
        self.assertTrue(found)
        self.assertAlmostEqual(found[0][0], L.LEHMER_NUMBER, places=10)
        self.assertIn(L.LEHMER, [p for m, p in found if abs(m - L.LEHMER_NUMBER) < 1e-9])
        for m, p in found:
            self.assertGreater(m, L.LEHMER_NUMBER - 1e-9)
            self.assertLess(m, 1.3)

    def test_smyth_bound(self):
        count, m, p = L.smyth_check(6)
        self.assertGreater(count, 300)
        self.assertAlmostEqual(m, L.THETA0, places=10)
        self.assertFalse(L.is_reciprocal(p))
        for q in L.nonreciprocal_family(5):
            self.assertGreater(L.mahler_measure(q), L.THETA0 - 1e-9, q)


class TestCli(unittest.TestCase):
    def test_commands(self):
        self.assertEqual(L.main(["measure", "1", "1", "0", "-1", "-1", "-1", "-1", "-1", "0", "1", "1"]), 0)
        self.assertEqual(L.main(["pierce", "5", "1", "-1", "-1"]), 0)
        self.assertEqual(L.main(["measure"]), 2)
        self.assertEqual(L.main(["nonsense"]), 2)


if __name__ == "__main__":
    unittest.main()
