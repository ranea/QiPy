import unittest
import doctest

from hypothesis import given, assume, settings
from hypothesis.strategies import integers, tuples

import qipy.quadratic_integer  # load_tests
from qipy.quadratic_integer import is_square_free
from qipy.quadratic_integer import QuadraticIntegerRing

min_value_d = -5
max_value_d = 5
min_value_int = -10
max_value_int = 10
min_value_exp = 0
max_value_exp = 3
samples = 10
UFD = [
    -163, -67, -43, -19, -11, -7, -3, -2, -1, 1, 2, 3, 5, 6, 7, 11, 13, 14,
    17, 19, 21, 22, 23, 29, 31, 33, 37, 38, 41, 43, 46, 47, 53, 57, 59, 61,
    62, 67, 69, 71, 73, 77, 83, 86, 89, 93, 94, 97, 101, 103, 107, 109,
    113, 118, 127, 129, 131, 133, 134, 137, 139, 141, 149
]


class TestQuadraticIntegers(unittest.TestCase):
    @given(
        tuples(integers(min_value=min_value_int, max_value=max_value_int),
               integers(min_value=min_value_int, max_value=max_value_int)),
        tuples(integers(min_value=min_value_int, max_value=max_value_int),
               integers(min_value=min_value_int, max_value=max_value_int)),
        tuples(integers(min_value=min_value_int, max_value=max_value_int),
               integers(min_value=min_value_int, max_value=max_value_int)),
        integers(min_value=min_value_d, max_value=max_value_d)
    )
    @settings(max_examples=samples, deadline=None)
    def test_ring_propierties(self, l1, l2, l3, d):
        assume(l1)
        assume(l2)
        assume(l3)
        assume(d != 0)
        assume(d != 1)
        assume(is_square_free(d))

        O = QuadraticIntegerRing(d)

        alpha = O(coeff_e=l1)
        beta = O(coeff_e=l2)
        gamma = O(coeff_e=l3)
        zero = 0
        one = 1

        print("test_ring_propierties | {0} | {1} | {2} | {3}".format(O.__name__, alpha, beta, gamma))

        assert alpha + (beta + gamma) == (alpha + beta) + gamma
        assert alpha + beta == beta + alpha
        assert alpha + zero == alpha == zero + alpha
        assert alpha + (-alpha) == zero
        assert alpha * (beta * gamma) == (alpha * beta) * gamma
        assert alpha * one == alpha
        assert alpha * (beta + gamma) == (alpha * beta) + (alpha * gamma)
        assert (alpha + beta) * gamma == (alpha * gamma) + (beta * gamma)
        assert alpha * beta == beta * alpha

    @given(
        tuples(integers(min_value=min_value_int, max_value=max_value_int), integers(min_value=min_value_int, max_value=max_value_int)),
        tuples(integers(min_value=min_value_int, max_value=max_value_int), integers(min_value=min_value_int, max_value=max_value_int)),
        integers(min_value=min_value_d, max_value=max_value_d)
    )
    @settings(max_examples=samples, deadline=None)
    def test_composition_add_sub(self, l1, l2, d):
        assume(l1)
        assume(l2)
        assume(d != 0)
        assume(d != 1)
        assume(is_square_free(d))

        O = QuadraticIntegerRing(d)

        alpha = O(coeff_e=l1)
        beta = O(coeff_e=l2)

        print("test_composition_add_sub | {0} | {1} | {2}".format(O.__name__, alpha, beta))

        assert alpha == (alpha - beta) + beta

    @given(
        tuples(integers(min_value=min_value_int, max_value=max_value_int), integers(min_value=min_value_int, max_value=max_value_int)),
        tuples(integers(min_value=min_value_int, max_value=max_value_int), integers(min_value=min_value_int, max_value=max_value_int)),
        integers(min_value=min_value_d, max_value=max_value_d)
    )
    @settings(max_examples=samples)
    def test_composition_mul_div(self, l1, l2, d):
        assume(l1)
        assume(l2)
        assume(l2 != (0, 0))
        assume(d != 0)
        assume(d != 1)
        assume(is_square_free(d))

        O = QuadraticIntegerRing(d)

        alpha = O(coeff_e=l1)
        beta = O(coeff_e=l2)

        print("test_composition_mul_div | {0} | {1} | {2}".format(O.__name__, alpha, beta))

        assert not beta.is_unit() or (alpha == (alpha * beta) / beta)

    @given(
        tuples(integers(min_value=min_value_int, max_value=max_value_int), integers(min_value=min_value_int, max_value=max_value_int)),
        integers(min_value=min_value_int, max_value=max_value_int),
        integers(min_value=min_value_d, max_value=max_value_d)
    )
    @settings(max_examples=samples)
    def test_scalar_multiplication(self, l1, n, d):
        assume(l1)
        assume(n)
        assume(d != 0)
        assume(d != 1)
        assume(is_square_free(d))

        O = QuadraticIntegerRing(d)

        alpha = O(coeff_e=l1)

        print("test_scalar_multiplication | {0} | {1} | {2}".format(O.__name__, alpha, n))

        assert alpha * n == n * alpha

    @given(
        tuples(integers(min_value=min_value_int, max_value=max_value_int), integers(min_value=min_value_int, max_value=max_value_int)),
        integers(min_value=min_value_exp, max_value=max_value_exp),
        integers(min_value=min_value_exp, max_value=max_value_exp),
        integers(min_value=min_value_d, max_value=max_value_d)
    )
    @settings(max_examples=samples // 2, deadline=None)
    def test_exponentation_properties(self, l1, e, f, d):
        assume(l1)
        assume(e)
        assume(f)
        assume(d != 0)
        assume(d != 1)
        assume(is_square_free(d))

        O = QuadraticIntegerRing(d)

        alpha = O(coeff_e=l1)

        print("test_exponentation_properties | {0} | {1} | {2} | {3}".format(O.__name__, alpha, e, f))

        assert alpha ** e * alpha ** f == alpha ** (e + f)
        assert (alpha ** e) ** f == alpha ** (e * f)

    @given(
        integers(min_value=min_value_int, max_value=max_value_int),
        integers(min_value=min_value_d, max_value=max_value_d)
    )
    @settings(max_examples=samples // 2)
    def test_elements_with_norm(self, n, d):
        assume(n != 0)
        assume(d != 0)
        assume(d != 1)
        assume(is_square_free(d))

        O = QuadraticIntegerRing(d)

        elements = O.elements_with_norm(n)

        print("test_elements_with_norm | {0} | {1} | {2}".format(O.__name__, n, d))

        for e in elements:
            assert e.norm == n

    @given(
        tuples(integers(min_value=min_value_int, max_value=max_value_int), integers(min_value=min_value_int, max_value=max_value_int)),
        integers(min_value=min_value_d, max_value=max_value_d)
    )
    @settings(max_examples=samples // 2, deadline=None)
    def test_factor(self, l1, d):
        assume(l1)
        assume(d != 0)
        assume(d != 1)
        assume(is_square_free(d))
        assume(d in UFD)

        O = QuadraticIntegerRing(d)

        alpha = O(coeff_e=l1)

        assume(alpha != 0)
        assume(not alpha.is_unit())

        print("test_factor | {0} | {1} | {2}".format(O.__name__, alpha, d))

        factors = alpha.factor()

        product = O(1)
        for f in factors:
            product *= f

        assert alpha == product


def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(qipy.quadratic_integer))
    return tests


if __name__ == '__main__':
    unittest.main()
