import unittest
import doctest

from sympy.ntheory.generate import prime

from hypothesis import given, assume, settings
from hypothesis.strategies import integers, tuples

import qipy.ideal  # load_tests
from qipy.ideal import Ideal

from qipy.quadratic_integer import is_square_free, QuadraticIntegerRing

min_value_d = -5
max_value_d = 5
min_value_int = -10
max_value_int = 10
min_value_exp = 0
max_value_exp = 3
samples = 10
max_index_prime = 10
UFD = [
    -163, -67, -43, -19, -11, -7, -3, -2, -1, 1, 2, 3, 5, 6, 7, 11, 13, 14,
    17, 19, 21, 22, 23, 29, 31, 33, 37, 38, 41, 43, 46, 47, 53, 57, 59, 61,
    62, 67, 69, 71, 73, 77, 83, 86, 89, 93, 94, 97, 101, 103, 107, 109,
    113, 118, 127, 129, 131, 133, 134, 137, 139, 141, 149
]


class TestIdeal(unittest.TestCase):
    @given(
        tuples(integers(min_value=min_value_int, max_value=max_value_int),
               integers(min_value=min_value_int, max_value=max_value_int)),
        tuples(integers(min_value=min_value_int, max_value=max_value_int),
               integers(min_value=min_value_int, max_value=max_value_int)),
        integers(min_value=min_value_d, max_value=max_value_d)
    )
    @settings(max_examples=samples)
    def test_lattice_reduce(self, l11, l12, d):
        assume(l11)
        assume(l11 != (0, 0))
        assume(l12)
        assume(l12 != (0, 0))
        assume(d != 0)
        assume(d != 1)
        assume(is_square_free(d))

        O = QuadraticIntegerRing(d)
        I = Ideal(O(coeff_e=l11), O(coeff_e=l12))

        print("test_lattice_reduce | {} | {}".format(O.__name__, I))

        assert I.contain(O(coeff_e=l11)) and I.contain(O(coeff_e=l12))

    @given(
        tuples(integers(min_value=min_value_int, max_value=max_value_int),
               integers(min_value=min_value_int, max_value=max_value_int)),
        tuples(integers(min_value=min_value_int, max_value=max_value_int),
               integers(min_value=min_value_int, max_value=max_value_int)),
        integers(min_value=min_value_d, max_value=max_value_d),
        integers(min_value=1, max_value=max_index_prime)
    )
    @settings(max_examples=samples)
    def test_prime_divisors(self, l11, l12, d, index_p):
        assume(l11)
        assume(l11 != (0, 0))
        assume(l12)
        assume(l12 != (0, 0))
        assume(d != 0)
        assume(d != 1)
        assume(is_square_free(d))

        O = QuadraticIntegerRing(d)
        p = prime(index_p)
        divisors = Ideal.prime_divisors(p, d)
        pI = Ideal(O(p))

        print("test_prime_divisors | {} | {} | {}".format(O.__name__, p, divisors))

        for I in divisors:
            assert I.is_prime()
            assert I.divide(pI)

    @given(
        tuples(integers(min_value=min_value_int, max_value=max_value_int),
               integers(min_value=min_value_int, max_value=max_value_int)),
        tuples(integers(min_value=min_value_int, max_value=max_value_int),
               integers(min_value=min_value_int, max_value=max_value_int)),
        tuples(integers(min_value=min_value_int, max_value=max_value_int),
               integers(min_value=min_value_int, max_value=max_value_int)),
        tuples(integers(min_value=min_value_int, max_value=max_value_int),
               integers(min_value=min_value_int, max_value=max_value_int)),
        integers(min_value=min_value_d, max_value=max_value_d)
    )
    @settings(max_examples=samples)
    def test_multiplication(self, l11, l12, l21, l22, d):
        assume(l11)
        assume(l11 != (0, 0))
        assume(l12)
        assume(l12 != (0, 0))
        assume(l21)
        assume(l21 != (0, 0))
        assume(l22)
        assume(l22 != (0, 0))
        assume(d != 0)
        assume(d != 1)
        assume(is_square_free(d))

        O = QuadraticIntegerRing(d)

        I = Ideal(O(coeff_e=l11), O(coeff_e=l12))
        J = Ideal(O(coeff_e=l21), O(coeff_e=l22))
        unit = Ideal.unit_ideal(d)

        print("test_multiplication | {} | {} | {}".format(O.__name__, I, J))

        assert I * J == J * I
        assert I == I * unit

    @given(
        tuples(integers(min_value=min_value_int, max_value=max_value_int),
               integers(min_value=min_value_int, max_value=max_value_int)),
        tuples(integers(min_value=min_value_int, max_value=max_value_int),
               integers(min_value=min_value_int, max_value=max_value_int)),
        integers(min_value=min_value_d, max_value=max_value_d)
    )
    @settings(max_examples=samples)
    def test_factor(self, l11, l12, d):
        assume(l11)
        assume(l11 != (0, 0))
        assume(l12)
        assume(l12 != (0, 0))
        assume(d != 0)
        assume(d != 1)
        assume(is_square_free(d))

        O = QuadraticIntegerRing(d)
        I = Ideal(O(coeff_e=l11), O(coeff_e=l12))
        assume(I.is_proper())

        print("test_factor | {0} | {1}".format(O.__name__, I))

        factors = I.factor()

        product = Ideal.unit_ideal(O.d)
        for f in factors:
            product *= f

        assert I == product


def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(qipy.ideal))
    return tests


if __name__ == '__main__':
    unittest.main()
