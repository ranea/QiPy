import unittest
import doctest

from hypothesis import given, assume, settings
from hypothesis.strategies import integers, tuples

import qipy.class_group  # load_tests
from qipy.class_group import IdealClass, ClassGroup

from qipy.quadratic_integer import is_square_free, QuadraticIntegerRing
from qipy.ideal import Ideal

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
class_groups = [
    [2, 1],
    [3, 1],
    [13, 1],
    [17, 1],
    [-1, 1],
    [-2, 1],
    [-3, 1],
    [-7, 1],
    [-14, 4],
    [79, 3],
    [-30, 4],
    [-65, 8],
    [-255, 12],
    [-299, 8],
]


class TestIdealClass(unittest.TestCase):
    @given(
        tuples(integers(min_value=min_value_int, max_value=max_value_int),
               integers(min_value=min_value_int, max_value=max_value_int)),
        tuples(integers(min_value=min_value_int, max_value=max_value_int),
               integers(min_value=min_value_int, max_value=max_value_int)),
        integers(min_value=min_value_d, max_value=max_value_d)
    )
    @settings(max_examples=samples, deadline=None)
    def test_inverse(self, l11, l12, d):
        assume(l11)
        assume(l11 != (0, 0))
        assume(l12)
        assume(l12 != (0, 0))
        assume(d != 0)
        assume(d != 1)
        assume(is_square_free(d))

        O = QuadraticIntegerRing(d)
        a = IdealClass(Ideal(O(coeff_e=l11), O(coeff_e=l12)))

        print("test_inverse | {0} | {1}".format(O.__name__, a))

        inverse = a.inverse()

        assert (a * inverse).is_trivial()

    @given(
        tuples(integers(min_value=min_value_int, max_value=max_value_int),
               integers(min_value=min_value_int, max_value=max_value_int)),
        tuples(integers(min_value=min_value_int, max_value=max_value_int),
               integers(min_value=min_value_int, max_value=max_value_int)),
        integers(min_value=min_value_d, max_value=max_value_d)
    )
    @settings(max_examples=samples, deadline=None)
    def test_order(self, l11, l12, d):
        assume(l11)
        assume(l11 != (0, 0))
        assume(l12)
        assume(l12 != (0, 0))
        assume(d != 0)
        assume(d != 1)
        assume(is_square_free(d))

        O = QuadraticIntegerRing(d)
        a = IdealClass(Ideal(O(coeff_e=l11), O(coeff_e=l12)))

        print("test_order | {0} | {1}".format(O.__name__, a))

        order = a.order
        I = a.representative
        product = Ideal.unit_ideal(d)
        for i in range(order):
            product *= I
        assert product.is_principal()

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
    @settings(max_examples=samples, deadline=None)
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

        a = IdealClass(Ideal(O(coeff_e=l11), O(coeff_e=l12)))
        b = IdealClass(Ideal(O(coeff_e=l21), O(coeff_e=l22)))
        trivial_class = IdealClass.trivial_class(d)

        print("test_multiplication | {} | {} | {}".format(O.__name__, a, b))

        assert a * b == b * a
        assert a == a * trivial_class


class TestClassGroup(unittest.TestCase):
    def test_class_group(self):
        for d, h in class_groups:
            print("test_class_group | {}".format(d))
            G = ClassGroup(d)
            assert G.class_number == h


def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(qipy.class_group))
    return tests


if __name__ == '__main__':
    unittest.main()
