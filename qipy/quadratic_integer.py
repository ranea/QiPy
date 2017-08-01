"""This module allows to work with `quadratic integers`_.

To compute operations with quadratic integers within the same quadratic ring,
all operands must be created as :any:`QuadraticInteger` objects. To do so:

 1. Create the quadratic integer ring :math:`\mathcal{O}_{\mathbb{Q}[\sqrt{d}]}`
    with the function :any:`QuadraticIntegerRing`.

    >>> Zi = QuadraticIntegerRing(-1)

 2. Use the returned factory to create quadratic integers of the form
    :math:`a + b \sqrt{d}`.

    >>> alpha = Zi("1 + I")
    >>> beta = Zi("1 - I")

 3. Then, compute the desired operations with the available operators and methods.

    >>> alpha / beta
    I
    >>> (alpha / beta).is_unit()
    True

Note that this module, ``quadratic_integer``, need to be imported to use
its classes and functions. There are basically two ways to import it:

 1. Import all functions and classes of QiPy: ::

    >>> from qipy import *
    >>> Zi = QuadraticIntegerRing(-1)

 2. Import only the package, and use the package's identifier to
    access the classes and functions: ::

    >>> import qipy
    >>> Zi = qipy.QuadraticIntegerRing(-1)

"""
from functools import lru_cache

from sympy import sqrt, simplify, Rational, roots, Abs, isprime, factorint
from sympy import denom, expand
from sympy.polys.numberfields import minimal_polynomial
from sympy.solvers.diophantine import diop_DN

from qipy.utilities import is_square_free


@lru_cache()
def QuadraticIntegerRing(d):
    """Return the class of quadratic integers of the form
    :math:`a + b \sqrt{d}`.

    The returned class is used as a factory, that is, to create quadratic
    integer.

        >>> Zi = QuadraticIntegerRing(-1)
        >>> Zi(sqrt(-1))
        I

    Args:
        d: a square-free integer.

    Return:
        The factory of quadratic integers.
    """
    class QuadraticInteger(object):
        """Represent a quadratic integer.

            >>> Zi = QuadraticIntegerRing(-1)
            >>> alpha = Zi(1 - sqrt(-1))
            >>> beta = Zi(sqrt(-1))
            >>> alpha
            1 - I
            >>> beta
            I
            >>> alpha + beta
            1
            >>> alpha * beta
            1 + I
            >>> beta ** (2)
            -1

        This class supports the operators  ``+``, ``-``, ``*``, ``/`` y ``**``
        with their natural meaning.

        There are four ways to create a quadratic integer:

        - With a python expression.

            >>> O = QuadraticIntegerRing(-3)
            >>> alpha = O(Rational(1,2) * (-1 + sqrt(-3)))
            >>> alpha
            -1/2 + sqrt(3)*I/2

        - With a SymPy expression surrounded by ``"``.

            >>> O = QuadraticIntegerRing(-3)
            >>> alpha = O("1/2 * (-1 + I * sqrt(3))")
            >>> alpha
            -1/2 + sqrt(3)*I/2

        - With a pair of integers that represents the coefficients respect
          the `integral basis`_ :math:`\{1, e\}` (see below).

            >>> O = QuadraticIntegerRing(-3)
            >>> alpha = O(coeff_e=(-1, 1))
            >>> alpha
            -1/2 + sqrt(3)*I/2

        - With a pair of rationals that represents the coefficients respect
          the basis :math:`\{1, \sqrt{d}\}` (see below).

            >>> O = QuadraticIntegerRing(-3)
            >>> alpha = O(coeff_d=(Rational(-1, 2), Rational(1, 2)))
            >>> alpha
            -1/2 + sqrt(3)*I/2

        Warning:
            If the input contains fractions, each fraction must be written
            as a Rational object or the input must be surround by ``"``.

            More info at `SymPy's documentation`_.

        Note:
            The element :math:`e` of the `integral basis`_ is defined as follows:

            .. math::

                \sqrt{d}, \quad \mathrm{if} \ d \ \\ne 1 \mod 4

                (1 + \sqrt{d})/2, \quad \mathrm{if} \ d = 1 \mod 4

            Notice that in the first case, the basis :math:`\{1, e\}` and
            :math:`\{1, \sqrt{d}\}` are the same and so the coefficients.

        Attributes:
            d: the non-square free integer that defines the ring of quadratic
                integers (class attribute).
            e: the element of the integral basis (class attribute, see below).
            norm: the `norm`_ of the quadratic integer.
            trace: the `trace`_ of the quadratic integer.
            conjugate: the `algebraic conjugate`_ of the quadratic integer.

        Note:
            Given :math:`a + b \sqrt{d}`, its conjugate is :math:`a - b \sqrt{d}`,
            which is not the complex conjugate when :math:`d > 0`.

                >>> Zi = QuadraticIntegerRing(-1)
                >>> alpha = Zi(1 + sqrt(-1))
                >>> alpha
                1 + I
                >>> alpha.conjugate
                1 - I
                >>> O = QuadraticIntegerRing(5)
                >>> beta = O(1 + sqrt(5))
                >>> beta
                1 + sqrt(5)
                >>> beta.conjugate
                -sqrt(5) + 1

        .. _integral basis: https://en.wikipedia.org/wiki/Algebraic_number_field#Integral_basis
        .. _SymPy's documentation: http://docs.sympy.org/dev/gotchas.html#python-numbers-vs-sympy-numbers
        .. _norm: https://en.wikipedia.org/wiki/Field_norm
        .. _trace: https://en.wikipedia.org/wiki/Field_trace
        .. _algebraic conjugate: https://en.wikipedia.org/wiki/Conjugate_element_(field_theory)
        """
        d = None
        e = None

        def __init__(self, value=None, coeff_e=None, coeff_d=None):
            if value is not None:
                self.value = simplify(value)
            elif coeff_e is not None:
                self.value = simplify(coeff_e[0] + coeff_e[1] * QuadraticInteger.e)
                self._coeff_e = coeff_e
            elif coeff_d is not None:
                self.value = simplify(coeff_d[0] + coeff_d[1] * sqrt(QuadraticInteger.d))
                self._coeff_d = coeff_d
            else:
                raise ValueError("Invalid input.")

            if QuadraticInteger.d < 0:
                self.conjugate = self.value.conjugate()
            else:
                polynml = minimal_polynomial(self.value, "x")
                poly_roots = list(roots(polynml))  # length 1 or 2
                if simplify(self.value - poly_roots[0]) == 0:
                    self.conjugate = poly_roots[-1]
                else:
                    self.conjugate = poly_roots[0]

            self.norm = simplify(self.value * self.conjugate)
            self.trace = simplify(self.value + self.conjugate)

            if not(self.norm.is_integer and self.trace.is_integer):
                raise ValueError("It isn't a quadratic integer.")
            else:
                for v in self.coeff_d:
                    sv = simplify(v)
                    if not sv.is_rational or denom(sv) not in [1, 2]:
                        raise ValueError("It isn't a quadratic integer.")

        # ----- properties ------
        # properties are used to enable cache for attributes calculated

        @property
        def coeff_e(self):
            """Tuple: the integer coefficients with respect to the
            `integral basis`_ :math:`\{1, e\}`.

                >>> O = QuadraticIntegerRing(-3)
                >>> alpha = O(Rational(1,2) * (-1 + sqrt(-3)))
                >>> alpha
                -1/2 + sqrt(3)*I/2
                >>> alpha.coeff_e
                (-1, 1)

            .. _integral basis: https://en.wikipedia.org/wiki/Algebraic_number_field#Integral_basis
            """
            try:
                return self._coeff_e
            except:
                if d % 4 != 1:
                    self._coeff_e = self.coeff_d
                    return self.coeff_e
                else:
                    x, y = self.coeff_d
                    a = x - y
                    b = 2 * y
                    self._coeff_e = (a, b)
                    return self._coeff_e

        @property
        def coeff_d(self):
            """Tuple: the rational (`half-integer`_) coefficients with respect
            to the basis :math:`\{1, \sqrt{d}\}`.

                >>> O = QuadraticIntegerRing(-3)
                >>> alpha = O(Rational(1,2) * (-1 + sqrt(-3)))
                >>> alpha
                -1/2 + sqrt(3)*I/2
                >>> alpha.coeff_d
                (-1/2, 1/2)

            .. _half-integer: https://en.wikipedia.org/wiki/half-integer
            """
            try:
                return self._coeff_d
            except:
                value = expand(self.value)
                x = value.coeff(sqrt(QuadraticInteger.d), 0)
                y = value.coeff(sqrt(QuadraticInteger.d), 1)
                self._coeff_d = (x, y)
                return self._coeff_d

        # ----- class methods ------

        @classmethod
        def elements_with_norm(cls, n):
            """Return the quadratic integers of norm :math:`n`.

            If :math:`d` is negative, it returns all of them. Otherwise,
            it returns a list of generators (see *Structure of solutions to*
            :math:`x^2 − Dy^2 = N` in `Solving the general Pell equation`_
            by J.P. Robertson).

                >>> Zi = QuadraticIntegerRing(-1)
                >>> Zi.elements_with_norm(1)
                [1, -1, I, -I]
                >>> O = QuadraticIntegerRing(5)
                >>> O.elements_with_norm(1)
                [sqrt(5)/2 + 3/2, 3*sqrt(5)/2 + 7/2, 4*sqrt(5) + 9]

            Returns:
                List[QuadraticInteger]: the quadratic integers of norm :math:`n`.

            .. _Solving the general Pell equation: http://www.jpr2718.org/pell.pdf
            """
            def diop_DN_negative_d(d, n):
                """Solve the Pell equation for negatives values of d."""
                if n < 0:
                    return []
                if n == 0:
                    return [(0, 0)]

                solutions = []
                for y in range(0, int(sqrt(Rational(n, -d)) + 1)):
                    x = sqrt(n + d * y**2)
                    if x.is_integer:
                        # to evade duplicate solutions
                        if x == 0:
                            # y != 0
                            partial_solutions = [(0, y), (0, -y)]
                        elif y == 0:
                            # x != 0
                            partial_solutions = [(x, 0), (-x, 0)]
                        else:
                            # x != 0, y != 0
                            partial_solutions = [(x, y), (x, -y), (-x, y), (-x, -y)]
                        solutions += partial_solutions
                return solutions

            d = QuadraticInteger.d

            if d % 4 != 1:
                if d < 0:
                    sols = diop_DN_negative_d(d, n)
                else:
                    sols = diop_DN(d, n)
                return [QuadraticInteger(coeff_e=(a, b)) for (a, b) in sols]
            else:
                if d < 0:
                    sols = diop_DN_negative_d(d, 4 * n)
                else:
                    sols = diop_DN(d, 4 * n)
                elements_norm_n = []
                for (x, y) in sols:
                    if (x - y) % 2 == 0:
                        element = QuadraticInteger(coeff_e=((x - y) / 2, y))
                        elements_norm_n.append(element)
                return elements_norm_n

        # ----- algebraic method -----

        def euclidean_function(self):
            """The value of the `Euclidean function`_ (the absolute
            value of the norm).

                >>> Zi = QuadraticIntegerRing(-1)
                >>> alpha = Zi(1 + sqrt(-1))
                >>> alpha
                1 + I
                >>> alpha.norm
                2
                >>> alpha.euclidean_function()
                2
                >>> O = QuadraticIntegerRing(5)
                >>> beta = O(1 + sqrt(5))
                >>> beta
                1 + sqrt(5)
                >>> beta.norm
                -4
                >>> beta.euclidean_function()
                4

            Warning:
                The ring of quadratic integers must be an Euclidean domain
                with the norm as a Euclidean function, that is,  :math:`d`
                (seq. `A048981`_) must be one of the following:

                    -11, -7, -3, -2, -1, 2, 3, 5, 6, 7, 11, 13, 17, 19, 21, 29, 33, 37, 41,
                    57, 73

            .. _Euclidean function: https://en.wikipedia.org/wiki/Quadratic_integer#Euclidean_rings_of_quadratic_integers
            .. _A048981: https://oeis.org/A048981
            """
            if QuadraticInteger.d < 0:
                return self.norm
            else:
                return Abs(self.norm)

        def is_unit(self):
            """Test whether the quadratic integer is a unit.

                >>> Zi = QuadraticIntegerRing(-1)
                >>> alpha = Zi(sqrt(-1))
                >>> alpha
                I
                >>> alpha.is_unit()
                True
                >>> beta = Zi(1 + sqrt(-1))
                >>> beta
                1 + I
                >>> beta.is_unit()
                False

            """
            return simplify(self.euclidean_function() - 1) == 0

        def is_irreducible(self):
            """Test whether the quadratic integer is irreducible.

                >>> O = QuadraticIntegerRing(-5)
                >>> alpha = O(3)
                >>> alpha.is_irreducible()
                True

            """
            if isprime(self.euclidean_function()):
                return True
            else:
                factors = factorint(self.euclidean_function())
                if len(factors) == 1 and list(factors.values())[0] == 2:
                    l1 = QuadraticInteger.elements_with_norm(list(factors.keys())[0])
                    l2 = QuadraticInteger.elements_with_norm(-(list(factors.keys())[0]))
                    return l1 + l2 == []
                else:
                    return False

        def factor(self):
            """Factor the quadratic integer as a product of irreducible ones.

                >>> Zi = QuadraticIntegerRing(-1)
                >>> alpha = Zi(2)
                >>> alpha
                2
                >>> alpha.factor()
                [1 + I, 1 - I]

            Returns:
                List[QuadraticInteger]: the irreducible factors.

            Warning:
                The ring of quadratic integers must be a unique factorization
                domain. Examples of :math:`d`  (seq. `A048981`_) that makes the
                ring of quadratic integers an UFD are the following:

                    -163, -67, -43, -19, -11, -7, -3, -2, -1, 1, 2, 3, 5, 6, 7, 11, 13, 14,
                    17, 19, 21, 22, 23, 29, 31, 33, 37, 38, 41, 43, 46, 47, 53, 57, 59, 61,
                    62, 67, 69, 71, 73, 77, 83, 86, 89, 93, 94, 97, 101, 103, 107, 109,
                    113, 118, 127, 129, 131, 133, 134, 137, 139, 141, 149, ...

                Otherwise, this method may not work.

            .. _A061574: https://oeis.org/A061574
            """
            if self.is_unit():
                raise ValueError("{0} is an unit".format(self))

            alpha = QuadraticInteger(self.value)  # a copy
            irreducible_elements = []

            while not alpha.is_unit() and not alpha.is_irreducible():
                norm_factors = factorint(alpha.norm)
                if -1 in norm_factors:
                    del norm_factors[-1]
                p1 = sorted(norm_factors.keys())[0]
                L = QuadraticInteger.elements_with_norm(p1)
                L += QuadraticInteger.elements_with_norm(-p1)

                if L == []:
                    alpha1 = QuadraticInteger(p1)
                    try:
                        alpha / alpha1
                    except ValueError:
                        raise ValueError("{} could not be factored.".format(self))
                    else:
                        irreducible_elements.append(alpha1)
                        alpha = alpha / alpha1
                else:
                    for alpha1 in L:
                        try:
                            alpha / alpha1
                        except ValueError:
                            continue
                        else:
                            irreducible_elements.append(alpha1)
                            alpha = alpha / alpha1
                            break
                    else:
                        raise ValueError("{} could not be factored.".format(self))

            if alpha.is_irreducible():
                irreducible_elements.append(alpha)
            else:
                print("{} is the product of the obtained factors".format(self),
                      "and the unit: {}".format(alpha))

            return irreducible_elements

        # ----- special methods -----

        def __eq__(self, other):
            return simplify(self - other) == 0

        def __ne__(self, other):
            return not self.__eq__(other)

        def __add__(self, other):
            try:
                return QuadraticInteger(self.value + other.value)
            except:
                return QuadraticInteger(self.value + other)

        __radd__ = __add__

        def __neg__(self):
            return QuadraticInteger(-self.value)

        def __sub__(self, other):
            try:
                return QuadraticInteger(self.value - other.value)
            except:
                return QuadraticInteger(self.value - other)

        def __rsub__(self, other):
            return -self.__sub__(other)

        def __mul__(self, other):
            try:
                return QuadraticInteger(self.value * other.value)
            except:
                return QuadraticInteger(self.value * other)

        __rmul__ = __mul__

        def __truediv__(self, other):
            try:
                quotient = other.value
            except:
                quotient = other

            return QuadraticInteger(self.value / quotient)

        def __rtruediv__(self, other):
            return QuadraticInteger(other).__truediv__(self)

        def __pow__(self, other):
            return QuadraticInteger(self.value ** other)

        def __str__(self):
            return self.value.__str__()

        __repr__ = __str__

    if d == 0 or d == 1 or not is_square_free(d):
        raise ValueError("d can't be 0, 1 or non square-free.")
    QuadraticInteger.d = d
    QuadraticInteger.__name__ = "O(sqrt({}))".format(QuadraticInteger.d)

    if d % 4 != 1:
        QuadraticInteger.e = sqrt(d)
    else:
        QuadraticInteger.e = Rational(1, 2) + Rational(1, 2) * sqrt(d)

    return QuadraticInteger
