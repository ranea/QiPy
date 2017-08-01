"""This module allows to work with ideals of quadratic integer rings.

To compute operations with ideal within the same quadratic ring,
all ideals must be created as :class:`Ideal` objects. To do so:

 1. Create the quadratic integer ring :math:`\mathcal{O}_{\mathbb{Q}[\sqrt{d}]}`
    with the function :any:`QuadraticIntegerRing`.

    >>> O = QuadraticIntegerRing(-5)

 2. Use the returned factory to create the generators of the ideal.

    >>> generator1 = O(3)
    >>> generator2 = O("1 + sqrt(-5)")

 3. Create the ideal object with the generators as arguments and use the
    available operators and methods.

    >>> I = Ideal(generator1, generator2)
    >>> I.factor()
    [<1 + sqrt(5)*I,3*sqrt(5)*I>]

Note that this module, ``ideal``, need to be imported to use
its classes and functions. There are several ways to import it:

 1. Import all functions and classes of QiPy: ::

    >>> from qipy import *
    >>> Zi = QuadraticIntegerRing(-1)
    >>> I = Ideal(Zi(3))

 2. Import only the package, and use the package's identifier to
    access the classes and functions: ::

    >>> import qipy
    >>> Zi = qipy.QuadraticIntegerRing(-1)
    >>> I = qipy.Ideal(Zi(3))

"""
from itertools import product

from sympy import sqrt, simplify, Rational, Abs, isprime, factorint
from sympy import Matrix, symbols, solve, poly
from sympy.polys.numberfields import minimal_polynomial

from qipy.quadratic_integer import QuadraticIntegerRing
from qipy.utilities import lattice_reduce


class Ideal(object):
    """Represent an ideal of a quadratic integer ring.

        >>> Zi = QuadraticIntegerRing(-1)
        >>> I = Ideal(Zi(3))
        >>> I
        <3>
        >>> J = Ideal(Zi("1 + I"))
        >>> J
        <1 + I>
        >>> I * J
        <3 + 3*I,6*I>

    This class supports the operators ``*`` and ``/`` (division only by prime
    ideals) with their natural meaning.

    Note:
        An ideal is represented by a reduced basis of two elements (as
        abelian group, not as ideal).

        However, if at some moment it is known that the ideal is principal,
        its representation is changed with a generator (as ideal).

            >>> Zi = QuadraticIntegerRing(-1)
            >>> I = Ideal(Zi("1 + I"), Zi("2*I"))
            >>> I
            <1 + I,2*I>
            >>> I.is_principal()
            True
            >>> I
            <1 + I>

    Args:
        generators: a sequence of quadratic integers (of the same quadratic ring)
                    that span the ideal.

    Attributes:
        O: the quadratic integer ring related to the ideal.
        norm: the norm of the ideal.
        basis: a list of two quadratic integers that span the ideal.
    """
    def __init__(self, *generators):
        d = generators[0].d
        for g in generators[1:]:
            if g.d != d:
                raise ValueError("Generators must belong to the same quadratic ring")

        self.O = QuadraticIntegerRing(d)

        if len(generators) == 1:
            self._generator = generators[0]
            self.norm = Abs(self._generator.norm)
            self.basis = [self._generator, self._generator * self.O.e]
        else:
            abelian_group_gen = []
            for g in generators:
                if g == 0:
                    raise ValueError("A generator can't be zero.")
                abelian_group_gen.append(g)
                abelian_group_gen.append(g * self.O.e)

            coeff_e = [g.coeff_e for g in abelian_group_gen]
            matrix = Matrix(coeff_e).T
            lattice_reduce(matrix)

            a, b, c = matrix[0, 0], matrix[1, 0], matrix[1, 1]
            self.norm = Abs(a * c)
            self.basis = [self.O(coeff_e=(a, b)), self.O(coeff_e=(0, c))]

    # ----- properties -----

    @property
    def generator(self):
        """A quadratic integer that spans the ideal if it
        is principal. If it isn't, it returns ``None``.

            >>> Zi = QuadraticIntegerRing(-1)
            >>> Ideal(Zi("1 + I"), Zi("2*I")).generator
            1 + I
            >>> O = QuadraticIntegerRing(-5)
            >>> Ideal(O(3), O("1 + sqrt(-5)")).generator is None
            True

        """
        try:
            return self._generator
        except AttributeError:
            O = self.O
            n = self.norm

            if simplify(Abs(n) - 1) == 0:
                self._generator = O(1)
                return self._generator

            L = O.elements_with_norm(n)
            for element in L:
                if self.contain(element):
                    self._generator = element
                    return self._generator

            L = O.elements_with_norm(-n)
            for element in L:
                if self.contain(element):
                    self._generator = element
                    return self._generator

            self._generator = None
            return self._generator

    # ----- class methods ------

    @classmethod
    def prime_divisors(cls, p, d):
        """Return the prime ideal(s) that divides :math:`\langle p \\rangle`
        (the ideal generated by :math:`p`) in the quadratic integer ring
        defined by :math:`d`.

            >>> divisors = Ideal.prime_divisors(3, -1)
            >>> divisors
            [<3>]
            >>> divisors = Ideal.prime_divisors(5, -1)
            >>> divisors
            [<-1 + 2*I,5*I>, <-1 + 3*I,5*I>]

        Returns:
            List[Ideal]: the prime divisors.
        """
        O = QuadraticIntegerRing(d)
        f = minimal_polynomial(O.e, "x")
        fp = poly(f, modulus=p)
        roots = list(fp.ground_roots().keys())
        if roots:
            if len(roots) == 1:
                I = Ideal(O(p), O(O.e - roots[0]))
                return [I]
            else:
                I1 = Ideal(O(p), O(O.e - roots[0]))
                I2 = Ideal(O(p), O(O.e - roots[1]))
                return [I1, I2]
        else:
            return [Ideal(O(p))]

    @classmethod
    def unit_ideal(cls, d):
        """Return the unit ideal of the quadratic integer ring
        defined by :math:`d`.

        Returns:
            Ideal: the unit ideal.
        """
        return Ideal(QuadraticIntegerRing(d)(1))

    # ----- algebraic method -----

    def is_proper(self):
        """Test whether the ideal is not the total.

            >>> Zi = QuadraticIntegerRing(-1)
            >>> Ideal(Zi(1)).is_proper()
            False
            >>> Ideal(Zi(3)).is_proper()
            True
        """
        result = simplify(self.norm - 1) != 0
        if not result:
            self._generator = self.O(1)
        return result

    def is_prime(self):
        """Test whether the ideal is prime.

            >>> Zi = QuadraticIntegerRing(-1)
            >>> Ideal(Zi(2)).is_prime()
            False
            >>> Ideal(Zi(3)).is_prime()
            True
        """
        if not self.is_proper():
            return False
        elif isprime(self.norm):
            return True
        else:
            factors = factorint(self.norm)
            if len(factors) == 1 and list(factors.values())[0] == 2:
                p = list(factors.keys())[0]
                fp = poly(minimal_polynomial(self.O.e, "x"), modulus=p)
                roots = fp.ground_roots().keys()
                return len(roots) == 0
            else:
                return False

    def is_principal(self):
        """Test whether the ideal is principal.

            >>> Zi = QuadraticIntegerRing(-1)
            >>> Ideal(Zi("1 + I"), Zi("2*I")).is_principal()
            True
            >>> O = QuadraticIntegerRing(-5)
            >>> Ideal(O(3), O("1 + sqrt(-5)")).is_principal()
            False
        """
        return self.generator is not None

    def contain(self, element):
        """Test whether the ideal contains an element.

            >>> Zi = QuadraticIntegerRing(-1)
            >>> I = Ideal(Zi(3))
            >>> I.contain(Zi(3))
            True
            >>> I.contain(Zi(1))
            False

        Note that the element must be a quadratic integer of the same ring
        as the ideal.
        """
        x, y = symbols("x, y")
        a, b = self.basis[0].coeff_e
        c = self.basis[1].coeff_e[1]
        a_, b_ = element.coeff_e
        sol = solve((a * x - a_, b * x + c * y - b_), x, y)
        a__, b__ = sol[x], sol[y]

        return a__.is_integer and b__.is_integer

    def divide(self, other):
        """Test whether the ideal divides other ideal.

            >>> Zi = QuadraticIntegerRing(-1)
            >>> I = Ideal(Zi("1 + I"))
            >>> I.divide(Ideal(Zi(2)))
            True
        """
        for g in other.basis:
            if not self.contain(g):
                return False
        return True

    def factor(self):
        """Factor the ideal as a product of prime ideals.

            >>> Zi = QuadraticIntegerRing(-1)
            >>> I = Ideal(Zi(2))
            >>> I
            <2>
            >>> factors = I.factor()
            >>> factors
            [<1 + I,2*I>, <1 + I,2*I>]
            >>> factors[0] * factors[1]
            <-2,2*I>

        Returns:
            List[Ideal]: the prime factors.
        """
        if not self.is_proper():
            raise ValueError("{0} is not a proper ideal".format(self))
        elif self.is_prime():
            return [self]

        I = Ideal(*self.basis)  # a copy
        prime_ideals = []

        while I.is_proper() and not I.is_prime():
            norm_factors = factorint(I.norm)
            p1 = sorted(norm_factors.keys())[0]
            L = Ideal.prime_divisors(p1, self.O.d)

            for P in L:
                if P.divide(I):
                    prime_ideals.append(P)
                    I = I / P
                    break

        if I.is_prime():
            prime_ideals.append(I)

        return prime_ideals

    # ----- special methods -----

    def __eq__(self, other):
        if simplify(self.norm - other.norm) != 0:
            return False
        else:
            return self.divide(other) and other.divide(self)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __mul__(self, other):
        generators = []
        for alpha, beta in product(self.basis, other.basis):
            generators.append(alpha * beta)
        return Ideal(*generators)

    def __truediv__(self, other):
        if not other.is_prime():
            raise NotImplementedError("Only supported for prime ideal denominators.")

        if sqrt(other.norm).is_integer:
            other_inverse_gen = [Rational(1, sqrt(other.norm))]
        else:
            p = other.norm
            fp = poly(minimal_polynomial(self.O.e, "x"), modulus=p)
            roots = list(fp.ground_roots().keys())

            if other.contain(self.O(self.O.e - roots[0])):
                g = self.O.e - roots[-1]
            else:
                g = self.O.e - roots[0]

            other_inverse_gen = [1, Rational(1, p) * g]

        generators = []
        for alpha, beta in product(self.basis, other_inverse_gen):
            generators.append(alpha * beta)
        return Ideal(*generators)

    def __str__(self):
        try:
            principal = self._generator
        except AttributeError:
            principal = None

        if principal is not None:
            return "<{}>".format(self.generator)
        else:
            g1, g2 = self.basis
            return "<{},{}>".format(g1, g2)

    __repr__ = __str__

    def __hash__(self):
        alpha, beta = self.basis
        return hash((alpha.coeff_e, beta.coeff_e))
