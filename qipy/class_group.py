"""This module allows to work with `ideal classes`_ and compute the
`class numbers`_ of quadratic integer rings.

To compute operations with ideal classes within the same quadratic ring,
all ideal classes must be created as :any:`IdealClass` objects. To do so:

 1. Create the quadratic integer ring :math:`\mathcal{O}_{\mathbb{Q}[\sqrt{d}]}`
    with the function :any:`QuadraticIntegerRing`.

    >>> O = QuadraticIntegerRing(-5)

 2. Create the representative ideal of the ideal class.

    >>> generator1 = O(3)
    >>> generator2 = O("1 + sqrt(-5)")
    >>> I = Ideal(generator1, generator2)

 3. Create the ideal class object with the ideal as argument and use the
    available operators and methods.

    >>> a = IdealClass(I)
    >>> a.order
    2

To compute the `class group`_ of a quadratic integer ring (i.e. the class number
and the generators of the class group), simply instance the class
:any:`ClassGroup` and access its attributes :any:`class_number` and
:any:`generators`.

Note that this module, ``class_group``, need to be imported to use
its classes and functions. There are several ways to import it:

 1. Import all functions and classes of QiPy: ::

    >>> from qipy import *
    >>> O = QuadraticIntegerRing(-5)
    >>> a = IdealClass( Ideal(O(2), O("1 + sqrt(-5)")) )

 2. Import only the package, and use the package's identifier to
    access the classes and functions: ::

    >>> import qipy
    >>> O = qipy.QuadraticIntegerRing(-5)
    >>> I = qipy.Ideal(O(2), O("1 + sqrt(-5)"))
    >>> a = qipy.IdealClass(I)

.. _class numbers: http://mathworld.wolfram.com/ClassNumber.html
.. _class group: https://en.wikipedia.org/wiki/Ideal_class_group
.. _ideal classes: https://en.wikipedia.org/wiki/Ideal_class_group
"""
from itertools import product
from itertools import combinations

from sympy import floor, primerange, simplify, igcd, factorint

from qipy.ideal import Ideal
from qipy.quadratic_integer import QuadraticIntegerRing
from qipy.utilities import minkowski_bound


class IdealClass(object):
    """Represent an ideal class of a quadratic integer ring.

        >>> Zi = QuadraticIntegerRing(-1)
        >>> I = Ideal(Zi(3))
        >>> a = IdealClass(I)
        >>> a
        [<3>]
        >>> a == a * a
        True
        >>> O = QuadraticIntegerRing(-5)
        >>> J = Ideal(O(2), O("1 + sqrt(-5)"))
        >>> b = IdealClass(J)
        >>> b
        [<1 + sqrt(5)*I,2*sqrt(5)*I>]
        >>> b ** 2
        [<-2,2*sqrt(5)*I>]

    This class supports the operators ``*`` and ``**`` with their natural meaning.

    Args:
        ideal: a representative of the ideal class
        name: an optional name for the representative.

    Attributes:
        O: the quadratic integer ring related to the ideal class.
        representative: a representative of the ideal class.
    """

    def __init__(self, ideal, name=None):
        self.representative = ideal
        self.O = ideal.O

        self._powers = [self.representative]

        if name is not None:
            self.name = name

    # ----- properties -----

    @property
    def order(self):
        """The `order`_ of the ideal class, that is, the smallest positive
        intenger such that ``self ** order`` is the trivial class.

            >>> Zi = QuadraticIntegerRing(-1)
            >>> a = IdealClass(Ideal(Zi(3)))
            >>> a
            [<3>]
            >>> a.order
            1
            >>> O = QuadraticIntegerRing(-5)
            >>> b = IdealClass(Ideal(O(2), O("1 + sqrt(-5)")))
            >>> b
            [<1 + sqrt(5)*I,2*sqrt(5)*I>]
            >>> b.order
            2

        .. _order: https://en.wikipedia.org/wiki/Order_(group_theory)
        """
        try:
            return self._order
        except:
            if self.representative.is_principal():
                self.representative = Ideal.unit_ideal(self.O.d)
                self._order = 1
                return self._order

            I = self.representative
            i = 2  # order candidate
            while True:
                product = self._powers[(i - 1) - 1] * I

                if product.is_principal():
                    self._order = i
                    return self._order

                self._powers.append(product)
                i += 1

    # ----- class methods ------

    @classmethod
    def trivial_class(cls, d):
        """Return the trivial class (with the unit ideal as representative) of
        the quadratic integer ring defined by :math:`d`.

            >>> IdealClass.trivial_class(-1)
            [<1>]

        Returns:
            IdealClass: the trivial ideal class.
        """
        return cls(Ideal.unit_ideal(d))

    # ----- algebraic method -----

    def is_trivial(self):
        """Test whether the ideal class is the trivial class.

            >>> Zi = QuadraticIntegerRing(-1)
            >>> a = IdealClass(Ideal(Zi(3)))
            >>> a
            [<3>]
            >>> a.is_trivial()
            True
            >>> O = QuadraticIntegerRing(-5)
            >>> b = IdealClass(Ideal(O(2), O("1 + sqrt(-5)")))
            >>> b
            [<1 + sqrt(5)*I,2*sqrt(5)*I>]
            >>> b.is_trivial()
            False

        """
        try:
            return self._order == 1
        except AttributeError:
            result = self.representative.is_principal()
            if result:
                self._order = 1
            return result

    def inverse(self):
        """Return the inverse of the ideal class.

            >>> Zi = QuadraticIntegerRing(-1)
            >>> a = IdealClass(Ideal(Zi(3)))
            >>> a
            [<3>]
            >>> a.inverse()
            [<1>]
            >>> O = QuadraticIntegerRing(-5)
            >>> b = IdealClass(Ideal(O(2), O("1 + sqrt(-5)")))
            >>> b
            [<1 + sqrt(5)*I,2*sqrt(5)*I>]
            >>> b.inverse()
            [<1 + sqrt(5)*I,2*sqrt(5)*I>]
        """
        return self ** (self.order - 1)

    # ----- special methods -----

    def __eq__(self, other):
        if self.is_trivial():
            return other.is_trivial()
        elif other.is_trivial():
            return False

        try:
            order = self._order
        except AttributeError:
            order = None

        if order is None:
            return (self * other.inverse()).is_trivial()
        else:
            return (self.inverse() * other).is_trivial()

    def __mul__(self, other):
        I = self.representative
        J = other.representative
        return IdealClass(I * J)

    def __pow__(self, exponent):
        if exponent < 0:
            raise ValueError("Exponent must be non negative.")
        elif exponent == 0:
            return IdealClass.trivial_class(self.O.d)
        elif exponent == 1:
            return self

        try:
            order = self._order
        except AttributeError:
            order = None

        if order is not None:
            if self.is_trivial():
                return self

            exponent = exponent % order
            if exponent == 0:
                return IdealClass.trivial_class(self.O.d)
            if exponent == 1:
                return self
            else:
                power = IdealClass(self._powers[exponent - 1])

                # extract power's order from self
                power_order = 2
                new_exponent = exponent
                while True:
                    new_exponent = (new_exponent + exponent) % order
                    if new_exponent == 0:
                        power._order = power_order
                        break
                    else:
                        power._powers.append(self._powers[new_exponent - 1])
                        power_order += 1
                return power
        else:
            I = self.representative
            product = I * I
            for i in range(2, exponent):
                product *= I
            return IdealClass(product)

    def __str__(self):
        try:
            return "[{0}]".format(self.name)
        except AttributeError:
            return "[{0}]".format(self.representative)

    __repr__ = __str__


class ClassGroup(object):
    """Represent the class group of a quadratic integer ring.

        >>> G = ClassGroup(-1)
        >>> G.class_number
        1
        >>> G.generators
        [[<1>]]
        >>> H = ClassGroup(-14, verbose=True)  # doctest: +SKIP
        List of generators (including dependent): [[p2], [p3]]
        Testing whether [p3] is dependent:
            [p3] is independent
        Testing whether [p2] is dependent:
            Dependency relation found: [<1>] == [p2] * [p3]^2
            [p2] removed from the list of generators
        Class number: 4
        Generators: [[p3]]

    The non-trivial generators are represented as ``[pX]`` where ``pX`` is
    a prime ideal that divides :math:`\langle X \\rangle` (the ideal
    generated by :math:`X`).

    Warning:
        In the of the object, the generators and the class number are computed.
        Therefore, the creation of the object may take some time.

    Args:
        d: a non-square free integer that defines the quadratic integer ring.
        verbose: if ``True``, information about the computation process
            of the class group is printed.

    Attributes:
        O: the quadratic integer ring related to the ideal class.
        generators: a list of ideal classes that span the class group.
        class_number: the order of the class group.
    """
    def __init__(self, d, verbose=False):
        self.O = QuadraticIntegerRing(d)
        self.generators = None
        self.class_number = None
        self.verbose = verbose
        self._compute_class_group()

    def _find_generators(self):
        """Compute a set of generators (probably dependent) of the class group."""
        bound = minkowski_bound(self.O.d)

        if bound < 2:
            self.generators = [IdealClass.trivial_class(self.O.d)]
            return

        stop = floor(bound) + 1

        self.generators = []
        primes = list(primerange(2, stop))
        for p in primes:
            I = Ideal.prime_divisors(p, self.O.d)[0]
            if not I.is_principal():
                g = IdealClass(I, "p{0}".format(p))
                self.generators.append(g)

        if len(self.generators) == 0:
            self.generators.append(IdealClass.trivial_class(self.O.d))

    def _is_generator_dependent(self, index_generator):
        """Test whether self.generator[index_generator] is dependent, that is,
        the group generated by self.generator[index_generator] (or a subgroup)
        can be generated with the rest of generators.
        """
        generator = self.generators[index_generator]
        others = []
        for i in range(len(self.generators)):
            if i != index_generator:
                others.append(self.generators[i])

        order_factors = factorint(generator.order).keys()
        exponent_combinations = product(*[range(0, g.order) for g in others])
        next(exponent_combinations)  # skip (0, ..., 0)

        if self.verbose:
            print("Testing whether {} is dependent:".format(generator))

        for exp_combination in exponent_combinations:
            class_product = IdealClass.trivial_class(self.O.d)
            msg = ""
            for g, e in zip(others, exp_combination):
                class_product *= g ** e
                msg += " * {}^{}".format(g, e)

            if (class_product * generator).is_trivial():
                if self.verbose:
                    msg = "[<1>] == {}{}".format(generator, msg)
                    print("\tDependency relation found: " + msg)
                    print("\t{} removed from the list of generators".format(
                        generator))
                return True, None

            if len(order_factors) > 1:
                for f in order_factors:
                    generator_power = generator ** (generator.order - f)
                    if (class_product * generator_power).is_trivial():
                        new_generators = []
                        for new_f in order_factors:
                            if new_f != f:
                                new_gen = generator ** new_f
                                new_gen.name = "{}^{}".format(generator.name, new_f)
                                new_generators.append(new_gen)
                        if self.verbose:
                            msg = "{}^{} == {}".format(generator, f, msg[3:])
                            print("\tDependency relation found: " + msg)
                            print("\t{} removed from the list of generators".format(
                                generator))
                            for new_gen in new_generators:
                                print("\t{} added to the list of generators".format(
                                    new_gen))
                        return True, new_generators

        if self.verbose:
            print("\t{} is independent".format(generator))
        return False, None

    def _are_generator_independent(self):
        """Do a one-way test to check if the generators are independent,
        that is, if the test successful, the generators are independent.
        Otherwise, they may or may not be dependent."""
        if len(self.generators) == 1:
            return True

        orders = [g.order for g in self.generators]
        for combination in combinations(orders, 2):
            are_independent = simplify(igcd(*combination) - 1) == 0
            if not are_independent:
                return False
        return True

    def _compute_class_group(self):
        """Compute the generators and the class number of the class group.

        First, it obtains a set of generators. Then, it removes those
        which are dependent.
        """
        self._find_generators()

        if self.verbose:
            print("List of generators (including dependent): {}".format(
                self.generators
            ))

        i = 1
        while i <= len(self.generators) and not self._are_generator_independent():
            index = len(self.generators) - i
            is_dependent, new_generators = self._is_generator_dependent(index)
            if is_dependent:
                if new_generators is None:
                    del self.generators[index]
                else:
                    self.generators[index:index + 1] = new_generators
            else:
                i += 1

        self.class_number = 1
        for g in self.generators:
            self.class_number *= g.order

        if self.verbose:
            print("Class number: {}".format(self.class_number))
            print("Generators: {}".format(self.generators))
