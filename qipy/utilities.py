from sympy import sqrt, simplify, Abs, pi, factorial, Matrix
from sympy.functions import sign
from sympy.ntheory.factor_ import core


def is_square_free(integer):
    """Test whether the input is a square-free integer.

        >>> is_square_free(2 * 3)
        True
        >>> is_square_free(2 * 2)
        False
    """
    square_free_part = core(Abs(integer), 2) * sign(integer)
    return simplify(square_free_part - integer) == 0


def minkowski_bound(d):
    """Compute the `Minkowski bound`_ for a given :math:`d`.

    .. _Minkowski bound: https://en.wikipedia.org/wiki/Minkowski%27s_bound
    """
    if d == 0:
        raise ValueError("d can't be zero.")
    elif d > 0:
        s = 2
        t = 0
    else:
        s = 0
        t = 1

    n = s + 2 * t
    M = ((4 / pi)**t) * (factorial(n) / (n**n))

    if d % 4 == 1:
        discriminant = d
    else:
        discriminant = 4 * d

    return M * sqrt(Abs(discriminant))


def lattice_reduce(matrix):
    """Modify the input matrix to obtain a reduced basis of cardinal 2.

    Given a matrix of the form (where (ai,bi) is a vector))
        Matrix([
            [a1, a2, ..., an],
            [b1, b2, ..., bn]
        ])
    , it returns a matrix with the form
        Matrix([
            [a, 0, ..., 0],
            [b, c, ..., 0]
        ])
    where the reduced basis is {(a,b), (0,c)}.
    """
    def basic_reduction(matrix):
        """Basic reduction step.

        Given a matrix of the form (where (ai,bi) is a quadratic integer)):
            Matrix([
                [a1, a2, ..., an],
                [b1, b2, ..., bn]
            ])
        it returns a matrix with the form
            Matrix([
                [a1', a2', ..., bn'],
                [b1', b2', ..., bn']
            ])
        where a1'=min(a1, ..., an) and ai' = ai % a1'.
        """
        # search the column with the minimun element (w.r.t abs)
        minimun = None
        index_min = 0
        for index, element in enumerate(matrix.row(0)):
            if element == 0:
                continue
            if minimun is None:
                minimun = element
                index_min = index
            elif Abs(element) < Abs(minimun):
                minimun = element
                index_min = index

        # swap this column with the first column
        matrix.col_swap(0, index_min)

        # apply elemental transformations
        a1 = matrix[0, 0]
        b1 = matrix[1, 0]
        for i in range(1, matrix.shape[1]):
            ai = matrix[0, i]
            q, r = divmod(ai, a1)
            matrix[0, i] = r
            matrix[1, i] -= q * b1

    # get [[a, 0,..., 0],[b, b2',..., bn']]
    while any(matrix[0, 1:]):
        basic_reduction(matrix)

    # submatrix = [[b2',...,bn'],[0,...0]]
    submatrix = matrix[:, 1:]  # withouth a1, b1
    submatrix.row_swap(0, 1)

    # get [[c,...,0],[0,...0]]
    while any(submatrix[0, 1:]):
        basic_reduction(submatrix)

    # insert submatrix to matrix
    submatrix.row_swap(0, 1)
    matrix[:, 1:] = submatrix

    matrix[1, 0] = matrix[1, 0] % matrix[1, 1]  # |c| > |b|


if __name__ == '__main__':
    M = Matrix([[-3, 7, 6, 5], [3, 4, 1, 6]])
    lattice_reduce(M)
    print(M)  # [[1, 0, 0, 0], [0, 1, 0, 0]]

    M = Matrix([[0, 4, 0, 2], [2, 0, 1, 0]])
    lattice_reduce(M)
    print(M)  # [[0, 1, 0, 0], [2, 0, 0, 0]]
