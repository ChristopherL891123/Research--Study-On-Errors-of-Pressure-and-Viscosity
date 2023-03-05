# Christopher Lama
# Special thanks to Dr. John Starner

import MatrixGeneration

# Decomposition algorithm
def DECOMP(A, n, SHOW_LU):
    """
    Description:
        Decomposes the A matrix into L and U matrices.

    Parameters:
        A: list
            The generated matrix (the system of linear equations)

        n: int
            Size of matrix; user-given in main() function as well as in the Graphical User Interface.

        SHOW_LU: bool
            Determines whether the resulting L and U matrices will be printed out.

    Returns:
        U: list
            The U matrix

        L: list
            The L matrix
    """

    # set up U and L
    L = []

    for i in range(n):
        L.append([])
        for j in range(n):
            L[i].append(0)

    for i in range(n):
        L[i][i] = 1

    U = A.copy() # can be made multithread

    for j in range(0, n - 1):  # j is column and represents the diagonal element
        for i in range(j + 1, n):  # i is row
            factor = U[i][j] / U[j][j]
            L[i][j] = factor
            # adjust the row
            for k in range(0, n):
                U[i][k] = U[i][k] - (factor * U[j][k])

    return U, L


# forward substitution algorithm
def FORWARD_SUB(n, L, B, SHOW_y):
    """
    Description:
        Performs forward substitution on the L matrix using the b vector to find the y vector.

    Parameters:
        L: list
            The calculated L matrix.

        n: int
            Size of matrix; user-given in main() function as well as in the Graphical User Interface.

        B: list
            B is the known right-hand side of the original matrix equation.

        SHOW_y: bool
            Detemrines whether the calculated y vector is to be displayed.

    Returns:
        y: list
            The y vector
        """


    # set up y
    y = []
    for i in range(n):
        y.append(0)

    for i in range(0, n):  # i is the row
        sum_row = 0
        for j in range(0, i):  # j is the column; starts from index 0 and ends the index before diagonal element
            sum_row += L[i][j] * y[j]  # get sum of elements from index 0 to element before the diagonal line
        y[i] = (-1 * sum_row) + B[i]  # calcualte the unknown value

    # print y vector
    if SHOW_y:
        print("y = ", y)

    return y


# backward substitution algorithm
def BACKWARD_SUB(y, n, U, SHOW_x):
    """
    Description:
            Performs backward substitution on the U matrix using the y vector to find the x vector

    Parameters:
        U: list
            The calculated U matrix.

        n: int
            Size of matrix; user-given in main() function as well as in the Graphical User Interface.

        y: list
            y vector calculated through FORWARD_SUB()

        SHOW_x: bool
            Detemrines whether the calculated x vector is to be displayed.

    Returns:
        x: list
            The x vector
    """



    x = [0 for i in range(n + 1)]  # n+1 because of boundary conditions: velocity is zero at the walls

    for i in range(-1, -n - 1, -1):  # start from last row
        sum_row = 0
        for j in range(i, 0):  # start from diagonal element and move to last element of the row
            sum_row += U[i][j] * x[j]  # get sum of row
        x[i] = ((-1 * sum_row) + y[i]) / U[i][i]  # calculate the unknown x value for the row
    x.append(0)  # boundary conditions

    return x



