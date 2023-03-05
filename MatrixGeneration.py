# Christopher Lama
# Special thanks to Dr. Charles Thangaraj

# Matrix generation algorithm
def GENERATE(n):
    """
    Description:
            Generates a NxN matrix composed of 2s in the diagonal line and -1s in the lower and upper diagonal

    Parameters:

        n: int
            Size of matrix; user-given in main() function as well as in the Graphical User Interface.

    Returns:
        Matrix: list
            The generated A matrix
            """

    Matrix = []
    # set up the matrix
    for a in range(n):
        Matrix.append([])
        for b in range(n):
            Matrix[a].append(0)

    # insert values
    i = 0
    for j in range(0, n):
        Matrix[j][j] = 2
        if i == n-1:
            break
        i = j + 1
        Matrix[j][i] = -1
        Matrix[i][j] = -1

    return Matrix

