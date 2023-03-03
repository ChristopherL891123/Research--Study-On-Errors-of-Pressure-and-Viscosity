# Christopher Lama

import LU
import MatrixGeneration as m
import matplotlib.pyplot as plt

# console executable
def engine(n:int,h:int,l:int,nu:int,deltaP:int):

    """
   Description:
       Console executable that gathers the dimensions of the matrix form the user and solves for the matrix. Then graphs the results using matplotlib.

   Parameters:
       none

   Returns:
       none
   """

    try:

        user_input = 'y'  # meant to be used as a switch that will end the program if the user does not want to make more plots

        while user_input == 'y':
            A = m.GENERATE(n)  # Generates the A matrix

            x, y_points, avg_absolute, avg_relative = LU.SOLVE(A, n,False, l, deltaP, h, nu, False, SHOW_table=True) # generates x,y points and the table

            # set up the graph
            # plt.margins(x=0, y=0)
            # plt.plot(x, y_points)
            # plt.title("Graph for {i} discrete points".format(i=n + 2))
            # plt.xlabel("Velocity", fontsize=12)
            # plt.ylabel("y", rotation="horizontal", fontsize=12)
            # plt.show()

            # user_input = input("Continue plotting? y/n : ")

    except:
        print("ERROR: values provided caused an error")
        import traceback as t
        t.print_exc()

    return avg_absolute,avg_relative

# engine()