# Christopher Lama
# Special thanks to Dr. Charles Thangaraj

import LU
import MatrixGeneration as m
import matplotlib.pyplot as plt

# console executable
def engine(n:int,h:int,l:int,nu:int,deltaP:int):

    """
   Description:
       Function call used to

   Parameters:
       none

   Returns:
       none
   """

    user_input = 'y'  # meant to be used as a switch that will end the program if the user does not want to make more plots


    A = m.GENERATE(n)  # Generates the A matrix

    x, y_points, avg_absolute, avg_relative = LU.SOLVE(A, n,False, l, deltaP, h, nu, False, SHOW_table=True) # generates x,y points and the table


    return avg_absolute,avg_relative

# engine()