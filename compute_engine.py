# Christopher Lama
# Special thanks to Dr. Charles Thangaraj
# Code borrowed and adapted from https://github.com/ChristopherL891123/Research--Code-Development-for-Flow-Simulation
# Code borrowed uses Length of artery in meters and Radius of artery is in centimeters
import matplotlib.pyplot as plt

import LU  # borrowed
import statistics as stats
import MatrixGeneration  # borrowed
import matplotlib as m
import random as r
import threading
import sys
import os

print("*=*=*=*=* INITIALIZING VALUES *=*=*=*=*")
n = 10
A = MatrixGeneration.GENERATE(n)
U, L = LU.DECOMP(A, n, False)
pressure_list = []  # keep track of normally distributed random values for pressure
viscosity_list = []  # keep track of normally distributed random values for viscosity
rel_error_list = [0,0,0,0]
abs_error_list = [0,0,0,0]
output = ["table", "table", "table", "table"]
progress_bar = ""  # append all progress of threads and print this repeatedly to the screen
print("*=*=*=*=* DONE INITIALIZING VALUES *=*=*=*=*")


def engine(l: float, r: float, Delta_P: float, Nu: float):
    """

    Args:
        l: length of artery
        r: radius of artery
        Delta_P: pressure
        Nu: viscosity

    Returns:
        Absolute_error: list with calculated absolute errors for the table
    """

    # make UNIVERSAL VAR of MATRIX AND LU

    global n
    global A
    global U
    global L

    # calculate values using formulas
    Delta_Y = (2 * r) / (n + 1)
    V_max = (Delta_P * r ** 2) / (2 * Nu * l)
    B = [(Delta_Y ** 2 * (2 * V_max)) / (r ** 2) for i in range(n)]
    EV = []
    Y_j_LIST = []

    for i in range(0, n + 2):
        Y_j = -r + i * Delta_Y
        Y_j_LIST.append(Y_j)
        EV.append(V_max * (1 - (Y_j / r) ** 2))

    # no need to generate A and calculate U,L as it will always be re-calculated every iteration of the program.
    # For size 1000, these matrices will always have the same values.

    y = LU.FORWARD_SUB(n, L, B, False)
    x = LU.BACKWARD_SUB(y, n, U, False)

    Absolute_error = [0.0]
    Relative_error = [0.0]
    for i in range(1, n + 1):
        Absolute_error.append(abs(x[i] - EV[i]))
        Relative_error.append(abs(Absolute_error[i] / EV[i]))
    Absolute_error.append(0.0)
    Relative_error.append(0.0)

    return stats.mean(Absolute_error), stats.mean(Relative_error)


# Part One: Generation of values 

def calc(outputs: list, index_writes: int, rel_write_index: int, abs_write_index: int, length: int, radius: int):
    global progress_bar
    global rel_error_list
    global abs_error_list

    count = 0  # to keep track of iteration for table
    # set up table
    header = "|{: ^8} | {: ^30} | {: ^30} | {: ^30} | {: ^30}|".format('Iteration', 'Nu', 'P', 'Avg Absolute Error',
                                                                       'Avg Relative Error')
    table = ""
    table += header + '\n'  # put newline
    table += ((len(header) + 2) * '-') + '\n'  # make division between header and body of table
    prg_count = 0  # for the progress bar; reset and prints "|" char every 1000 iterations

    temp_rel_err = []
    temp_abs_err = []

    for i in range(10000):

        if prg_count == 1000:
            progress_bar += "|"
            prg_count = 0
        prg_count += 1

        count += 1
        Delta_P = r.randrange(1, 1000)  # unit: Pascals
        pressure_list.append(Delta_P)
        Nu = r.randrange(1, 1000)  # unit: Pa*s
        viscosity_list.append(Nu)

        avg_abs, avg_rel = engine(l=length, r=radius, Delta_P=Delta_P, Nu=Nu)  # unit of l: meters, unit of r: centimeters

        temp_abs_err.append(avg_abs)
        temp_rel_err.append(avg_rel)

        # add row of values to table
        table += "|{: ^8} | {: ^30} | {: ^30} | {: ^30} | {: ^30}|".format(count, Nu, Delta_P, avg_abs, avg_rel) + '\n'

    rel_error_list[rel_write_index] = temp_rel_err
    abs_error_list[abs_write_index] = temp_abs_err

    outputs[index_writes] = table



# --> Run code <--

# code is multithreaded, all results are written to a list. After all threads have written resulting table of values to the list,
# each item of the list is appended to a file, item by item with individual markers and information about readings for pressure
# and viscosity and their corresponding calculated errors.

task1 = threading.Thread(target=calc, args=[output, 0, 0,0, 0.1, 10], name='calc1').start() #for Radius = 10 centimeters, Length = 0.1 meters
task2 = threading.Thread(target=calc, args=[output, 1, 1,1  ,0.1, 100], name='calc2').start() # for Radius = 10 centimeters, Length = 1 meter
task3 = threading.Thread(target=calc, args=[output, 2, 2,2 ,0.3, 10], name='calc3').start() # for Radius = 30 centimeters, Length = 0.1 meters
task4 = threading.Thread(target=calc, args=[output, 3, 3,3 , 0.3, 100], name='calc4').start() #for Radius = 30 centimeters, Length = 1 meter

 # argument for function must be inside list

# write outputs to file
while True:

    os.system("cls")

    print("\033[92m {}\033[00m".format("*=*=*=*=* PROGRESS: " + progress_bar + " *=*=*=*=*"))
# TODO CHECK IF ALL OF THE VALUES IN TEH LIST HAVE BEEN MODIFIED
    if output[0] != "table1" and output[1] != "table2" and output[2] != "table3" and output[3] != "table4":

        file = open("outputs.txt", "a")
        file.write("*****TABLE 1*****\n\n LENGTH=10 CENTIMETERS AND RADIUS=0.1 METERS\n")
        file.write(output[0])
        file.write("\n\n")
        file.close()

        file = open("outputs.txt", "a")
        file.write("*****TABLE 2*****\n\n LENGTH=10 CENTIMETERS AND RADIUS=1 METER\n")
        file.write(output[1])
        file.write("\n\n")
        file.close()

        file = open("outputs.txt", "a")
        file.write("*****TABLE 3*****\n\n LENGTH=30 CENTIMETERS AND RADIUS=0.1 METERS\n")
        file.write(output[2])
        file.write("\n\n")
        file.close()



#TODO run 10 diamteres (choose range) and 10 lengths (100 tables) , save the plots

        if input("Exit(y/n): ").lower() == "y":
            sys.exit(0)

    else:
        continue

def graph(r:int, l:int, output_list:list, output_index, table_num:int):
    file = open("outputs.txt", "a")
    file.write("*****{TABLE}*****\n\n LENGTH= {a} CENTIMETERS AND RADIUS={b} METER\n".format(TABLE=str(table_num), a=l, b=r))
    file.write(output_list[output_index])
    file.write("\n\n")
    file.close()


