# Christopher Lama
# Special thanks to Dr. Charles Thangaraj
# Code borrowed and adapted from https://github.com/ChristopherL891123/Research--Code-Development-for-Flow-Simulation
# Code borrowed uses Length of artery in meters and Radius of artery is in centimeters

import LU #borrowed
import statistics as stats
import MatrixGeneration #borrowed
import random as r
import threading


n = 1000
A = MatrixGeneration.GENERATE(1000)
U, L = LU.DECOMP(A, n, False)
pressure_list = [] # keep track of normally distributed random values for pressure
viscosity_list = [] # keep track of normally distributed random values for viscosity

output = ["table1", "table2", "table3", "table4"]

def GENERATE(n): # borowed
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
        if i == n - 1:
            break
        i = j + 1
        Matrix[j][i] = -1
        Matrix[i][j] = -1

    return Matrix


def engine(l:float, r:float, Delta_P:float, Nu:float):

#  make UNIVERSAL VAR of MATRIX AND LU FACTORIZATION

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

    return stats.mean(Absolute_error),stats.mean(Relative_error)


# Part One: Generation of values for Radius = 10 centimeters, Length = 0.1 meters

def calc1(outputs):

    count = 0 # to keep track of iteration for table
    # set up table
    header = "|{: ^8} | {: ^30} | {: ^30} | {: ^30} | {: ^30}|".format('Iteration', 'Nu', 'P', 'Avg Absolute Error', 'Avg Relative Error')
    table = ""
    table += header + '\n'  # put newline
    table += ((len(header) + 2) * '-') + '\n'  # make division between header and body of table

    for i in range(10000):
        count += 1
        Delta_P = r.randrange(1,1000) # unit: Pascals
        pressure_list.append(Delta_P)
        Nu = r.randrange(1,1000) # unit: Pa*s
        viscosity_list.append(Nu)

        avg_abs, avg_rel = engine(l=0.1, r=10, Delta_P=Delta_P, Nu=Nu) # unit of l: meters, unit of r: centimeters

        # add row of values to table
        table += "|{: ^8} | {: ^30} | {: ^30} | {: ^30} | {: ^30}|".format(count, Nu, Delta_P, avg_abs, avg_rel) + '\n'

    outputs[0] = table

# Part Two: Generation of values for Radius = 10 centimeters, Length = 0.1 meters
def calc2(outputs):

    count = 0 # to keep track of iteration for table
    # set up table
    header = "|{: ^8} | {: ^30} | {: ^30} | {: ^30} | {: ^30}|".format('Iteration', 'Nu', 'P', 'Avg Absolute Error', 'Avg Relative Error')
    table = ""
    table += header + '\n'  # put newline
    table += ((len(header) + 2) * '-') + '\n'  # make division between header and body of table

    for i in range(10000):
        count += 1
        Delta_P = r.randrange(1,1000) # unit: Pascals
        pressure_list.append(Delta_P)
        Nu = r.randrange(1,1000) # unit: Pa*s
        viscosity_list.append(Nu)

        avg_abs, avg_rel = engine(l=0.1, r=100, Delta_P=Delta_P, Nu=Nu) # unit of l: meters, unit of r: centimeters

        # add row of values to table
        table += "|{: ^8} | {: ^30} | {: ^30} | {: ^30} | {: ^30}|".format(count, Nu, Delta_P, avg_abs, avg_rel) + '\n'

    outputs[1] = table

def calc3(outputs):

    # to keep track of iteration for table
    count = 0
    # set up table
    header = "|{: ^8} | {: ^30} | {: ^30} | {: ^30} | {: ^30}|".format('Iteration', 'Nu', 'P', 'Avg Absolute Error', 'Avg Relative Error')
    table = ""
    table += header + '\n'  # put newline
    table += ((len(header) + 2) * '-') + '\n'  # make division between header and body of table

    for i in range(10000):
        count += 1
        Delta_P = r.randrange(1,1000) # unit: Pascals
        pressure_list.append(Delta_P)
        Nu = r.randrange(1,1000) # unit: Pa*s
        viscosity_list.append(Nu)

        avg_abs, avg_rel = engine(l=0.3, r=10, Delta_P=Delta_P, Nu=Nu) # unit of l: meters, unit of r: centimeters

        # add row of values to table
        table += "|{: ^8} | {: ^30} | {: ^30} | {: ^30} | {: ^30}|".format(count, Nu, Delta_P, avg_abs, avg_rel) + '\n'

    outputs[2] = table

def calc4(outputs):

    # to keep track of iteration for table
    count = 0
    # set up table
    header = "|{: ^8} | {: ^30} | {: ^30} | {: ^30} | {: ^30}|".format('Iteration', 'Nu', 'P', 'Avg Absolute Error', 'Avg Relative Error')
    table = ""
    table += header + '\n'  # put newline
    table += ((len(header) + 2) * '-') + '\n'  # make division between header and body of table

    for i in range(10000):
        count += 1
        Delta_P = r.randrange(1,1000) # unit: Pascals
        pressure_list.append(Delta_P)
        Nu = r.randrange(1,1000) # unit: Pa*s
        viscosity_list.append(Nu)

        avg_abs, avg_rel = engine(l=0.1, r=100, Delta_P=Delta_P, Nu=Nu) # unit of l: meters, unit of r: centimeters

        # add row of values to table
        table += "|{: ^8} | {: ^30} | {: ^30} | {: ^30} | {: ^30}|".format(count, Nu, Delta_P, avg_abs, avg_rel) + '\n'

    outputs[3] = table

# --> Run code <--

# code is multithreaded, all results are written to a list. After all threads have written resulting table of values to the list,
# each item of the list is appended to a file, item by item with individual markers and information about readings for pressure
# and viscosity and their corresponding calculated errors.

task1 = threading.Thread(target=calc1,  args=[output], name='calc1').start()
task2 = threading.Thread(target=calc2, args=[output], name='calc1').start() # argument for function must be inside list
task3 = threading.Thread(target=calc3, args=[output], name='calc1').start()
task4 = threading.Thread(target=calc4, args=[output] , name='calc1').start()

# write outputs to file
while True:

    if output[0] != "table1" and output[1] != "table2" and output[2] != "table3" and output[3] != "table4":

        f = open("outputs.txt", "a")
        f.write("*****TABLE 1*****\n\n LENGTH=10 CENTIMETERS AND RADIUS=0.1 METERS")
        f.write(output[0])
        f.write("\n\n")
        f.close()


        f = open("outputs.txt", "a")
        f.write("*****TABLE 2*****\n\n LENGTH=10 CENTIMETERS AND RADIUS=1 METER")
        f.write(output[1])
        f.write("\n\n")
        f.close()


        f = open("outputs.txt", "a")
        f.write("*****TABLE 3*****\n\n LENGTH=30 CENTIMETERS AND RADIUS=0.1 METERS")
        f.write(output[3])
        f.write("\n\n")
        f.close()


        f = open("outputs.txt", "a")
        f.write("*****TABLE 4*****\n\n LENGTH=30 CENTIMETERS AND RADIUS=1 METER")
        f.write(output[4])
        f.write("\n\n")
        f.close()

        break

    else:
        continue



