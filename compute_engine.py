# Christopher Lama
# Special thanks to Dr. Charles Thangaraj and Dr. Arati Pati
# Code borrowed and adapted from https://github.com/ChristopherL891123/Research--Code-Development-for-Flow-Simulation
# Code borrowed uses Length of artery in meters and Radius of artery is in centimeters


import LU  # borrowed
import statistics as stats
import MatrixGeneration  # borrowed
import matplotlib.pyplot as plt
import random as r
import threading


print("*=*=*=*=* INITIALIZING VALUES *=*=*=*=*")
n = 1000 # size of matrix
A = MatrixGeneration.GENERATE(n) # generate the matrix
U, L = LU.DECOMP(A, n, False) # decompose it
pressure_list = []  # keep track of normally distributed random values for pressure
viscosity_list = []  # keep track of normally distributed random values for viscosity
rel_error_list = [0 for i in range(100)] # list of placeholders for average relative error for each run. Will be made up of nested lists where each nested list is the average for each run
abs_error_list = [0 for i in range(100)] # list of placeholders for average absolute error for each run. Will be made up of nested lists where each nested list is the average for each run
output = [" " for i in range(100)]
repeat = [] # used in cont()
print("*=*=*=*=* DONE INITIALIZING VALUES *=*=*=*=*")

# Part One: Generation of values
def engine(l: float, r: float, Delta_P: float, Nu: float):

    """
    Description: used to compute and return two float values: average relative and absolute error.
    Will be called repeatedly to calculate list of absolute and relative error averages for plotting.
    Uses code borrowed and outfitted from "Code Development for Flow Simulation" by Christopher Lama, Advisor: Dr. Arati Nanda Pati

    Args:
        l: length of artery
        r: radius of artery
        Delta_P: pressure
        Nu: viscosity

    Returns:
        Absolute_error: list with calculated absolute errors for the table
    """

    # use A and U,L as global variables to avoid unnecessary re-computation

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



def calc(outputs: list, index_writes: int, rel_write_index, abs_write_index, length: int, radius: int):
    """
        Description:

        Args:
            outputs: used to store table outputs.
            index_writes: tells the program at what index of outputs var to write the generated table to.
            rel_write_index: tells the program at what index of rel_error_list var to write the generated average relative errors to.
            abs_write_index: tells the program at what index of abs_error_list var to write the generated average absolute errors to.
            length: length of the artery.
            radius: radius of the artery.

        Returns: None

        """

    global progress_bar
    global rel_error_list
    global abs_error_list

    count = 0  # to keep track of iteration for table

    # set up table
    header = "|{: ^8} | {: ^30} | {: ^30} | {: ^30} | {: ^30}|".format('Iteration', 'Nu', 'P', 'Avg Absolute Error',
                                                                       'Avg Relative Error')
    table = "***** {a} *****\n\n LENGTH={b} CENTIMETERS AND RADIUS={c} METERS\n\n\n".format(a="TABLE "+str(index_writes), b=str(length), c=str(radius))
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

# *=*=*=*=* Run code *=*=*=*=*

# code is multithreaded, all results are written to a list. After all threads have written resulting table of values to the list,
# each item of the list is appended to a file, item by item with individual markers and information about readings for pressure
# and viscosity and their corresponding calculated errors. Then Absolute vs Relative errors are plotted and saved to an image file each

def main():
    """
        Description: main executable program. Makes threads to efficiently iterate many times through the calc function
        to generate the tables and data for plots for different values for lengths and radius.
        Writes to disk the tables and saves the plots as images.

        Args:
            None

        Returns:
            None
        """

    write_index = 0
    temp_index_tables = 0
    c2 = 0
    length_range = [i for i in range(10, 110, 10)]
    radius_range = [(i + 1) / 10.00 for i in range(0, 10)]  # avoid the initial 0.0 # can't do floats in range() https://stackoverflow.com/questions/7267226/range-for-floats

    count = 0
    for i in range(10):
        l_ = length_range[i]
        for j in range(10):
            r_ = radius_range[j]
            task1 = threading.Thread(target=calc, args=[output, temp_index_tables, write_index, write_index, l_, r_], name='thread').start()
            count += 1
            write_index += 1
            temp_index_tables += 1

    index = 0
    while True:
        if index >= 100:
            exit(0)
            print([char for char in "FINISHED"])

        if output[index] != " ":

            file = open("outputs.txt", "a")
            file.write(output[index])
            file.close()
            print("wrote to file")
            plt.plot(abs_error_list[index],rel_error_list[index])
            plt.savefig("{a}.png".format(a=str(c2)), bbox_inches="tight")  # https://stackoverflow.com/questions/9622163/save-plot-to-image-file-instead-of-displaying-it-using-matplotlib
            plt.clf() #https://www.tutorialspoint.com/how-do-i-close-all-the-open-pyplot-windows-matplotlib

            c2 += 1
            index += 1

        else:
            continue

main()
