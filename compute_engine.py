# Christopher Lama
# Special thanks to Dr. Charles Thangaraj and Dr. Arati Pati
# Code borrowed and adapted from https://github.com/ChristopherL891123/Research--Code-Development-for-Flow-Simulation
# Code borrowed uses Length of artery in meters and Radius of artery is in centimeters

# The units of viscosity are Pa · s = (N/m2) · s or Poise (dynes · s/cm2) with 1 Pa · s = 10 Poise
#normal blood viscosity is η = (3 ∼ 4) · 10−3 Pa · s
#v = (2.8 to 3.8)*(10)^-2 m^2/s
#mmhg t pa s
#sigma = v_nom/6
import LU  # borrowed
import statistics as stats
import MatrixGeneration  # borrowed
import matplotlib.pyplot as plt
import random
import threading
import pprint


print("*=*=*=*=* INITIALIZING VALUES *=*=*=*=*")
n = 10000 # size of matrix
A = MatrixGeneration.GENERATE(n) # generate the matrix
print("*=*=*=*=* INITIALIZED A MATRIX *=*=*=*=*")
U, L = LU.DECOMP(A, n, False) # decompose it
print("*=*=*=*=* INITIALIZED U , L MATRICES *=*=*=*=*")

pressure_list = []  # keep track of normally distributed random values for pressure
viscosity_list = []  # keep track of normally distributed random values for viscosity
rel_error_list = [0 for i in range(n)] # list of placeholders for average relative error for each run. Will be made up of nested lists where each nested list is the average for each run
abs_error_list = [0 for i in range(n)] # list of placeholders for average absolute error for each run. Will be made up of nested lists where each nested list is the average for each run
output = [" " for i in range(n)]
repeat = [] # used in cont()
progress_bar = ""
event = threading.Event()
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
        pressure_random_reading = random.gauss(0.32 / 1.5, 0.32 * 1.5)
        pressure_list.append(pressure_random_reading)
        viscosity_random_reading = random.gauss(0.035 / 1.5,0.035 * 1.5)  # viscosity nominal/1.5 , viscosity nominal * 1.5
        viscosity_list.append(viscosity_random_reading)

        avg_abs, avg_rel = engine(l=length, r=radius, Delta_P=pressure_random_reading, Nu=viscosity_random_reading)  # unit of l: meters, unit of r: centimeters

        temp_abs_err.append(avg_abs)
        temp_rel_err.append(avg_rel)

        # add row of values to table
        table += "|{: ^8} | {: ^30} | {: ^30} | {: ^30} | {: ^30}|".format(count, viscosity_random_reading, pressure_random_reading, avg_abs, avg_rel) + '\n'

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

    subsets_length = []
    subsets_radius = []

    tempLength = []  # where to store the temporary subset
    tempRadius = []
    cut = 0  # tell loop when to cut the flow of data, to append only two numbers into a subset.
    diff = 2

    for i in range(10):

        if cut == 2:
            subsets_length.append(tempLength)
            subsets_radius.append(tempRadius)
            tempLength = []
            tempRadius = []
            tempLength.append(length_range[i])
            tempRadius.append(radius_range[i])
            cut = 1
            # data is cut, but when the loop restarts it moves on with the next value, does not append index 2 as an operation was done at i = 2.
        else:
            tempLength.append(length_range[i])
            tempRadius.append(radius_range[i])
            cut += 1

    subsets_length.append(tempLength)  # append last values
    subsets_radius.append(tempRadius)  # append last values
    del tempLength
    del tempRadius

    # have 10 threads running at any given time; pc is at least 10 threads.
    thread_count = 0
    thread_table = {}
    for subsets_l in subsets_length:
        for length_readings in subsets_l:
            for subsets_r in subsets_radius:
                for radius_readings in subsets_r:
                    if thread_count < 10:  # if less than 10 treads running, run new thread
                        thread = threading.Thread(target=calc, name="thread{}".format(thread_count),
                                                  args=[output, temp_index_tables, write_index,write_index, length_readings,
                                                        radius_readings] )
                        thread.start() # start the thread asynchronously
                        thread_count += 1

                        thread_table.update({"thread{}".format(thread_count): [thread.name, thread.ident, thread.is_alive()]})

                        write_index += 1
                        temp_index_tables += 1

                        # if event.wait() == False: # wait until event.set() is called in one of the threads(thread finishes)
                        #     thread_count -= 1
                        #     print("*****CURRENT THREAD COUNT: ", thread_count, "*****")

                        # # debugging
                        # print("************************\n\n")
                        # pprint.pprint(thread_table)
                        # print("************************\n\n")
                        # # debugging

                        # debugging

                        # look for inactive threads
                        for thread in thread_table:
                            if False in thread_table[thread]:
                                del thread_table[thread]
                                thread_count -= 1
                        # debugging




    index = 0
    while True:

        if index == n:
            print(char for char in "FINISHED")
            exit(0)

        if output[index] != " ":

            file = open("outputs.txt", "a")
            file.write(output[index])
            file.close()
            print("wrote to file")
            plt.hist(rel_error_list[index]) # wants arrays
            plt.savefig("{a}rel.png".format(a=str(c2)),
                        bbox_inches="tight")
            plt.hist(abs_error_list[index])
            plt.savefig("{b}abs.png".format(b=str(c2)),
                        bbox_inches="tight")  # https://stackoverflow.com/questions/9622163/save-plot-to-image-file-instead-of-displaying-it-using-matplotlib
            plt.clf()  # https://www.tutorialspoint.com/how-do-i-close-all-the-open-pyplot-windows-matplotlib

            c2 += 1
            index += 1

        else:
            continue

main()
