import matplotlib.pyplot as plt

mean = []
stdev = []
max = []
median = []
min = []

with open("stats_outputs_numpy.txt", 'r') as f:
    while True:
        line = f.readline()

        if "MEAN" in line:
            l = line.split(":")
            mean.append(float(l[-1]))

        if "STANDARD DEVIATION" in line:
            l = line.split(":")
            stdev.append(float(l[-1]))

        if "MAX" in line:
            l = line.split(":")
            max.append(float(l[-1]))

        if "MEDIAN" in line:
            l = line.split(":")
            median.append(float(l[-1]))

        if "MIN" in line:
            l = line.split(":")
            min.append(float(l[-1]))

        if "END OF FILE" in line:
            break



plt.hist(mean, bins=20)
plt.title("Mean")
plt.savefig("mean.png")
plt.show()

plt.hist(stdev, bins=20)
plt.title("Standard Deviation")
plt.savefig("stdev.png")
plt.show()

plt.hist(median, bins=20)
plt.title("Median")
plt.savefig("median.png")
plt.show()

plt.hist(max, bins=20)
plt.title("Max")
plt.savefig("max.png")
plt.show()

plt.hist(min, bins=20)
plt.title("Min")
plt.savefig("min.png")
plt.show()
