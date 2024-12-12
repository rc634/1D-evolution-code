import matplotlib.pyplot as plt
import numpy as np
import csv

def read_data(file_name):
    data = []
    with open(file_name, 'r') as file:
        file = csv.reader(file,delimiter=',')
        for line in file:
            for number in line:
            # Assuming data is space-separated
                data.append(list(map(float, number.strip().split())))
    return data


# Read data from files
data_x = read_data("x.dat")
data_Vuu = read_data("V_UU.dat")
data_Vuh = read_data("V_UH.dat")
data_Vhh = read_data("V_HH.dat")

# Extract x and y data
x = [row[0] for row in data_x]
Vuu = [row[0] for row in data_Vuu]
Vuh = [row[0] for row in data_Vuh]
Vhh = [row[0] for row in data_Vhh]

plt.plot(x, Vuu, label=r'$V_{UU}$', color='red')
plt.plot(x, Vuh, label=r'$V_{UH}$', color='green')
plt.plot(x, Vhh, label=r'$V_{HH}$', color='blue')



plt.xlabel('x')
#plt.ylabel('a or c')
# plt.title(r'a06_Q10130NEW')
plt.legend()
plt.xlim([40,185])
# #plt.ylim([-40,0]) # was minus 62
plt.grid(True)
# plt.savefig("a06_Q10130NEW.png")
# #plt.savefig("phi-psi_r=60.png")
plt.show()
