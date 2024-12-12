import matplotlib.pyplot as plt
import numpy as np

def read_data(file_name):
    data = []
    with open(file_name, 'r') as file:
        for line in file:
            # Assuming data is space-separated
            data.append(list(map(float, line.strip().split())))
    return data


# Read data from files
data_a = read_data("U_of_time.dat")
data_b = read_data("H_of_time.dat")
data_c = read_data("time.dat")
data_d = read_data("log_U_of_time.dat")
data_e = read_data("log_H_of_time.dat")
#data_hq =read_data('../DATA/log_phi_of_time.dat')
#data_hqt=read_data("../DATA/time.dat")
# data_25 =read_data('../DATA/log_phi_of_time_25.dat')
# data_time_25=read_data("../DATA/time_25.dat")

# Extract x and y data
phi = [row[0] for row in data_a]
psi = [row[0] for row in data_b]
time = [row[0] for row in data_c]
log_phi = [row[0] for row in data_d]
log_psi = [row[0] for row in data_e]

#log_phi_hq = [row[0] for row in data_hq]
#time_hq = [row[0] for row in data_hqt]
# log_phi_25 = [row[0] for row in data_25]
# time_25 = [row[0] for row in data_time_25]

# SC
wr = 0.110455
wi = -0.104896
A = 30

# # RN
# wr = 0.133250
# wi = - 0.101575
A = 6

fit_data = np.zeros(len(time))
for i in range(len(time)):
    fit_data[i] = 0.5 * np.log( np.power(A * np.sin(wr * time[i]),2) ) + wi * time[i]

# Plotting

#plt.plot(time, phi, label=r'$\phi$', color='red')
#plt.plot(time, psi, label=r'$\psi$', color='blue')

plt.plot(time, log_phi, label=r'$\log\|\phi\|$', color='red')
plt.plot(time, log_psi, label=r'$\log\|\psi\|$', color='blue')


plt.xlabel('Time')
#plt.ylabel('a or c')
plt.title(r'a06_Q10130NEW')
plt.legend()
#plt.xlim([time[0],1900])
#plt.ylim([-40,0]) # was minus 62
plt.grid(True)
plt.savefig("a06_Q10130NEW.png")
#plt.savefig("phi-psi_r=60.png")
plt.show()
