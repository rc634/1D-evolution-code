import numpy as np
import matplotlib.pyplot as plt

# Load data from field.dat and time.dat
field_data = np.fromfile('phi_of_time.dat', dtype=np.float64, count=-1, sep=' ')
time_data = np.fromfile('time.dat', dtype=np.float64, count=-1, sep=' ')


nums = 10000
end_t = 400
dt = float(end_t)/float(nums-1)
A = 0.8
B = 0.5
w1 = 3.47
w2 = 4.54
wtau = 0.2

test_phi = np.zeros(nums)
test_time = np.zeros(nums)

for i in range(nums):
    t = i * dt
    test_phi[i] = (A * np.sin(w1 * t) + B * np.cos(w2 * t)) * np.exp(-wtau * t)
    test_time[i] = t

# Calculate Fourier transform
freqs = np.fft.fftfreq(len(time_data), time_data[1] - time_data[0])
field_fft = np.fft.fft(field_data)

# Calculate Fourier transform
test_freqs = np.fft.fftfreq(len(test_time), test_time[1] - test_time[0])
test_field_fft = np.fft.fft(test_phi)

print(type(field_fft[0]))

# Plot the Fourier-transformed data

#simulation data
plt.plot(2.*np.pi*freqs, (field_fft))

#test function data
# plt.plot(2.*np.pi*test_freqs, np.abs(test_field_fft))

plt.xlabel('Frequency')
plt.ylabel('Amplitude')
plt.title('Fourier Transform of Field Data')
plt.grid(True)
plt.show()

