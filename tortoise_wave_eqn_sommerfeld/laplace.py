import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define the test function to generate time series data
def laplace_integrand(t,sig,om):
    return np.exp(-sig*t) * (np.cos(-om*t) + 1.j*np.sin(-om*t))   # Modify this function as needed

def laplace_transform(time, phi, dt, sig_, om_):
    integral = 0.
    for i in range(len(time)):
        t = time[i]
        integral += phi[i] * laplace_integrand(t,sig_,om_) * dt
    return np.abs(integral)

# resolution 
num_sig = 80
num_om = 80
num_time = 2000
end_time = 40

# QNM input
sigma = 1.41
omega = 2.78

# Generate time series data
time = np.linspace(0, end_time, num_time)
dt = time[1] - time[0]


phi = np.zeros(num_time)
for i in range(num_time):
    t = time[i]
    phi[i] = 2.2 * np.cos(omega*t + 1.1) * np.exp(-sigma*t)

# Generate 2D grid
x = np.linspace(-0., 3., num_sig)
y = np.linspace(-0., 3., num_om)
X, Y = np.meshgrid(x, y)
Z=np.zeros((num_sig, num_om))


plt.plot(time, phi)
plt.grid(True)
plt.show()

for i in range(num_sig):
    for j in range(num_om):
        sig_ = x[i]
        om_ = y[j]
        Z[i][j] = laplace_transform(time, phi, dt, sig_, om_)


#3d height plot

# Plot the 2D function# Add a color bar which maps values to colors
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, Z, cmap='viridis')

# Add labels and title
ax.set_xlabel(r'$\sigma$')
ax.set_ylabel(r'$\omega$')
ax.set_zlabel(r'Laplace of $f(t)$ : $\Lambda(\sigma, \omega)$')
plt.title('Ringdown laplace transform')

fig.colorbar(surf, label='Function Value')


plt.show()


for x in range(10):
    pass


# flat heatmap


# Plot the 2D function as a flat plot
plt.figure(figsize=(8, 6))
plt.imshow(Z, extent=[np.min(X), np.max(X), np.min(Y), np.max(Y)], cmap='viridis', origin='lower')
plt.colorbar(label=r'Laplace of $f(t)$ : $\Lambda(\sigma, \omega)$')
plt.xlabel(r'$\sigma$')
plt.ylabel(r'$\omega$')
plt.title('Ringdown laplace transform')
plt.show()

