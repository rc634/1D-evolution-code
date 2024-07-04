
# import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation
# import numpy as np

# # Read the .dat files and extract rows of data
# phi_file = 'phi.dat'  # Replace 'phi.dat' with the path to your phi.dat file
# psi_file = 'psi.dat'  # Replace 'psi.dat' with the path to your psi.dat file
# x_file = 'x.dat'      # Replace 'x.dat' with the path to your x.dat file

# # Assuming each row of the .dat files contains comma-separated values
# phi_data = np.loadtxt(phi_file, delimiter=',')
# psi_data = np.loadtxt(psi_file, delimiter=',')
# xdata = np.loadtxt(x_file, delimiter=',')

# # Create a figure and axis for plotting
# fig, ax = plt.subplots()
# lines = [ax.plot([], [], lw=2, label='phi.dat')[0], ax.plot([], [], lw=2, label='psi.dat')[0]]  # Create two lines for phi and psi data

# # Add legend
# ax.legend()

# # Add title
# plt.title("Animation of Gaussian Wavepacket")

# # Label x-axis
# plt.xlabel('Spatial Gridpoints')

# # Function to initialize the plot
# def init():
#     ax.set_xlim(min(x_data), max(x_data))  # Set the x-axis limit to the maximum length of both data sets
#     ax.set_ylim(min(np.min(phi_data), np.min(psi_data)), max(np.max(phi_data), np.max(psi_data)))  # Set the y-axis limit to encompass both data sets
#     return lines

# # Function to update the plot for each frame of the animation
# def update(frame):
#     lines[0].set_data(range(len(phi_data[frame])), phi_data[frame])  # Update phi data
#     lines[1].set_data(range(len(psi_data[frame])), psi_data[frame])  # Update psi data
#     return lines

# # Create the animation
# ani = FuncAnimation(fig, update, frames=min(len(phi_data), len(psi_data)), init_func=init, blit=True, interval=40)

# # Show the animation
# plt.show()



import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# Read the .dat files and extract rows of data
phi_file = 'phi.dat'  # Replace 'phi.dat' with the path to your phi.dat file
psi_file = 'psi.dat'  # Replace 'psi.dat' with the path to your psi.dat file
x_file = 'x.dat'      # Replace 'x.dat' with the path to your x.dat file

# Assuming each row of the .dat files contains comma-separated values
phi_data = np.loadtxt(phi_file, delimiter=',')
psi_data = np.loadtxt(psi_file, delimiter=',')
x_data = np.loadtxt(x_file, delimiter=',')

# Create a figure and axis for plotting
fig, ax = plt.subplots()
lines = [ax.plot([], [], lw=2, label='phi.dat')[0], ax.plot([], [], lw=2, label='psi.dat')[0]]  # Create two lines for phi and psi data

# Add legend
ax.legend()

# Add title
plt.title("Animation of Gaussian Wavepacket")

# Label x-axis
plt.xlabel('X coordinate')

# Function to initialize the plot
def init():
    ax.set_xlim(np.min(x_data), np.max(x_data))  # Set the x-axis limit to the maximum length of both data sets
    ax.set_ylim(min(np.min(phi_data), np.min(psi_data)), max(np.max(phi_data), np.max(psi_data)))  # Set the y-axis limit to encompass both data sets
    return lines

# Function to update the plot for each frame of the animation
def update(frame):
    lines[0].set_data(x_data[frame], phi_data[frame])  # Update phi data
    lines[1].set_data(x_data[frame], psi_data[frame])  # Update psi data
    return lines

# Create the animation
ani = FuncAnimation(fig, update, frames=min(len(phi_data), len(psi_data)), init_func=init, blit=True, interval=40)

# Show the animation
plt.show()


