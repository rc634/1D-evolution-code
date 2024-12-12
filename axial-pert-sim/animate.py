
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
U_file = 'U.dat'  # Replace 'phi.dat' with the path to your phi.dat file
H_file = 'H.dat'  # Replace 'psi.dat' with the path to your psi.dat file
x_file = 'x.dat'      # Replace 'x.dat' with the path to your x.dat file
v_UU_file = 'V_UU.dat'
v_UH_file = 'V_UH.dat'
v_HH_file = 'V_HH.dat'

# Assuming each row of the .dat files contains comma-separated values
U_data = np.loadtxt(U_file, delimiter=',')
H_data = np.loadtxt(H_file, delimiter=',')
x_data = np.loadtxt(x_file, delimiter=',')
v_UU_data = np.loadtxt(v_UU_file, delimiter=',')
v_UH_data = np.loadtxt(v_UH_file, delimiter=',')
v_HH_data = np.loadtxt(v_HH_file, delimiter=',')

# Create a figure and axis for plotting
fig, ax = plt.subplots()
lines = []
lines += [ax.plot([], [], lw=0.5, color='k')[0]]
lines += [ax.plot([], [], lw=2, label='V_UU', color='coral')[0]]
lines += [ax.plot([], [], lw=2, label='V_UH', color='limegreen')[0]]
lines += [ax.plot([], [], lw=2, label='V_HH', color='darkturquoise')[0]]
lines += [ax.plot([], [], lw=2, label=r'$U(x)$', color='red')[0]]
lines += [ax.plot([], [], lw=2, label=r'$H(x)$', color='blue')[0]]

# Add legend
ax.legend()

# Add title
plt.title("Animation of Gaussian Wavepacket")

# Label x-axis
plt.xlabel('X coordinate')

# Function to initialize the plot
def init():
    ax.set_xlim(np.min(x_data), np.max(x_data))  # Set the x-axis limit to the maximum length of both data sets
    ax.set_ylim(min(np.min(U_data), np.min(H_data)), max(np.max(U_data), np.max(H_data)))  # Set the y-axis limit to encompass both data sets
    return lines

# Function to update the plot for each frame of the animation
def update(frame):
    lines[0].set_data([x_data[frame][0],x_data[frame][-1]],[0.,0.] ) # grey line on y=0
    lines[1].set_data(x_data[frame], v_UU_data[frame])
    lines[2].set_data(x_data[frame], v_UH_data[frame])
    lines[3].set_data(x_data[frame], v_HH_data[frame])
    lines[4].set_data(x_data[frame], U_data[frame])  # Update psi data
    lines[5].set_data(x_data[frame], H_data[frame])  # Update phi data
    # Save each frame as an image
    # plt.savefig(f'../png/emdbh_tort_{frame:04d}.png', dpi=100)  # Adjust dpi as needed
    return lines

# Create the animation
ani = FuncAnimation(fig, update, frames=min(len(U_data), len(H_data)), init_func=init, blit=True, interval=40)

# Show the animation
plt.show()
