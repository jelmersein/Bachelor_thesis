import numpy as np
from amuse.units import nbody_system, units
from amuse.community import ph4
from amuse.lab import Particles
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Masses
masses = np.array([0.87, 0.8, 1.0, 1.0]) | nbody_system.mass

# Initial positions 
x1 = -0.1855174644
positions = np.array([
    [x1, 0, 0],           # Body 1
    [1, 0, 0],          # Body 2
    [0, 0, 0],             # Body 3
    [-6,-6,0]            # Incoming Body
]) | nbody_system.length

# Initial velocities 
v1 = 2.0221546880
v2 = 0.3968976468
velocities = np.array([
    [0, v1, 0],           # Body 1
    [0, v2, 0],           # Body 2
    [0, -(masses[0]*v1 + masses[1]*v2)/masses[2], 0],  # Body 3 (momentum conservation)
    [2,2,0]           # Incoming Body
]) | nbody_system.speed

# Create particles
bodies = Particles(4)
bodies.mass = masses
bodies.position = positions
bodies.velocity = velocities

# Set up the integrator
gravity = ph4.Ph4()
gravity.particles.add_particles(bodies)
channel = gravity.particles.new_channel_to(bodies)
# Evolve and store positions
dt = 0.01 | nbody_system.time
x, y = [], []
for i in range(601):
    gravity.evolve_model(i*dt)
    x.append(gravity.particles.x.value_in(nbody_system.length))
    y.append(gravity.particles.y.value_in(nbody_system.length))
    channel.copy()
gravity.stop()

x = np.array(x)
y = np.array(y)

# Plot
plt.figure(figsize=(6,6))
plt.plot(x[:,0], y[:,0], 'b', label='Body 1')
plt.plot(x[:,1], y[:,1], 'r', label='Body 2')
plt.plot(x[:,2], y[:,2], 'k', label='Body 3')
plt.plot(x[:,3], y[:,3], 'y', label='Incoming Body')
plt.xlabel('x')
plt.ylabel('y')
plt.title("Example of Single-Triple Scattering")
plt.legend()
plt.axis('equal')
plt.savefig("scatter.pdf")  # Save each plot
plt.close()  # Close the plot to save memory###


# Optional animation code 

#fig, ax = plt.subplots(figsize=(6,6))
#lines = [
#    ax.plot([], [], color, label=label)[0]
#    for color, label in zip(['b', 'r', 'k', 'y'], ['Body 1', 'Body 2', 'Body 3', 'Incoming Body'])
#]
#ax.set_xlim(np.min(x)-1, np.max(x)+1)
#ax.set_ylim(np.min(y)-1, np.max(y)+1)
#ax.set_title('Example of Single-Triple Scattering')
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.legend()
#ax.axis('equal')
#
#def init():
#    for line in lines:
#        line.set_data([], [])
#    return lines
#
#def animate(i):
#    for j, line in enumerate(lines):
#        line.set_data(x[:i, j], y[:i, j])
#    return lines
#
#ani = animation.FuncAnimation(
#    fig, animate, frames=len(x), init_func=init, blit=True, interval=10
#)
#
#ani.save('single_triple_scatter.mp4', writer='ffmpeg', fps=60)

# Initializing reversed scattering
bodies.velocity = -bodies.velocity
gravity = ph4.Ph4()
gravity.particles.add_particles(bodies)
channel = gravity.particles.new_channel_to(bodies)
#Evolve and store positions again
x, y = [], []
for i in range(2001):
    gravity.evolve_model(i*dt)
    x.append(gravity.particles.x.value_in(nbody_system.length))
    y.append(gravity.particles.y.value_in(nbody_system.length))
    channel.copy()
gravity.stop()
x = np.array(x)
y = np.array(y)
# Plot again
plt.figure(figsize=(6,6))
plt.plot(x[:,0], y[:,0], 'b', label='Body 1')
plt.plot(x[:,1], y[:,1], 'r', label='Body 2')
plt.plot(x[:,2], y[:,2], 'k', label='Body 3')
plt.plot(x[:,3], y[:,3], 'y', label='Outgoing Body')
plt.xlabel('x')
plt.ylabel('y')
plt.title("Reversed Single-Triple Scattering")
plt.legend()
plt.xlim(-4, 4)
plt.ylim(-4, 4)
plt.savefig("bbinteraction.pdf")  # Save each plot
plt.close()  


# Zoomed in plot to see the braid-like product
plt.figure(figsize=(6,6))
plt.plot(x[:,0], y[:,0], 'b', label='Body 1')
plt.plot(x[:,1], y[:,1], 'r', label='Body 2')
plt.plot(x[:,2], y[:,2], 'k', label='Body 3')
plt.plot(x[:,3], y[:,3], 'y', label='Incoming Body')
plt.xlabel('x')
plt.ylabel('y')
plt.title("Reversed Single-Triple Scattering (Zoomed In)")
plt.legend()
plt.xlim(-1, 1.5)
plt.ylim(-1.25, 1.25)
plt.savefig("bbzoomed.pdf")  # Save each plot
plt.close()  # Close the plot to save memory###

# Optional animation code again

#fig, ax = plt.subplots(figsize=(6,6))
#lines = [
#    ax.plot([], [], color, label=label)[0]
#    for color, label in zip(['b', 'r', 'k', 'y'], ['Body 1', 'Body 2', 'Body 3', 'Incoming Body'])
#]
#ax.set_xlim(np.min(x)-1, np.max(x)+1)
#ax.set_ylim(np.min(y)-1, np.max(y)+1)
#ax.set_title('Binary-Binary scattering yields Braid')
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.legend()
#ax.axis('equal')
#
#def init():
#    for line in lines:
#        line.set_data([], [])
#    return lines
#
#def animate(i):
#    for j, line in enumerate(lines):
#        line.set_data(x[:i, j], y[:i, j])
#    return lines
#
#ani = animation.FuncAnimation(
#    fig, animate, frames=len(x), init_func=init, blit=True, interval=10
#)
#
#ani.save('Get_braid.mp4', writer='ffmpeg', fps=60)