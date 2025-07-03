import numpy as np
from amuse.units import nbody_system, units
from amuse.community import ph4
from amuse.lab import Particles
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Masses
masses = np.array([0.87, 0.8, 1.0]) | nbody_system.mass

# Initial positions 
x1 = -0.1855174644
positions = np.array([
    [x1, 0, 0],           # Body 1
    [1, 0, 0],          # Body 2
    [0, 0, 0]          # Body 3
]) | nbody_system.length

# Initial velocities (from symmetry and table)
v1 = 2.0221546880
v2 = 0.3968976468
velocities = np.array([
    [0, v1, 0],           # Body 1
    [0, v2, 0],           # Body 2
    [0, -(masses[0]*v1 + masses[1]*v2)/masses[2], 0] # Body 3 (momentum conservation)
]) | nbody_system.speed

# Create particles
bodies = Particles(3)
bodies.mass = masses
bodies.position = positions
bodies.velocity = velocities

# Set up the integrator
gravity = ph4.Ph4()
gravity.particles.add_particles(bodies)
channel = gravity.particles.new_channel_to(bodies)

# Evolve and store positions
dt = 0.01 | nbody_system.time
x, y = np.zeros((2001, 3)), np.zeros((2001, 3))
for i in range(2001):
    gravity.evolve_model(i*dt)
    x[i] = gravity.particles.x.value_in(nbody_system.length)
    y[i] = gravity.particles.y.value_in(nbody_system.length)
    channel.copy()
gravity.stop()

# Creating a plot
plt.figure(figsize=(6,6))
plt.plot(x[:,0], y[:,0], 'b', label='Body 1')
plt.plot(x[:,1], y[:,1], 'r', label='Body 2')
plt.plot(x[:,2], y[:,2], 'k', label='Body 3')
plt.xlabel('x')
plt.ylabel('y')
plt.title("Plot of the used braid")
plt.legend()
plt.axis('equal')
plt.savefig("braid.pdf")  # Save each plot
plt.close()  # Close the plot to save memory###