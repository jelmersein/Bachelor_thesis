import numpy as np
from amuse.units import nbody_system, units
from amuse.community import ph4
from amuse.lab import Particles
import matplotlib.pyplot as plt

# Function to check if a pair of bodies is bound
def check_bound(particles):
    centered_momentum_velocity = particles.velocity - particles.center_of_mass_velocity()
    kinetic_energy = (0.5 * particles.mass * (centered_momentum_velocity.lengths() ** 2).sum()).sum()
    potential_energy = -(2*particles.mass[0]*particles.mass[1] / ((particles.x[0] - particles.x[1])**2 + (particles.y[0] - particles.y[1])**2).sqrt())*nbody_system.G
    total_energy = kinetic_energy + potential_energy
    return total_energy < 0 |nbody_system.length**2 * nbody_system.time**-2 * nbody_system.mass

# Masses
masses = np.array([0.87, 0.8, 1.0, 1.0]) | nbody_system.mass

# For loop of initial positions incoming body
double_binaries, not_bb = [], []
for phi in np.union1d(
    np.union1d(np.union1d(np.union1d(np.linspace(0, 2*np.pi*0.097, 98),
    np.linspace(0.099*2*np.pi, 2*np.pi*0.210, 112)),
    np.linspace(0.212*2*np.pi, 2*np.pi*0.347, 136)),
    np.union1d(np.linspace(0.349*2*np.pi, 2*np.pi*0.522, 174),
    np.linspace([0.524*2*np.pi], 2*np.pi*0.526, 3))),
    np.union1d(np.union1d(np.union1d(np.linspace([0.528*2*np.pi], 2*np.pi*0.592, 65),
    np.linspace([0.594*2*np.pi], 2*np.pi*0.708, 115)),
    np.union1d(np.linspace([0.710*2*np.pi], 2*np.pi*0.781, 72),
    np.linspace([0.783*2*np.pi], 2*np.pi*0.982, 200))),
    np.linspace([0.984*2*np.pi], 2*np.pi*0.999, 16))): #removing the stalling cases
    print("phi:", phi)
    x_inc = 6*np.sqrt(2) * np.cos(phi)
    y_inc = 6*np.sqrt(2) * np.sin(phi)

    # Initial positions 
    x1 = -0.1855174644
    positions = np.array([
        [x1, 0, 0],           # Body 1
        [1, 0, 0],          # Body 2
        [0, 0, 0],             # Body 3
        [x_inc,y_inc,0]            # Incoming Body
    ]) | nbody_system.length

    # Initial velocities 
    v1 = 2.0221546880
    v2 = 0.3968976468
    vx_inc = -2*np.sqrt(2) * np.cos(phi)
    vy_inc = -2*np.sqrt(2) * np.sin(phi)
    velocities = np.array([
        [0, v1, 0],           # Body 1
        [0, v2, 0],           # Body 2
        [0, -(masses[0]*v1 + masses[1]*v2)/masses[2], 0],  # Body 3 (momentum conservation)
        [vx_inc,vy_inc,0]           # Incoming Body
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
    x, y = np.zeros((2001, 4)), np.zeros((2001, 4))
    for i in range(2001):
        gravity.evolve_model(i*dt)
        #x[i] = gravity.particles.x.value_in(nbody_system.length)
        #y[i] = gravity.particles.y.value_in(nbody_system.length)
        channel.copy()
    gravity.stop()

    # Check if the bodies form a double binary
    if check_bound(bodies[[0,1]]) and check_bound(bodies[[2,3]]):
        double_binaries.append((x_inc, y_inc))
    elif check_bound(bodies[[0,2]]) and check_bound(bodies[[1,3]]):
        double_binaries.append((x_inc, y_inc))
    elif check_bound(bodies[[0,3]]) and check_bound(bodies[[1,2]]):
        double_binaries.append((x_inc, y_inc))
    else:
        not_bb.append((x_inc, y_inc))
    # Plot individual cases (optional)
    # Uncomment the following lines to visualize each case
    #plt.figure(figsize=(6,6))
    #plt.plot(x[:,0], y[:,0], 'b', label='Body 1')
    #plt.plot(x[:,1], y[:,1], 'r', label='Body 2')
    #plt.plot(x[:,2], y[:,2], 'k', label='Body 3')
    #plt.plot(x[:,3], y[:,3], 'y', label='Incoming Body')
    #plt.xlabel('x')
    #plt.ylabel('y')
    #plt.title(f"Incoming Body Angle: ({phi:.2f})")
    #plt.legend()
    #plt.axis('equal')
    #plt.savefig(f"Documents/BRP/forloop/plot_{phi:.2f}.png")  # Save each plot
    #plt.close()  # Close the plot to save memory###

    print("done with phi:", phi)

# Preparing arrays for plotting
double_binaries = np.array(double_binaries)
not_bb = np.array(not_bb)

# Print the amount of double binaries found
print("Double binaries found:", len(double_binaries))

# Plotting the results
plt.figure(figsize=(6, 6))
plt.scatter(not_bb[:, 0], not_bb[:, 1], color='r', label='Other Scattering products', s=2)
plt.scatter(double_binaries[:,0], double_binaries[:, 1],marker='+', color='b',  label='Double Binaries', s=20)
plt.title('Double Binaries as scattering product')
plt.xlabel('x position of incoming body')
plt.ylabel('y position of incoming body')
plt.legend()
plt.axis('equal')
plt.savefig("double_binaries_results.pdf")
plt.show()

