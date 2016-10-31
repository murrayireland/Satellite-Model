# Executable for orbital and attitude dynamics simulation
# Murray Ireland
# October 2016


# Load packages
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm

# Load data
from OrbitModel import OrbitModel
import UniversalProps, EarthProps, SatProps

# Set parameters of orbit
OrbitalElements = {'Semi-major axis' : 6978.1e3,
    'Eccentricity' : 0.0,
    'Inclination' : 97.7,
    'Right ascension' : 0,
    'Argument of perigee' : 0}

# Number of orbits to simulate
NumOrbits = 4

# Initialise simulation object
Sim = OrbitModel( OrbitalElements, UniversalProps, EarthProps, SatProps, NumOrbits )

# Preallocate save structures
NSamps = Sim.RunTime/Sim.SampleTime
NSamps = int( np.floor( NSamps ) + 1 )
NStates = Sim.States.shape
Time = np.zeros(NSamps)
States = np.zeros(( NStates[0], NSamps))

# Simulation loop
k = 0
t = Sim.Time
ts = Sim.Time
NumSteps = int( Sim.RunTime/Sim.TimeStep ) + 1
SampFreq = int( Sim.SampleTime/Sim.TimeStep )
States[:,0] = Sim.States.A1

pbar = tqdm( total = Sim.RunTime )

for i in range(1,NumSteps) :

    # Update model
    Sim.UpdateModel()

    # Save data
    if ( i % SampFreq == 0 ) :
        k = k + 1
        Time[k] = Sim.Time
        States[:,k] = Sim.States.A1

    # Update progress bar
    pbar.update( Sim.TimeStep )

pbar.close()

fig = plt.figure( dpi=267, figsize=[5,5], facecolor='black' )
ax = fig.add_subplot( 111, projection='3d', axisbg='black' )
plt.plot( States[0,:], States[1,:], States[2,:], 'red' )
plt.axis( 'equal' )
plt.axis( 'off' )
plt.show()
