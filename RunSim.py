# Executable for orbital and attitude dynamics simulation
# Murray Ireland
# October 2016


# Load packages
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.basemap import Basemap
from tqdm import tqdm

# Load data
from OrbitModel import OrbitModel
import UniversalProps, EarthProps, SatProps
import Data, Plots
from Presets import PresetOrbits

# Set parameters of orbit
OrbitalElements = {'Semi-major axis' : 26562.1e3,
    'Eccentricity' : 0.74,
    'Inclination' : 0,
    'Right ascension' : 0,
    'Argument of perigee' : -90}

# Or load orbit presets from file
PresetOrbit = 'Molniya'
OrbitalElements = PresetOrbits[PresetOrbit]

# Number of orbits to simulate
NumOrbits = 1

# Initialise simulation object
Sim = OrbitModel( OrbitalElements, UniversalProps, EarthProps, SatProps, NumOrbits )

# Preallocate save structures
NSamps = Sim.RunTime/Sim.SampleTime
NSamps = int( np.floor( NSamps ) + 1 )
NStates = Sim.States.shape
Data.Time = np.zeros(NSamps)
Data.States = np.zeros(( NStates[0], NSamps ))
Data.StatesECEF = np.zeros(( 6, NSamps ))
Data.LatLongAlt = np.zeros(( 3, NSamps ))

# Simulation loop
k = 0
NumSteps = int( Sim.RunTime/Sim.TimeStep ) + 1
SampFreq = int( Sim.SampleTime/Sim.TimeStep )
Data.Time[0] = Sim.Time
Data.States[:,0] = Sim.States
Data.StatesECEF[:,0] = Sim.StatesECEF
Data.LatLongAlt[:,0] = Sim.LatLongAlt

pbar = tqdm( total = Sim.RunTime )

for i in xrange(1,NumSteps) :

    # Update model
    Sim.UpdateModel()

    # Save data
    if ( i % SampFreq == 0 ) :
        k = k + 1
        Data.Time[k] = Sim.Time
        Data.States[:,k] = Sim.States
        Data.StatesECEF[:,k] = Sim.StatesECEF
        Data.LatLongAlt[:,k] = Sim.LatLongAlt

    # Update progress bar
    pbar.update( Sim.TimeStep )

pbar.close()

print(Data.States.shape)

# Concatenate Data
# Labels = ['Time', 'x', 'y', 'z', 'u', 'v', 'w', 'x_ecef', 'y_ecef', 'z_ecef', 'u_ecef', 'v_ecef', 'w_ecef', 'Lat', 'Long', 'Alt' ]
# AllData = np.column_stack( (Data.Time, Data.States.T, Data.StatesECEF.T, Data.LatLongAlt.T ) )

# Save Data
# np.savetxt( 'Data.csv', AllData, delimiter=',' )

# Plots
Plots.Plot3DPath( Data )
Plots.PlotGroundTrack( Data )
Plots.PlotAttitude( Data )
plt.show()