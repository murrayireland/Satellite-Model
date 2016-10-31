# Class definition for orbital and attitude dynamics model
# Murray Ireland
# October 2016

from math import sin, cos, tan, atan2, sqrt, pi
import numpy as np

class OrbitModel :

      Solver = 'RK4'

      def __init__( self, Elements, UniversalProps, EarthProps, SatProps, NumOrbits ) :
            # Class constructor
            
            # Assign physical properties
            self.OrbitalElements = Elements
            self.UniversalProps = UniversalProps
            self.EarthProps = EarthProps
            self.SatProps = SatProps

            # Initialise sim variables
            self.Time = 0
            self.TimeStep = 1
            self.SampleTime = 1
            self.States = self.Initialise()
            self.Inputs = 0
            self.RunTime = NumOrbits*self.OrbitalElements['Orbital period']
      
      def Initialise( self ) :
            # Obtain initial position and velocity from orbital Elements

            # Orbital elements
            a = self.OrbitalElements['Semi-major axis']
            e = self.OrbitalElements['Eccentricity']
            i = self.OrbitalElements['Inclination']
            Omega = self.OrbitalElements['Argument of perigee']
            omega = self.OrbitalElements['Right ascension']

            # Universal properties
            G = self.UniversalProps.GravitationalConstant

            # Earth properties
            M = self.EarthProps.Mass

            # Satellite properties
            m = self.SatProps.Mass

            # Gravitational parameter
            mu = G*(M + m)    

            # Mean angular motion
            n = sqrt(mu/a**3)

            # True anomaly
            TA = 0*pi/180
                  
            # Eccentric anomaly
            EA = atan2(sqrt(1 - e**2)*sin(TA), e + cos(TA))

            # Position in intermediate frame
            rQ = np.array([[a*(cos(EA) - e)],
                        [a*sqrt(1 - e**2)*sin(EA)],
                        [0]])
            vf = n*a/(1 - e*cos(EA))
            vQ = vf*np.array([[-sin(EA)],
                              [sqrt(1 - e**2)*cos(EA)],
                              [0]])
            
            # Rotate to ECI frame
            ROmega = np.array([[cos(Omega), -sin(Omega), 0],
                               [sin(Omega),  cos(Omega), 0],
                               [0,                0,               1]])
            Ri = np.array([[1, 0,            0],
                           [0, cos(i), -sin(i)],
                           [0, sin(i),  cos(i)]])
            Romega = np.array([[cos(omega), -sin(omega), 0],
                              [sin(omega),  cos(omega), 0],
                              [0,                0,               1]])

            # Initial position and velocity in ECI frame
            R = np.mat(ROmega)*np.mat(Ri)*np.mat(Romega)
            Position = R*rQ
            Velocity = R*vQ
            States = np.concatenate( (Position, Velocity), axis=0 )

            # Save useful properties
            self.OrbitalElements['Gravitational parameter'] = mu
            self.OrbitalElements['Orbital period'] = 2*pi*sqrt(a**3/mu)

            # Outputs
            return States

      def UpdateModel( self ) :
            # Update model dynamics and guidance

            self.States, self.Time = self.Integrate( self.States, self.Inputs, self.Time )
      
      def OrbitalDynamics( self, X, U, t ) :
            # Orbital dynamics model

            # Gravitational acceleration
            r = np.linalg.norm(X[0:3])
            mu = self.OrbitalElements['Gravitational parameter']
            g = -mu*X[0:3]/r**3

            # State derivatives
            Xdot = np.concatenate( ( X[3:7], g ), axis=0 )
            return Xdot

      def Integrate( self, X, U, t ) :
            # Integrate dynamic model

            dt = self.TimeStep

            if self.Solver == 'Euler' :
                  Xdot = self.OrbitalDynamics( X, U, t )
                  X = X + dt*Xdot
            elif self.Solver == 'RK4' :
                  k1 = self.OrbitalDynamics( X, U, t )
                  k2 = self.OrbitalDynamics( X + dt*k1/2, U, t + dt/2 )
                  k3 = self.OrbitalDynamics( X + dt*k2/2, U, t + dt/2 )
                  k4 = self.OrbitalDynamics( X + dt*k3, U, t + dt )
                  X = X + dt*(k1 + 2*k2 + 2*k3 + k4)/6
            
            t = round( (t + dt)/dt, 0 )*dt

            return ( X, t )
