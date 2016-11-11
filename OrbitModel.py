# Class definition for orbital and attitude dynamics model
# Murray Ireland
# October 2016

from math import sin, cos, tan, atan2, sqrt, pi
import numpy as np
from numpy.linalg import inv

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
            self.TimeStep = 5
            self.SampleTime = 5
            self.States = self.Initialise()
            self.Inputs = 0
            self.RunTime = NumOrbits*self.OrbitalElements['Orbital period']
            #self.RunTime = 3
            self.StatesECEF = self.ECI2ECEF( self.States, self.Time )
            self.LatLongAlt = self.ECEF2LatLong( self.StatesECEF )      

      def UpdateModel( self ) :
            # Update model dynamics and guidance

            # Update states
            self.States, self.Time = self.Integrate( self.States, self.Inputs, self.Time )

            # Update ECEF position
            self.StatesECEF = self.ECI2ECEF( self.States, self.Time )

            # Update latitude and longitude
            self.LatLongAlt = self.ECEF2LatLong( self.StatesECEF )

      def Initialise( self ) :
            # Obtain initial position and velocity from orbital Elements

            # Orbital elements
            a = self.OrbitalElements['Semi-major axis']
            e = self.OrbitalElements['Eccentricity']
            i = self.OrbitalElements['Inclination']*pi/180
            omega = self.OrbitalElements['Argument of perigee']*pi/180
            Omega = self.OrbitalElements['Right ascension']*pi/180

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
            p = a*(1 - e**2)
            r = p/(1 + e*cos(TA))
            rQ = np.array([[r*cos(TA)],
                           [r*sin(TA)],
                           [0]])
            vf = sqrt(mu/p)
            vQ = vf*np.array([[-sin(TA)],
                              [e + cos(TA)],
                              [0]])
            
            # Rotate to ECI frame
            ROmega = np.array([[cos(Omega), -sin(Omega), 0],
                               [sin(Omega),  cos(Omega), 0],
                               [0,           0,          1]])
            Ri = np.array([[1, 0,       0],
                           [0, cos(i), -sin(i)],
                           [0, sin(i),  cos(i)]])
            Romega = np.array([[cos(omega), -sin(omega), 0],
                               [sin(omega),  cos(omega), 0],
                               [0,           0,          1]])

            # Initial position and velocity in ECI frame
            R = np.mat(ROmega)*np.mat(Ri)*np.mat(Romega)
            # R = R.T
            Position = R*rQ
            Velocity = R*vQ

            # Satellite properties
            Attitude = np.array([[1], [0], [0], [0]])
            AngVel = np.array([[0.0], [0.0], [0.0]])

            States = np.concatenate( (Position, Velocity, Attitude, AngVel), axis=0 )

            # Save useful properties
            self.OrbitalElements['Gravitational parameter'] = mu
            self.OrbitalElements['Orbital period'] = 2*pi*sqrt(a**3/mu)

            # Outputs
            return States.A1
      
      def Dynamics( self, X, U, t ) :
            # Orbital dynamics model

            # Convert arrays to column vectors
            X = np.matrix(X).T
            #U = np.matrix(U).T

            # Orbital dynamics
            Xdot_orb = self.OrbitalDynamics( X, t )

            # Attitude dynamics
            Xdot_att = self.SatelliteDynamics( X, U, t )

            # Combine states
            Xdot = np.concatenate( (Xdot_orb, Xdot_att), axis=0 )

            # State derivatives
            return Xdot.A1

      def OrbitalDynamics( self, X, t ) :
            # Orbital dynamics

            # Gravitational acceleration
            r = np.linalg.norm(X[0:3])
            mu = self.OrbitalElements['Gravitational parameter']
            g = -mu*X[0:3]/r**3

            # State derivatives
            Xdot = np.concatenate( ( X[3:6], g ), axis=0 )

            return Xdot

      def SatelliteDynamics( self, X, U, t ) :
            # Satellite attitude Dynamics

            # Satellite properties
            I = self.SatProps.Inertia
            mm = self.SatProps.MagneticMoment

            # Moment
            M = np.array([[0], [0], [0]])

            # Angular acceleration
            AngAcc = inv(np.mat(I))*(M - np.cross(X[10:13], I*X[10:13], axis=0))

            # Questernion rates
            RateMatrix = np.array([[0, -X[10], -X[11], -X[12]],
                                    [X[10], 0, X[12], -X[11]],
                                    [X[11], -X[12], 0, X[10]],
                                    [X[12], X[11], -X[10], 0]])
            
            Quaternions = 0.5*np.mat(RateMatrix)*X[6:10]

            # State derivatives
            Xdot = np.concatenate( (Quaternions, AngAcc), axis=0 )

            return Xdot

      def Integrate( self, X, U, t ) :
            # Integrate dynamic model

            dt = self.TimeStep

            if self.Solver == 'Euler' :
                  Xdot = self.Dynamics( X, U, t )
                  X = X + dt*Xdot
            elif self.Solver == 'RK4' :
                  k1 = self.Dynamics( X, U, t )
                  k2 = self.Dynamics( X + dt*k1/2, U, t + dt/2 )
                  k3 = self.Dynamics( X + dt*k2/2, U, t + dt/2 )
                  k4 = self.Dynamics( X + dt*k3, U, t + dt )
                  X = X + dt*(k1 + 2*k2 + 2*k3 + k4)/6
            
            t = round( (t + dt)/dt, 0 )*dt

            return ( X, t )

      def ECI2ECEF( self, XECI, t ) :
            # Convert properties in ECI frame to ECEF frame

            # Get rotational speed of Earth
            Omega = self.EarthProps.AngularRate

            # Get current rotation of Earth
            theta = Omega*t

            # Convert position and velocity arrays to column vectors
            rECI = np.matrix(XECI[0:3]).T
            vECI = np.matrix(XECI[3:6]).T
            
            # Rotate current position and velocity through theta
            R = np.array([[ cos(theta), sin(theta), 0],
                          [-sin(theta), cos(theta), 0],
                          [ 0,          0,          1]])
            rECEF = R*rECI
            rot = np.array([[Omega*rECI[1]], [-Omega*rECI[0]], [0]])
            vECEF = R*vECI + rot

            # Combine
            XECEF = np.concatenate( (rECEF, vECEF), axis=0 )

            return XECEF.A1

      def ECEF2LatLong( self, XECEF ) :
            # Convert ECEF properties to latitude and longitude

            # Longitude
            Long = atan2( XECEF[1], XECEF[0] ) * 180/pi

            # Latitude
            a = self.EarthProps.SemiMajorAxis
            b = self.EarthProps.SemiMinorAxis
            e = self.EarthProps.Eccentricity
            e2 = sqrt( (a**2 - b**2)/b**2 )
            p = sqrt( XECEF[0]**2 + XECEF[1]**2 )
            theta = atan2( a*XECEF[2], p*b )
            Lat = atan2( XECEF[2] + e2**2*b*(sin(theta))**3, p - e**2*a*(cos(theta))**3 ) * 180/pi

            # Altitude
            N = a/sqrt( 1 - e**2*(sin(theta))**2 )
            Alt = p/cos(theta) - N

            # Combine
            LatLongAlt = np.array( [[Lat], [Long], [Alt]] )

            return LatLongAlt[:,0]

      def ECEF2Geomagnetic( self, XECEF ) :
            # Convert ECEF properties to geomagnetic

            return 0