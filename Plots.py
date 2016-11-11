# Plotting functions for orbital and attitude dynamics model
# Murray Ireland
# October 2016

import PIL
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.basemap import Basemap

import EarthProps

def Plot3DPath( Data ) :
    # Plot 3D orbital path of satellite around Earth

    # Get axes limits
    rmax = np.max(np.abs( Data.States[0:3,:] ))

    # Plot path
    fig = plt.figure( dpi=267, figsize=[7,5], facecolor='black' )
    ax = fig.add_subplot( 111, projection='3d', axisbg='black' )
    ax.set_aspect( 'equal' )
    plt.plot( Data.States[0,:], Data.States[1,:], Data.States[2,:], 'red' )
    # plt.xlabel('x [m]')
    ax.set_xlim([-rmax, rmax])
    ax.set_ylim([-rmax, rmax])
    ax.set_zlim([-rmax, rmax])
    EarthSphere( ax )
    plt.axis( 'off' )
    plt.tight_layout()
    plt.show( block=False )

def PlotGroundTrack( Data ) :
    # Plot ground track of satellite path over Earth map

    # Test figure
    # testfig = plt.figure( dpi=267, figsize=[7,5], facecolor='white' )
    # testax = testfig.add_subplot( 111, axisbg='white' )
    # plt.plot( Data.LatLongAlt[1,:], Data.LatLongAlt[0,:] )
    # plt.axis( 'equal' )
    # plt.grid( 'on' )

    # Map figure
    fig = plt.figure( dpi=267, figsize=[7,5], facecolor='white' )
    map = Basemap( projection='mill', resolution='l' )
    # map = Basemap( projection='ortho', lat_0=45, lon_0=0, resolution='l' )
    map.drawcoastlines( linewidth=0.5 )
    map.fillcontinents( color='coral', lake_color='aqua' )
    map.drawmapboundary( fill_color='aqua' )
    map.drawmeridians( np.arange(0,360,30) )
    map.drawparallels( np.arange(-90,90,30) )

    # Split path into segments
    segment = np.vstack( (Data.LatLongAlt[1,:], Data.LatLongAlt[0,:]) )
    thresh = 90.
    isplit = np.nonzero( np.abs( np.diff( segment[0] ) ) > thresh )[0]
    subsegs = np.split( segment, isplit+1, axis=+1 )

    for seg in subsegs :
        x, y = map( seg[0], seg[1] )
        map.plot( x, y, c='red' )
    
    plt.show( block=False )


def PlotAttitude( Data ) :
    # Plot attitude response of satellite

    Euler = Quaternion2Euler( Data.States[6:10,:] )

    fig = plt.figure( dpi=150, figsize=[12,10], facecolor='white' )
    ax1 = fig.add_subplot( 321, axisbg='white' )
    plt.plot( Data.Time, Euler[0,:] )
    plt.xlabel( 'Time [s]' )
    plt.ylabel( 'Roll [deg]' )

    ax2 = fig.add_subplot( 322, axisbg='white' )
    plt.plot( Data.Time, Data.States[10,:] )
    plt.xlabel( 'Time [s]' )
    plt.ylabel( 'Roll rate [rad/s]')

    ax3 = fig.add_subplot( 323, axisbg='white' )
    plt.plot( Data.Time, Euler[1,:] )
    plt.xlabel( 'Time [s]' )
    plt.ylabel( 'Pitch [deg]')

    ax4 = fig.add_subplot( 324, axisbg='white' )
    plt.plot( Data.Time, Data.States[11,:] )
    plt.xlabel( 'Time [s]' )
    plt.ylabel( 'Pitch rate [rad/s]')

    ax3 = fig.add_subplot( 325, axisbg='white' )
    plt.plot( Data.Time, Euler[2,:] )
    plt.xlabel( 'Time [s]' )
    plt.ylabel( 'Yaw [deg]')

    ax4 = fig.add_subplot( 326, axisbg='white' )
    plt.plot( Data.Time, Data.States[12,:] )
    plt.xlabel( 'Time [s]' )
    plt.ylabel( 'Yaw rate [rad/s]')

    return 0

def EarthSphere( ax ) :
    # Create Earth spheroid for 3D plots

    lons = np.linspace( -180, 180, 100 ) * np.pi/180
    lats = np.linspace( -90, 90, 100 ) * np.pi/180

    x = EarthProps.SemiMajorAxis*np.outer( np.cos(lons), np.cos(lats) ).T
    y = EarthProps.SemiMajorAxis*np.outer( np.sin(lons), np.cos(lats) ).T
    z = EarthProps.SemiMinorAxis*np.outer( np.ones( np.size(lons) ), np.sin(lats) ).T
    ax.plot_surface( x, y, z, rstride=4, cstride=4, color='blue', edgecolors='cyan' )

def Quaternion2Euler( Q ) :
    # Convert attitude in quaternions to Euler angles

    Euler = np.zeros( (3, Q.shape[1]) )

    Euler[0,:] = np.arctan( 2*(Q[2,:]*Q[3,:] + Q[0,:]*Q[1,:]) / ( Q[0,:]**2 - Q[1,:]**2 - Q[2,:]**2 + Q[3,:]**2 ) )
    Euler[1,:] = np.arcsin( 2*(Q[1,:]*Q[3,:] - Q[0,:]*Q[2,:]) ) * 180/np.pi
    Euler[2,:] = np.arctan( 2*(Q[1,:]*Q[2,:] + Q[0,:]*Q[3,:]) / ( Q[0,:]**2 + Q[1,:]**2 - Q[2,:]**2 - Q[3,:]**2 ) )

    return Euler