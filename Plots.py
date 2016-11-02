# Plotting functions for orbital and attitude dynamics model
# Murray Ireland
# October 2016

import PIL
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.basemap import Basemap
import EarthProps

def EarthSphere( ax ) :
    # Create Earth spheroid for 3D plots

    lons = np.linspace( -180, 180, 100 ) * np.pi/180
    lats = np.linspace( -90, 90, 100 ) * np.pi/180

    x = EarthProps.SemiMajorAxis*np.outer( np.cos(lons), np.cos(lats) ).T
    y = EarthProps.SemiMajorAxis*np.outer( np.sin(lons), np.cos(lats) ).T
    z = EarthProps.SemiMinorAxis*np.outer( np.ones( np.size(lons) ), np.sin(lats) ).T
    ax.plot_surface( x, y, z, rstride=4, cstride=4, color='blue' )

def Plot3DPath( Data ) :
    # Plot 3D orbital path of satellite around Earth

    fig = plt.figure( dpi=267, figsize=[7,5], facecolor='black' )
    ax = fig.add_subplot( 111, projection='3d', axisbg='black' )
    plt.plot( Data.States[0,:], Data.States[1,:], Data.States[2,:], 'red' )
    EarthSphere( ax )
    plt.axis( 'equal' )
    plt.axis( 'off' )
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
    map.drawcoastlines( linewidth=0.5 )
    map.fillcontinents( color='coral', lake_color='aqua' )
    map.drawmapboundary( fill_color='aqua' )
    map.drawmeridians( np.arange(0,360,30) )
    map.drawparallels( np.arange(-90,90,30) )
    x, y = map( Data.LatLongAlt[1,:], Data.LatLongAlt[0,:] )
    map.plot( x, y )
    plt.show( block=False )