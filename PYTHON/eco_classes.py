#!/usr/bin/python
# Filename: eco_classes.py

class LicenseError(Exception):
    pass

class SpatialRefProjError (Exception):
    pass


import arcpy as ap
import numpy as np
import os, sys, time, glob, math, string, shutil
from math import pi, sin, cos, tan, sqrt


#    
# ASSORTMENT OF TERRAIN DERIVATIVES BASED ON DEM
#

def calc_roughness(dem):
    '''Calculate and reclassifiy the roughness of the terrain

        @type   dem: c{str}
        @param  dem: the path of the elevation dem
        @rtype  returns an in memory raster
    '''

    # START SPATIAL ANALYST AND GET DEM METADATA
    dem_desc = ap.Describe(dem)
    ext = dem_desc.extent
    rect = str(str(ext.XMin) + " " + str(ext.YMin) + " " + str(ext.XMax) + " " + str(ext.YMax))
    
    # CHECK RASTER ROJECTION
    if dem_desc.spatialReference.type != "Projected":
        raise SpatialRefProjError
    
    # CALC DEM STD AND TAKE THE SQ ROOT
    focalstats_tmp = ap.sa.FocalStatistics(dem, "#", "STD")
    roughness = ap.sa.SquareRoot(focalstats_tmp)

    return roughness


def calc_hli(dem):
    '''Calculates the Heat Load Index

        @type   dem: c{str}
        @param  dem: the path of the elevation dem
        @rtype  returns an in memory raster to be saved
    '''

    desc_raster = ap.Describe(dem)
    extent = desc_raster.extent
    
    lat = utm_to_lat(dem)
    l = lat * 57.296

    cl = cos(lat)
    sl = cos(l)

    slope = ap.sa.Slope(dem) * 0.017453293
    aspect = ap.sa.Aspect(dem)
    abs_aspect = ap.sa.Abs(180-ap.sa.Abs(aspect)) * 0.017453293
    cos_slope = ap.sa.Cos(slope)
    sin_slope = ap.sa.Sin(slope)
    cos_abs_aspect = ap.sa.Cos(abs_aspect)
    sin_abs_aspect = ap.sa.Sin(abs_aspect)

    calc_raster = ap.sa.Exp(-1.467 + 1.582 * cl * cos_slope - 1.5 * cos_abs_aspect * sin_slope * sl - 0.262 * sl * sin_slope + 0.607 * sin_abs_aspect * sin_slope)

    m = ap.GetRasterProperties_management(calc_raster, "MAXIMUM")
    m_ras = ap.sa.CreateConstantRaster(m, "FLOAT", desc_raster.MeanCellHeight, extent)
    
    out_raster = ap.sa.Int(((calc_raster/m_ras)*10000) + 0.5)

    return out_raster


#
# HELPER FUNCTIONS
#

def checkout_ext(ext_type):
    if ap.CheckExtension(ext_type) == 'Available':
        ap.CheckOutExtension(ext_type)
        print "Checking out " + ext_type + " Extension"
    else:
        raise LicenseError
        print "Cannot checkout " + ext_type + " Extension"

def checkin_ext(ext_type):
    ap.CheckInExtension(ext_type)
    print "Checking in " + ext_type + " Extension"

# CALCULATE AND RETURN LATITUDE OF RASTER CENTER FROM WGS84, UTM NAD83
def utm_to_lat(raster):
    '''Convert a raster's UTM Northing coordinate center to Latitude'''
    # GET UTM CENTER AND SETUP K CONSTANT
    extent = ap.Describe(raster).extent
    
    y = (extent.YMax - extent.YMin)/2 + extent.YMin
    x = ((extent.XMax - extent.XMin)/2 + extent.XMin) - 500000.0
    k0 = float(0.9996)
    _rad2deg = 180.0 / pi
    
    e2 = float(0.00669438)
    e1 = (1-sqrt(1-e2))/(1+sqrt(1-e2))
    
    a = 6378137
    b = float(6356752.3142)
    n = (a-b)/(a+b)

    eprime2 = (e2)/(1-e2)
    
    M = y/k0
    mu = M/(a*(1-e2/4-3*e2*e2/64-5*e2*e2*e2/256))

    phi1_rad = (mu + (3*e1/2-27*e1*e1*e1/32)*sin(2*mu) 
               + (21*e1*e1/16-55*e1*e1*e1*e1/32)*sin(4*mu)
               +(151*e1*e1*e1/96)*sin(6*mu))

    phi1 = phi1_rad*_rad2deg

    N1 = a/sqrt(1-e2*sin(phi1_rad)*sin(phi1_rad))
    T1 = tan(phi1_rad)*tan(phi1_rad)
    C1 = eprime2*cos(phi1_rad)*cos(phi1_rad)
    R1 = a*(1-e2)/pow(1-e2*sin(phi1_rad)*sin(phi1_rad), 1.5)
    D = x/(N1*k0)

    Lat = phi1_rad - (N1*tan(phi1_rad)/R1)*(D*D/2-(5+3*T1+10*C1-4*C1*C1-9*eprime2)*D*D*D*D/24
                    +(61+90*T1+298*C1+45*T1*T1-252*eprime2-3*C1*C1)*D*D*D*D*D*D/720)

    lat = Lat * _rad2deg
    
    return float(lat)

# NORMALIZE AND STRETCH DATA BETWEEN 0 AND 255
def stretch(raster, out_path):
    ''' Normalize and stretch raster data from 0 to 255
    '''

    print "Stretching " + ap.Describe(raster).baseName
    # GET MIN AND MAX VALUES FROM SELECTED RASTER
    max_val = ap.GetRasterProperties_management(raster, "MAXIMUM")
    min_val = ap.GetRasterProperties_management(raster, "MINIMUM")

    # GET RASTER EXTENT and Cell Size
    desc = ap.Describe(raster)
    ap.env.cellSize = raster
    ap.env.extent = raster

    # CREATE RASTER CONSTANTS
    max_raster = ap.sa.CreateConstantRaster(max_val, "FLOAT", desc.MeanCellHeight, ap.env.extent)
    min_raster = ap.sa.CreateConstantRaster(min_val, "FLOAT", desc.MeanCellHeight, ap.env.extent)

    # SET THE NEW RASTER RANGE
    INPUT_MAX = 255
    INPUT_MIN = 0

    input_max_raster = ap.sa.CreateConstantRaster(INPUT_MAX, "FLOAT", desc.MeanCellHeight, ap.env.extent)
    input_min_raster = ap.sa.CreateConstantRaster(INPUT_MIN, "FLOAT", desc.MeanCellHeight, ap.env.extent)

    out_raster = ((raster - min_raster) * input_max_raster) / ((max_raster - min_raster) + input_min_raster)
    output_name = desc.baseName + ".img"
##    out_path = os.path.join(ap.env.workspace, name)

    if os.path.exists(out_path) != True:
        os.mkdir(out_path)
    
    out_raster.save(os.path.join(out_path, output_name))

# CLIP RASTER TO AOI BOUNDARY DERIVED FROM THE boundary() function
def clip(raster, out_path, name, boundary):
    '''Clip raster to designated boundary
    '''

    if os.path.exists(out_path) != True:
        os.mkdir(out_path)
        
    output = os.path.join(out_path, name)
    ap.Clip_management(raster, boundary, output, "#", "-999", "NONE")

# RETURNS XY MIN/MAX FROM A RASTER OR VECTOR
def boundary(raster):
    ''' Get the extent boundary of a raster
    '''

    ext = ap.Describe(raster).extent
    rect = str(str(ext.XMin) + " " + str(ext.YMin) + " " + str(ext.XMax) + " " + str(ext.YMax))

    return rect


##
##  CREATING THE OUTPUT DATA
##

def run_dem(dem, out_path):
    ''' Calculate all of the Terrain DEM Derivatives and normalize the clipped directory
    '''

    # Checkout Extension and Set workspace
    ap.env.workspace = out_path
    ap.env.overwriteOutput = True
    
    if os.path.exists(out_path) != True:
       os.mkdir(out_path)

    # Calc Slope
    print "Calculating Slope"
    slope = ap.sa.Slope(dem, "PERCENT_RISE")
    slope.save(os.path.join(out_path, "slope.img"))

    # Calc roughness
    print "Calculating Roughness"
    roughness = calc_roughness(dem)
    roughness.save(os.path.join(out_path, "roughness.img"))

    # Calc solar radiation
    print "Calculating Heat Load Index"
    hli = calc_hli(dem)
    hli.save(os.path.join(out_path, "hli.img"))


def aoi_data_prep(ndvi, satv, dem, aoi, out_path):
    '''Returns the unsupervised ecosite classification in a given AOI
    '''

    if os.path.exists(out_path) != True:
        os.mkdir(out_path)

    ap.env.workspace = out_path
    bndry = boundary(aoi)

    # CREATE TEMPORARY "TEMP" DIRECTORY
    temp = os.path.join(out_path, "temp")

    if os.path.exists(temp) != True:
        os.mkdir(temp)

    # CREATE FINAL OUTPUT "COMPLETE" DIRECTORY
    complete = os.path.join(out_path, "complete")

    if os.path.exists(complete) != True:
        os.mkdir(complete)

    checkout_ext("Spatial")
    
    # cLIP DEM TO AOI
    print "Clipping DEM to AOI"
    clip(dem, temp, "dem.img", bndry)

    # Clip NDVI to AOI
    print "Clipping NDVI to AOI"
    clip(ndvi, temp, "ndvi.img", bndry)

    # CLIP SATV TO AOI
    print "Clipping SATV to AOI"
    clip(satv, temp, "satv.img", bndry)

    # SET AOI DEM TO CALCULATE SURFACES
    aoi_dem = os.path.join(temp, "dem.img")

    # RUN THE DEM SURFACE DERIVATIVES
    run_dem(aoi_dem, temp)

    ap.env.workspace = temp
    rasters = ap.ListRasters()

    for img in rasters:
        name = ap.Describe(img).file
        raster_path = os.path.join(ap.env.workspace, name)
        stretch(raster_path, complete)

    checkin_ext("Spatial")

def rm_tmp_data(out_path):
    ''' Remove directories and data 
    '''
    temp = os.path.join(out_path, "temp")

    try:
        shutil.rmtree(temp)
    except:
        print "Unsuccessful deletion of temporary data."
    finally:
        print "Temporary data cleanup complete."

def run_iso(in_path, out_path=None):
    ''' Run an Unspervised classification on the list of rasters in the specified in_path directory
    '''
    # SET UP ENVIRONMENT AND VARIABLES
    ap.env.workspace = in_path
    rasters = ap.ListRasters()
    checkout_ext("Spatial")

    if out_path is None:
        out_path = os.path.abspath(os.path.join(in_path, os.pardir))

    out_name = os.path.join(out_path, "clusters.img")

    clusters = ap.sa.IsoClusterUnsupervisedClassification(rasters, 12, "#", "#")
    clusters.save(out_name)

    checkin_ext("Spatial")
    
# TEST RUN THE DATA

def output_data(ndvi, satv, dem, aoi, out_path):
    aoi_data_prep(ndvi, satv, dem, aoi, out_path)

    print "Removing Temp Data"
    rm_tmp_data(out_path)

    complete = os.path.join(out_path, "complete")
    print "Running IsoCluster"
    run_iso(complete)

dem = "F:/Elevation_Data/NED_10m/Lattices/Pinal/ned10m_1-1"
aoi = "c:/Users/andrew.burnes/Documents/GIS/response-units/test/data/aoi_2.shp"
ndvi = "c:/landsat/processed/toa/LC80370372014028LGN00_NDVI.img"
satv = "F:/Raster_Data/Cover/cover_spring11.img"
output = "c:/Users/andrew.burnes/Documents/GIS/response-units/test/data/output"

output_data(ndvi, satv, dem, aoi, output)


