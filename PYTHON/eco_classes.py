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

# Transform all terain derrevatives
def terrain_der(path):
    '''Transform the DEM elevation data into slope, radiation, and ruggedness

        @type    path: c{str}
        @param   path: Path to DEM file directory
        @rtype   raster
    '''

    # Calc Slope
    slope = ap.sa.Slope(path, "PERCENT_RISE")
    slope.save("slope.img")

    # Calc roughness
    roughness = calc_roughness(path)
    roughness.save("roughness.img")

    # Calc solar radiation
    hli = calc_hli(path)
    hli.save("hli.img")
    
# Individual Terrain Derrivavtive functions
def calc_roughness(dem):
    '''Calculate and reclassifiy the roughness of the terrain

        @type   dem: c{str}
        @param  dem: the path of the elevation dem
        @rtype  returns an in memory raster
    '''

    # START SPATIAL ANALYST AND GET DEM METADATA
    # checkout_ext("Spatial")
    dem_desc = ap.Describe(dem)
    ext = dem_desc.extent
    rect = str(str(ext.XMin) + " " + str(ext.YMin) + " " + str(ext.XMax) + " " + str(ext.YMax))
    
    # CHECK RASTER ROJECTION
    if dem_desc.spatialReference.type != "Projected":
        raise SpatialRefProjError
    
    # CALC DEM STD AND TAKE THE SQ ROOT
    focalstats_tmp = ap.sa.FocalStatistics(dem, "#", "STD")
    roughness = ap.sa.SquareRoot(focalstats_tmp)
    rough_slice = ap.sa.Slice(roughness, 3, "NATURAL_BREAKS")
    aggr_10 = ap.sa.Aggregate(rough_slice, 10, "MEDIAN", "TRUNCATE")
    filter_rough = ap.sa.MajorityFilter(aggr_10, "EIGHT", "MAJORITY")

    return filter_rough


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

        
# Helper functions
def checkout_ext(ext_type):
    if ap.CheckExtension(ext_type) == 'Available':
        ap.CheckOutExtension(ext_type)
        print "\nChecking out " + ext_type + " Extension"
    else:
        raise LicenseError
        print "\nCannot checkout " + ext_type + " Extension"

def checkin_ext(ext_type):
    ap.CheckInExtension(ext_type)
    print "\nChecking in " + ext_type + " Extension"

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

def stretch(raster, out_path):
    ''' Normalize and stretch raster data from 0 to 255
    '''
    # GET MIN AND MAX VALUES FROM SELECTED RASTER
    max_val = ap.GetRasterProperties_management(raster, "MAXIMUM")
    min_val = ap.GetRasterProperties_management(raster, "MINIMUM")

    # GET RASTER EXTENT and Cell Size
    extent = ap.Describe(raster).extent
    ap.env.cellSize = raster
    ap.env.extent = raster

    # CREATE RASTER CONSTANTS
    max_raster = ap.sa.CreateConstantRaster(max_val, "INTEGER")
    min_raster = ap.sa.CreateConstantRaster(min_val, "INTEGER")

    # SET THE NEW RASTER RANGE
    INPUT_MAX = 255
    INPUT_MIN = 0

    input_max_raster = ap.sa.CreateConstantRaster(INPUT_MAX, "INTEGER")
    input_min_raster = ap.sa.CreateConstantRaster(INPUT_MIN, "INTEGER")

    out_raster = (raster - min_raster) * input_max_raster / (max_raster - min_raster) + input_min_raster    
    output_name = ap.Describe(raster).baseName + ".img"
##    out_path = os.path.join(ap.env.workspace, name)

    if os.path.exists(out_path) != True:
        os.mkdir(out_path)
    
    out_raster.save(os.path.join(out_path, output_name))

def clip(raster, out_path, name, boundary):
    '''Clip raster to designated boundary
    '''

    if os.path.exists(out_path) != True:
        os.mkdir(out_path)

    print boundary
        
    output = os.path.join(out_path, name)

    ap.Clip_management(raster, boundary, output, "#", "-999", "NONE")

def boundary(raster):
    ''' Get the extent boundary of a raster
    '''

    ext = ap.Describe(raster).extent
    rect = str(str(ext.XMin) + " " + str(ext.YMin) + " " + str(ext.XMax) + " " + str(ext.YMax))

    return rect

    

# TEST RUN THE DATA

dem = "c:/Users/andrew.burnes/Documents/GIS/response-units/test/data/dem.img"
aoi = "c:/Users/andrew.burnes/Documents/GIS/response-units/test/data/aoi.shp"
output = "c:/Users/andrew.burnes/Documents/GIS/response-units/test/data/output"

def run_test(dem, aoi, out_path):

    # Checkout Extension and Set workspace
    checkout_ext("Spatial")
    ap.env.workspace = out_path
    ap.env.overwriteOutput = True

    out_bndry = boundary(aoi)
    
    if os.path.exists(out_path) != True:
       os.mkdir(out_path)

    final_output = os.path.join(ap.env.workspace, "final")

    if os.path.exists(final_output) != True:
        os.mkdir(final_output)

    # Calc Slope
    print "Calculating Slope"
    slope = ap.sa.Slope(dem, "PERCENT_RISE")
    slope.save(os.path.join(final_output, "slope.img"))

    # Calc roughness
    print "Calculating roughness"
    roughness = calc_roughness(dem)
    roughness.save(os.path.join(final_output, "roughness.img"))
    #clip(roughness, final_output, "rough.img", out_bndry)
    #roughness.save("roughness.img")

    # Calc solar radiation
    print "Calculating heat load index"
    hli = calc_hli(dem)
    hli.save(os.path.join(final_output, "hli.img"))

    # GETTING GARBAGE RASTER LIST
##    print "cOLLECTING temp Raster List"
##    tmp_list = ap.ListRasters()
    

##    for tmp in tmp_list:
##        bs_name = ap.Describe(tmp).baseName
##        tmp_path = os.path.join(ap.env.workspace, str(bs_name) + "*")
##        os.remove(tmp_path)

    ap.env.workspace = os.path.join(out_path, "final")
    rasters = ap.ListRasters()

    complete = os.path.join(os.path.abspath(os.path.join(output, os.pardir)), 'complete')
    clipped= os.path.join(os.path.abspath(os.path.join(output, os.pardir)), 'clipped')

    if os.path.exists(complete) != True:
        os.mkdir(complete)

    if os.path.exists(clipped) != True:
        os.mkdir(clipped)

    # CLIPPING OUPUT RASTERS TO AOI
    print "Clipping surface rasters to boundary"
    for img in rasters:
        print img
        name = ap.Describe(img).file
        img_path = os.path.join(ap.env.workspace, name)
        clip(img_path, clipped, name, out_bndry)

    ap.env.workspace = clipped
    clips = ap.ListRasters()

    # STRETCH THE CLIPPED DEM DEERIVATIVES
    print "Normalizing and Stretching Data to Values: 0 - 255"
    for raster in clips:
        name = ap.Describe(raster).file
        raster_path = os.path.join(ap.env.workspace, name)
        stretch(raster_path, complete)
    
    checkin_ext("Spatial")

run_test(dem, aoi, output)

print os.path.abspath(os.path.join(output, os.pardir))
