#!/usr/bin/env python
# -*- coding: utf-8 -*-
# python

import numpy as np
import shapefile

from . import cluster
from . import geometry
from . import geo2rot


def shp(sf_path, rlat_d, rlon_d, udir_d, uheight1_d, rlat_v, rlon_v,
        FR_URBAN, FR_ROOF):
    # Inizializations
    AREA_BLD = np.zeros((1, rlat_d, rlon_d))
    BUILD_W = np.zeros((1, udir_d, rlat_d, rlon_d))
    FR_STREETD = np.zeros((1, udir_d, rlat_d, rlon_d))
    LAMBDA_F = np.zeros((1, udir_d, rlat_d, rlon_d))
    LAMBDA_P = np.zeros((1, rlat_d, rlon_d))
    MEAN_HEIGHT = np.zeros((1, rlat_d, rlon_d))
    STREET_W = np.zeros((1, udir_d, rlat_d, rlon_d))
    VERT_AREA = np.zeros((1, udir_d, rlat_d, rlon_d))
    eps = 1e-4
    # Read the dataset
    sf = shapefile.Reader(sf_path)
    shapes = sf.shapes()
    # Calculations
    N = len(shapes)
    for x in range(0, N):
        # Reading Geometry
        shape = shapes[x]
        p = shape.points
        p = np.array(p)
        # Reading the Area (horizontal)
        area = sf.record(x)[4]
        # Calculating the coordinates of the centroid of the polygon
        lonC = np.mean(p[:, 0])
        latC = np.mean(p[:, 1])
        # Converting centroid to rotated coordinates
        [lonC_r, latC_r] = geo2rot.g2r(lonC, latC)
        # Calculating the index of the correspoding grid point
        lon_idx = np.abs(rlon_v - lonC_r).argmin()
        lat_idx = np.abs(rlat_v - latC_r).argmin()
        # Allocating the total building area
        AREA_BLD[0, lat_idx, lon_idx] += area
        # Clustering the geomerty heights
        hgt = sf.record(x)[0]
        MEAN_HEIGHT[0, lat_idx, lon_idx] += hgt * area
        hgt_class = cluster.height_old(hgt)
        # Looping over building segments (facades)
        for k in range(1, len(p)):
            vert_area = geometry.distWGS(p, k) * hgt
            ang = geometry.angle(p, k)
            ang_class = cluster.angle(ang)
            #
            VERT_AREA[0, ang_class, lat_idx, lon_idx] += vert_area
            # Allocating the building area by direction and height
            FR_ROOF[0, ang_class, hgt_class, lat_idx, lon_idx] += vert_area

    # Calculating the area densities (Grimmond and Oke, 1998)
    np.seterr(divide='ignore')  # disabled division by 0 warning
    area_grid = 77970  # m2, calculated in GIS. TO DO: calculate from grid
    lambda_p = AREA_BLD[0, :, :]/(area_grid*FR_URBAN[0, :, :])
    lambda_p[lambda_p > 0.9] = 0.9  # upper limit
    lambda_p[lambda_p < 0.1] = 0.1  # lower limit
    lambda_f = VERT_AREA[0, :, :, :]/(area_grid*FR_URBAN[0, :, :])
    lambda_f[lambda_f > 0.9] = 0.9  # upper limit
    lambda_f[lambda_f < 0.1] = 0.1  # lower limit
    LAMBDA_P[0, :, :] = lambda_p
    LAMBDA_F[0, :, :, :] = lambda_f
    # Calculating the street and building widths (Martilli, 2009)
    h_m = MEAN_HEIGHT[0, :, :] / AREA_BLD[0, :, :]
    h_m[h_m == 15] = 10.
    BUILD_W[0, :, :, :] = lambda_p[np.newaxis, :, :] / lambda_f[:, :, :] * h_m
    STREET_W[0, :, :, :] = (1 / lambda_p[np.newaxis, :, :] - 1) * lambda_p[np.newaxis, :, :] \
        / lambda_f * h_m
    # STREET_W[STREET_W<5] = 5  # min aspect ration LCZ 2
    # STREET_W[STREET_W>100] = 100  # max aspect ration LCZ 9

    # Calculating and normalizing the canyon direction distribution
    FR_STREETD = np.sum(FR_ROOF, 2)
    norm_streetd = np.sum(FR_STREETD, 1)
    FR_STREETD = FR_STREETD/norm_streetd

    # Normalizing FR_ROOF
    norm_fr_roof = np.sum(FR_ROOF, 1)
    for j in range(0, udir_d):
        for k in range(0, uheight1_d):
            for o in range(0, rlat_d):
                for z in range(0, rlon_d):
                    if norm_fr_roof[0, k, o, z] != 0:
                        FR_ROOF[0, j, k, o, z] = FR_ROOF[0, j,
                                                         k, o, z]/norm_fr_roof[0, k, o, z]

    FR_ROOF[FR_ROOF < 1e-8] = 0  # to avoid negative values

    return BUILD_W, STREET_W, FR_ROOF, FR_STREETD, shapes


def shp_alternate(sf_path, rlat_d, rlon_d, udir_d, uheight1_d, rlat_v, rlon_v,
                  lat_mid, dlat, dlon, FR_URBAN, FR_ROOF):
    # Inizializations
    AREA_BLD = np.zeros((1, rlat_d, rlon_d))
    BUILD_W = np.zeros((1, udir_d, rlat_d, rlon_d))
    FR_STREETD = np.zeros((1, udir_d, rlat_d, rlon_d))
    LAMBDA_F = np.zeros((1, udir_d, rlat_d, rlon_d))
    LAMBDA_P = np.zeros((1, rlat_d, rlon_d))
    MEAN_HEIGHT = np.zeros((1, rlat_d, rlon_d))
    STREET_W = np.zeros((1, udir_d, rlat_d, rlon_d))
    VERT_AREA = np.zeros((1, udir_d, rlat_d, rlon_d))
    # Read the dataset
    sf = shapefile.Reader(sf_path)
    shapes = sf.shapes()
    # Calculations
    area_grid = geometry.surface_element_earth(dlat, dlon, lat_mid)
    print("dlat: {}, dlon: {}, lat_mid: {}".format(dlat, dlon, lat_mid))
    print(area_grid)
    N = len(shapes)
    for x in range(0, N):
        print(x)
        # Reading Geometry
        shape = shapes[x]
        hgt = sf.record(x)[-5]
        if hgt < 1e-1:
            continue
        p = shape.points
        p = np.array(p)
        # Reading the Area (horizontal)
        area = sf.record(x)[-1]
        # Calculating the coordinates of the centroid of the polygon
        lonC = np.mean(p[:, 0])
        latC = np.mean(p[:, 1])
        # Converting centroid to rotated coordinates
        #[lonC_r,latC_r] = geo2rot.g2r(lonC,latC)
        # Calculating the index of the correspoding grid point
        lon_idx = np.abs(rlon_v - lonC).argmin()
        #print(lon_idx)
        lat_idx = np.abs(rlat_v - latC).argmin()
        #print(lat_idx)
        # Allocating the total building area
      
        AREA_BLD[0, lat_idx, lon_idx] += area
        # Clustering the geomerty heights
       
        MEAN_HEIGHT[0, lat_idx, lon_idx] += hgt * area
        hgt_class = cluster.height_old(hgt)
        # Looping over building segments (facades)
        for k in range(1, len(p)):
            vert_area = geometry.distWGS(p, k) * hgt
            ang = geometry.angle(p, k)
            ang_class = cluster.angle(ang)
            #
            VERT_AREA[0, ang_class, lat_idx, lon_idx] += vert_area
            # Allocating the building area by direction and height
            FR_ROOF[0, ang_class, hgt_class, lat_idx, lon_idx] += vert_area
    
    np.seterr(divide='ignore')  # disabled division by 0 warning
    FR_BUILD = (AREA_BLD[0, :, :] / area_grid)
    print("area_grid: {}".format(area_grid))
    print("FR_BULD MAX {}".format(FR_BUILD.max()))
    FR_BUILD[FR_BUILD>0.9] = 0.9
    FR_URBAN[0, :, :] = FR_BUILD + (1 - FR_BUILD)/2
    FR_URBAN[0, :, :][FR_BUILD < 0.1] = 0
    FR_URBAN[0, :, :][FR_URBAN[0, :, :] > 1] = 1


    # Calculating the area densities (Grimmond and Oke, 1998)
    lambda_p = AREA_BLD[0, :, :]/(area_grid*FR_URBAN[0, :, :])
    lambda_p[lambda_p > 0.9] = 0.9  # upper limit
    lambda_p[lambda_p < 0.1] = 0.1  # lower limit
    lambda_f = VERT_AREA[0, :, :, :]/(area_grid*FR_URBAN[0, :, :])
    lambda_f[lambda_f > 0.9] = 0.9  # upper limit
    lambda_f[lambda_f < 0.1] = 0.1  # lower limit
    LAMBDA_P[0, :, :] = lambda_p
    LAMBDA_F[0, :, :, :] = lambda_f
    # Calculating the street and building widths (Martilli, 2009)
    h_m = MEAN_HEIGHT[0, :, :] / AREA_BLD[0, :, :]
    h_m[h_m == 15] = 10.
    BUILD_W[0, :, :, :] = lambda_p[np.newaxis, :, :] / lambda_f[:, :, :] * h_m
    STREET_W[0, :, :, :] = (1 / lambda_p[np.newaxis, :, :] - 1) * lambda_p[np.newaxis, :, :] \
        / lambda_f * h_m
    # STREET_W[STREET_W<5] = 5  # min aspect ration LCZ 2
    # STREET_W[STREET_W>100] = 100  # max aspect ration LCZ 9

    # Calculating and normalizing the canyon direction distribution
    #FR_ROOF[0,:,-1:-2,:,:] = 0

    FR_STREETD = np.sum(FR_ROOF, 2)
    norm_streetd = np.sum(FR_STREETD, 1)
    FR_STREETD = FR_STREETD/norm_streetd

    # Normalizing FR_ROOF
    norm_fr_roof = np.sum(FR_ROOF, 1)
    for j in range(0, udir_d):
        for k in range(0, uheight1_d):
            for o in range(0, rlat_d):
                for z in range(0, rlon_d):
                    if norm_fr_roof[0, k, o, z] != 0:
                        FR_ROOF[0, j, k, o, z] = FR_ROOF[0, j, k, o, z]/norm_fr_roof[0, k, o, z]
#                        if FR_ROOF[0, j, :, o, z].sum() < 1e-3:
#                                print(FR_ROOF[0, j, :, o, z].sum())
#                                FR_URBAN[0, o, z] = 0
#                                FR_STREETD[0, j, o, z] = 0
    # for j in range(udir_d):
    #    for k in range(uheight1_d):
    #         FR_URBAN[0,:,:][FR_ROOF[0,j,k,:,:] < eps] = 0
#    for i in range (rlon_d):
#        for j in range (rlat_d):
#            if FR_ROOF[0,:,:,i,j].sum() < 1e-5:
#                FR_URBAN[0,i,j] = 0
    #FR_ROOF[0,:,0,:,:] = 0
    FR_ROOF[FR_ROOF < 0] = 0  # to avoid negative values
    print(FR_URBAN[0, :, :].max())
    

    return BUILD_W, STREET_W, FR_ROOF, FR_STREETD, FR_URBAN, shapes
