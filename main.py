#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 16:40:56 2020

@author: blobfishman
# """


import sys
import geopandas as gpd
import pandas as pd
pd.set_option('display.max_rows', 20)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import numpy as np
import time

def main(args):
    if len(args) != 5 :
        print("missing arguments. Usage:")
        print("python3 main.pyINFRA_SCORES.shp POINTS.shp POINTS_BACKGROUND.shp SCORES_FINAL.shp")
        return

    start_time = time.time()

    # Grid File no longer needed as it is also included in INFRA_SCORES File
    # print("Loading Grid Shape File at " + args[1])
    # try:
    #     grid = gpd.read_file(args[1])
    #     try:
    #         del grid['left']
    #         del grid['top']
    #         del grid['bottom']
    #         del grid['right']
    #     except:
    #         pass
    # except:
    #     print("Grid Shape File could not be located or is in the wrong format. Do you also have the .shx and .dbf files in the same directory?. Aborting..")
    #     return

    print("Loading Infra Scores Shape File at " + args[1])

    try:
        infra = gpd.read_file(args[1])
        grid  = gpd.read_file(args[1])
        del grid["infra_scor"]
    except:
        print("Points Infra Scores Shape File could not be located or is in the wrong format. Do you also have the .shx and .dbf files in the same directory?. Aborting..")
        return

    print("Loading Points Shape File at " + args[2])

    try:
        points = gpd.read_file(args[2])
        try:
            del points['field_1']
            del points['field_2']
        except:
            pass
    except:
        print("Points Shape File could not be located or is in the wrong format. Do you also have the .shx and .dbf files in the same directory?. Aborting..")
        return

    print("Loading Points Background Shape File at " + args[3])

    try:
        noise_points = gpd.read_file(args[3])
        try:
            del points['field_1']
            del points['field_2']
        except:
            pass
    except:
        print("Points Background Shape File could not be located or is in the wrong format. Do you also have the .shx and .dbf files in the same directory?. Aborting..")
        return

    #points = gpd.read_file(args[2])

    #print("points")
    #print(points)
    print("Calculating final scores...")
    # join points
    dfsjoin = gpd.sjoin(grid, points, how="left", op='contains')
    #print("join")
    #print(dfsjoin)

    # Super ekliges Workaround:
    # Zählfunktion ignoritert Geometry und nan
    # TODO pivot_table fixen
    def nanlen(lst):
        ct = 0
        for l in lst:
            if str(type(l)) == "<class 'float'>":
                if  not str(l) == "nan":
                    #print(l)
                    ct+=1
        return ct
    # dropna macht nichts?
    dfpivot = pd.pivot_table(dfsjoin,index='id', aggfunc={'index_right':nanlen}, dropna=True)
    #print("pivot")
    #print(dfpivot)

    dfmerge = grid.merge(dfpivot, how='left',on='id')
    #print("merge")
    #print(dfmerge)

    # load infrastructure
    #infra = gpd.read_file(args[3])
    #print("infra")
    #print(infra)
    # join noise points
    #noise_points = gpd.read_file(args[4])

    dfsjoin_noise = gpd.sjoin(grid, noise_points, how="left", op='contains')
    dfpivot_noise = pd.pivot_table(dfsjoin_noise,index='id', aggfunc={'index_right':nanlen})

    dfmerge_noise = grid.merge(dfpivot_noise, how='left',on='id')
    #print("merge noise")
    #print(dfmerge_noise)

    # # calculate value
    val_points = dfmerge["index_right"]
    val_points_noise = dfmerge_noise["index_right"]
    y  = infra["infra_scor"]

    #Grundlage
    '''x = np.maximum(val_points-val_points_noise, np.ones(len(y)))
    # normalisieren
    val = x**y + x # Funktion'''

    #Normalisieren

    #Problem: Exponentielle Ansteckungsgefahr bei lin. Wachstum der Menschenanzahl wird eliminiert.
    '''x = np.maximum(val_points-val_points_noise, np.ones(len(y)))
    x = np.log(x) + 1
    val = x**y + x # Funktion'''

    #Idee/Ansatz: Kleinere Basis als exp(x) um exponentielles Wachstum zu erhalten.
    x = np.maximum(val_points-val_points_noise, np.ones(len(y)))
    ''' Wurzel entfernt: x = 2**np.log(x)
    Warum sollten wir erhöhte Menschenaufkommen
    bestrafen? Genau darum geht es doch, oder?'''
    val = x**y

    '''Problem/Ansatz: Keine Unterscheidung zwischen guten und irrelevanten Aufkommen
     x = np.maximum(val_points-val_points_noise, np.ones(len(y)))
    val = x**y'''

    '''Idee: Vertausche Argumente -> gute Infrastrukturen werden von irrelevanten unterschieden.
    x = np.maximum(val_points-val_points_noise, np.ones(len(y)))
    x = x / np.max(x) + 1
    val = y**x'''

    # create return file
    dffinal = dfmerge_noise
    dffinal["index_right"] = val
    dffinal.columns=["id", "geometry", "val"]
    dffinal['pop'] = x
    dffinal['infra_score'] = y
    print("Resulting Grid:")
    print(dffinal)

    if args[4].endswith(".shp"):
        print("Saving grid to .shp file: " + args[4])
        dffinal.to_file(driver="ESRI Shapefile", filename=args[4])
    elif args[4].endswith(".geojson"):
        print("Saving grid to .geojson file: " + args[4])
        dffinal.to_file(driver="GeoJSON", filename=args[4])
    else:
        print("Unsupported Export File format:" + args[4] + ". Aborting....")
        return
    print("Finished")
    print("--- "+ str(round((time.time() - start_time), 2)) + " seconds ---")
    return
    
if __name__ == "__main__":
    main(sys.argv)
