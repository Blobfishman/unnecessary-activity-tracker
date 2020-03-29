 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 15:28:56 2020
@author: maxi
"""

import pandas as pd
pd.set_option('display.max_rows', 20)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import geopandas as gpd
from geopandas import GeoDataFrame
from shapely.geometry import Point
from shapely.geometry import Polygon
import numpy as np
import overpy
import sys
import time

### Used for displaying Progressbars
# Shoutout to eusoubrasileiro on Stackoverflow!
def progressbar(it, prefix="", size=60, file=sys.stdout):
    count = len(it)
    def show(j):
        x = int(size*j/count)
        file.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, count))
        file.flush()
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    file.write("\n")
    file.flush()

#Main method of generator
def main(args):
    if not (len(args) == 4 or len(args) == 7):
        print("Missing arguments. Usage:")
        print("python3 grid_infra_score_gen.py GRID.shp TAG_WEIGHTS.csv INFRA_SCORES.shp (BLUR_RADIUS GRID_WIDTH GRID_HEIGHT)")
        print("(BLUR RADIUS GRID_WIDTH GRID_HEIGHT) is optional, if you want to blur the outcoming infrastructure values over the grid.")
        return
    if len(args) == 7:
        try:
            blur_radius = float(args[4])
        except:
            print("Invalid blur radius. Input a float! Usage: ")
            print("python3 grid_infra_score_gen.py GRID.shp TAG_WEIGHTS.csv INFRA_SCORES.shp (BLUR_RADIUS GRID_WIDTH GRID_HEIGHT)")
            print("(BLUR RADIUS GRID_WIDTH GRID_HEIGHT) is optional, if you want to blur the outcoming infrastructure values over the grid.")
            return
        try:
            grid_width = int(args[5])
            grid_height = int(args[6])
        except:
            print("Invalid grid dimensions GRID_WIDTH and GRID_HEIGHT. Input Integers! Usage:")
            print("python3 grid_infra_score_gen.py GRID.shp TAG_WEIGHTS.csv INFRA_SCORES.shp (BLUR_RADIUS GRID_WIDTH GRID_HEIGHT)")
            print("(BLUR RADIUS GRID_WIDTH GRID_HEIGHT) is optional, if you want to blur the outcoming infrastructure values over the grid.")
            return

        if blur_radius <= 0 or blur_radius >= 1:
            print("Blur radius should be between 0 and 1. Aborting..")
            return

        if grid_width <= 1 or grid_height <= 1:
            print("Grid dimensions should be greater than 1. Aborting..")
            return

    start_time = time.time()

    print("Loading Grid Shape File at " + args[1])

    try:
        grid = gpd.read_file(args[1])
        minx, miny, maxx, maxy = grid.geometry.total_bounds
        rnd = 4
        bb = str(round(miny,rnd)) + "," + str(round(minx,rnd)) + "," + str(round(maxy,rnd)) + "," + str(round(maxx,rnd))
        print("Found boundary box (", bb,  ")")
        try:
            del grid['left']
            del grid['top']
            del grid['bottom']
            del grid['right']
        except:
            pass
    except:
        print("Grid shape file could not be located or is in the wrong format. Do you also have the .shx and .dbf files in the same directory?. Aborting..")
        return

    if not grid_width*grid_height==len(grid['geometry']):
        print("Given grid dimensons do not match given grid file.")
        print(str(grid_width)+"*"+str(grid_height)+"="+str(grid_width*grid_height)+"≠"+str(len(grid["geometry"])))
        print("Aborting..")
        return

    print("Loading Tagweights at " + args[2])
    try:
        df=pd.read_csv(args[2], sep=':',header=0)
    except:
        print("Tagweights file could not be located. Aborting..")
        return


    # Hier Timeout ändern, falls keine Rückgabe!
    #resultstring = "[out:json][timeout:200];area[name=\"Heidelberg\"]->.searchArea;\n("


    resultstring = "[out:json][timeout:1000];\n("
    crs = {'init': 'epsg:4326'} #4326
    #df=pd.read_csv('tagweights.csv', sep=':',header=0)
    dict = {}
    print("Collecting dict: ")
    n_cols = len(df.columns)-3
    for i in range(0,len(df.values)):
        for x in range(2, n_cols+2):
            string = str(df.values[i][x])
            if not string == "nan":
                liststring = string.split(", ")
                for value in liststring:
                    dict[value] = df.values[i][n_cols+2]
                    print(df.columns[x] + ":" + value + " [" + str(df.values[i][n_cols+2]) + "]")
                    for pre in ["node", "way"]:
                        resultstring += "\t" + pre + '[\"' + df.columns[x] + '\"=\"' + value + '\"](' + bb +');'

    resultstring += ");out; >; out;"
    #print(resultstring)
    api = overpy.Overpass()
    #print(dict["cafe"])
    print("Starting query. This may take a while(<1000s)...")
    try:
        result = api.query(resultstring)
    except:
        print("Query failed. Either there is no internet connection or the query timed out. Aborting..")
        return


    print("Finished query!")
    print(str(len(result.nodes)) + " nodes found.(also containing the unzipped ways)")
    print(str(len(result.ways )) + " ways found.")

    if len(result.nodes) == 0 and len(result.ways) == 0:
        print("Query returned no results. Maybe it took to long. Aborting...")
        return
    #print(dict['bakery'])

    def findWeight(tags):
        for poss_key in df.columns[2:5]:
            if poss_key in tags:
                if tags[poss_key] in dict.keys():
                    return dict[tags[poss_key]]
        return 0


    COLNAMES = ["point_id", "lon", "lat", "infra_score"]
    nodedf = pd.DataFrame(columns = COLNAMES)

    #### F+r Fortschrittsleiste
    # Try to get Terminal width
    try:
        _, terminal_columns = os.popen('stty size', 'r').read().split()
    except:
        terminal_columns = 60
    print("Analyzing nodes...")
    # Nodes
    nodes = result.get_nodes()
    ct = 0
    for i in progressbar(range(len(nodes)), size=terminal_columns):

        n = nodes[i]
        #print(n.tags)
        weight = findWeight(n.tags)
        nodedf.loc[ct] = [n.id, n.lon, n.lat, weight]
        ct+=1
        #print(n.id, n.lon, n.lat, weight)

    #Berechnung der Infrastrukturwerte der ways
    print("Analyzing ways(e.g. building outlines)...")
    #ct = len(nodedf.loc)-1
    ways = result.get_ways()
    for i in progressbar(range(len(ways)),size=terminal_columns):
        w = ways[i]

        nodes = w.get_nodes(resolve_missing=False)
        N = len(nodes)
        weight = findWeight(w.tags)
        for j in range(N):
            n = nodes[j]
            nodedf.loc[ct] = [n.id, n.lon, n.lat, weight]
            ct+=1

    #print(df)

    geometry = [Point(xy) for xy in zip(nodedf.lon, nodedf.lat)]

    geo_df = GeoDataFrame(nodedf, crs=crs, geometry=geometry)
    del geo_df['lat']
    del geo_df['lon']
    del geo_df['point_id']

    #print("Eingangs- Grid: ")
    #print(grid)

    dfsjoin = gpd.sjoin(grid, geo_df, how="left", op='contains')
    dfpivot = pd.pivot_table(dfsjoin, index="id", aggfunc={'infra_score': np.sum})
    # dfpivot.columns = dfpivot.columns.droplevel()

    #print(dfpivot)

    dfgridnew = grid.merge(dfpivot, how='left', on="id")
    #dfgridnew.drop(["left", "top", "right", "bottom"], axis=1)
    #print("Gemergtes Grid: ")
    #print(dfgridnew)

    #id, infra_wert, geometry


    # ways = result.get_ways()
    # weight_list = [0] * len(ways) * 200
    # meta_geo = [None] * len(ways) * 200
    # ct = 0
    # for i in range(len(ways)):
    #     w = ways[i]
    #     #print(w.tags)
    #
    #     nodes = w.get_nodes(resolve_missing=False)
    #     N = len(nodes)
    #     lat_list = [0] * N
    #     lon_list = [0] * N

        # Scheint nicht zu funktionieren -> siehe Discord
        # for j in range(N):
        #     n = nodes[j]
        #     lat_list[j] = n.lat
        #     lon_list[j] = n.lon
        #
        #
        #     #print(n.id, n.lon, n.lat,weight)
        # try:
        #     poly_geo = Polygon(zip(lon_list, lat_list))
        #     meta_geo[ct] = poly_geo
        #     weight_list[ct] = findWeight(w.tags)
        #     ct+=1
        # except:
        #     print("WARNING: Invalid way-shape fount @id=", w.id, "i=", i , " ignoring...")
        #     pass

        #Alternative Eckpunkte mit Bewertung/Anzahl der Eckpunkte
        # weight_w = findWeight(w.tags)
        # for j in range(N):
        #     n = nodes[j]
        #     lat_list[j] = n.lat
        #     lon_list[j] = n.lon
        #
        #     point_geo = Point(n.lat, n.lon)
        #     meta_geo[ct] = point_geo
        #     weight_list[ct] = weight_w/N
        #     ct+=1

    #arrs zuschneiden
    # weight_list = weight_list[:ct]
    # meta_geo = meta_geo[:ct]
    #
    # weight_df = pd.DataFrame(weight_list, columns=["infra_score"])
    #
    # poly_df = GeoDataFrame(weight_df, crs=crs, geometry=meta_geo)
    # # print("infra_sco")
    # # print(weight_df)
    #
    # polysjoin = gpd.sjoin(grid, poly_df, how="left", op='contains')
    # polypivot = pd.pivot_table(polysjoin, index="id", aggfunc={'infra_score': np.sum})
    # # dfpivot.columns = dfpivot.columns.droplevel()
    #
    # finalgrid = grid.merge(polypivot, how='left', on="id")
    # #print(finalgrid)
    # insc_ways = np.array(finalgrid["infra_score"])

    #BLUR
    if len(args) == 5:

        # Get Grid width and height
        N = len(dfgridnew["geometry"])
        #TODO: Width und Height Findung richtig implementieren

        # Die Lösung hier funktioniert nicht, da die Gridrechtecke da neben-
        # einanderliegende Kästchen in x-richtung andere y-Werte haben(wegen
        # Erdkrümmung)
        # _, miny, _, _ = dfgridnew["geometry"][0].bounds
        # width = 1
        # print(miny)
        # for i in range(1, N):
        #
        #     #print(dfgridnew["geometry"][i].exterior.coords[0][1])
        #     _, y, _, _ = dfgridnew["geometry"][i].bounds
        #     print(i, y, miny == round(miny,4) == np.round(y,4))
        #     if miny == y:
        #         width+=1

        # height = np.floor(N/weight)

        # Die Lösung hier ist sehr WIP und im Allgemeinen falsch
        #, da nicht alle grids quadratisch sind, aber in unserem Fall funtkioniert
        # sie
        #width = int(np.floor(np.sqrt(N)))
        #height= int( np.ceil(np.sqrt(N)))

        # Mögliche Lösungen:
        # I  : Der user gibt width und height an
        # II : Es gibt einen intelligenteren Weg width und height herauszufinden
        # III: Es gibt einen intelligenteren Weg zu bluren

        # Lösung I
        width = grid_width
        height = grid_height
        print("Found grid: (width:", width ,"height:", height,")")
        new_infra = np.array(dfgridnew["infra_score"])
        print("Bluring the grid with blur radius", blur_radius)
        for i in progressbar(range(N)):
            blur_val = dfgridnew["infra_score"][i]*blur_radius
            new_infra[j+1] +=  blur_val
            new_infra[j-1] +=  blur_val
            new_infra[j+width] +=  blur_val
            new_infra[j-width] +=  blur_val

        #apply blur
        dfgridnew["infra_score"] = new_infra



    insc = - np.array(dfgridnew["infra_score"])
    #print(insc_nodes)
    ct = 0
    for x in insc:
        if x==0:
            ct+=1
    print("Amount of unaffected cells: " + str(ct) + "/" + str(len(grid)) +
            "  ( " + str(round(ct/len(grid)*100,2)) + "% )")
    print("Distribution")
    # print("Nodes:")
    print("\tµ = ", np.mean(insc))
    print("\tσ = ", np.std (insc))
    # print("Ways:")
    # print("\tµ = ", np.mean(insc_ways))
    # print("\tσ = ", np.std (insc_ways))
    #normalisieren


    #Normalisieren
    '''--Problem mit Wurzeln: Ich hasse Wurzeln
    fac = 0.1
    sgn = np.sign(insc)
    insc = sgn*(sgn*insc/np.max(np.abs(insc)))**fac + 1
    dfgridnew["infra_score"] =  insc'''

    '''--Problem mit Extremwerten
    insc = -insc_nodes - insc_ways
    insc = insc/np.max(np.abs(insc)) + 1
    dfgridnew["infra_score"] =  insc'''

    #Idee: Log. skalieren
    insc = insc + abs(np.min(insc)) + 1
    insc = np.log(insc) + 1
    dfgridnew["infra_score"] =  insc

    #Nochmal ausgeben
    print("After applying a normilization function:" )
    print("\tµ = ", np.mean(insc))
    print("\tσ = ", np.std (insc))
    #print("Gebäudeumrisse mit Infrastrukturwert: ")
    #print(poly)

    print("Generated grid:")
    print(dfgridnew)
    dfgridnew.to_file(driver="ESRI Shapefile", filename=args[3])
    print("Saved grid to file: " + args[3])
    print("Finished")
    print("--- "+ str(round((time.time() - start_time), 2)) + " seconds ---")
    return

if __name__ == "__main__":
    main(sys.argv)
