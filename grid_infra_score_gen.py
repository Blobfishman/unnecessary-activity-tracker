 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 15:28:56 2020
@author: maxi
"""

import pandas as pd
pd.set_option('display.max_rows', 7000)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import geopandas as gpd
from geopandas import GeoDataFrame
from shapely.geometry import Point
from shapely.geometry import Polygon

import numpy as np
import overpy

#       pass
df=pd.read_csv('tagweights.csv', sep=':',header=0)
# Hier Timeout ändern, falls keine Rückgabe!
resultstring = "[out:json][timeout:100];area[name=\"Heidelberg\"]->.searchArea;\n("
crs = {'init': 'epsg:4326'} #4326
#df=pd.read_csv('tagweights.csv', sep=':',header=0)
dict = {}

for i in range(0,len(df.values)):
    for x in range(2, 5):
        string = str(df.values[i][x])
        if not string == "nan":
            liststring = string.split(", ")
            for value in liststring:
                dict[value] = df.values[i][6]
                print(value + " [" + str(df.values[i][6]) + "]")
                for pre in ["node", "way"]:
                    resultstring += "\t" + pre + '[\"' + df.columns[x] + '\"=\"' + value + '\"](area.searchArea);'

resultstring += "); out center;"
#print(resultstring)
api = overpy.Overpass()
#print(dict["cafe"])
result = api.query(resultstring)

print("Finished Query!")
print(str(len(result.nodes)) + " Nodes found.")
print(str(len(result.ways )) + " Ways found.")

#print(dict['bakery'])

def findWeight(tags):
    for poss_key in df.columns[2:5]:
        if poss_key in tags:
            if tags[poss_key] in dict.keys():
                return dict[tags[poss_key]]
    return 0


COLNAMES = ["point_id", "lon", "lat", "infra_score"]
nodedf = pd.DataFrame(columns = COLNAMES)

print("nodes werden analyisert...")
# Nodes
nodes = result.get_nodes()
for i in range(len(nodes)):
    n = nodes[i]
    #print(n.tags)
    weight = findWeight(n.tags)
    nodedf.loc[i] = [n.id, n.lon, n.lat, weight]
    #print(n.id, n.lon, n.lat, weight)

#print(df)

geometry = [Point(xy) for xy in zip(nodedf.lon, nodedf.lat)]

geo_df = GeoDataFrame(nodedf, crs=crs, geometry=geometry)
del geo_df['lat']
del geo_df['lon']
del geo_df['point_id']
print("Infrastrukturwerte der Nodes(Punkte): ")
print(geo_df)
grid = gpd.read_file("gitter_wgs84.shp")
del grid['left']
del grid['top']
del grid['bottom']
del grid['right']
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


#Berechnung der Infrastrukturwerte der ways
print("ways(Umrisse) werden analysiert...")
ways = result.get_ways()
weight_list = [0] * len(ways)
meta_geo = [None] * len(ways)
for i in range(len(ways)):
    w = ways[i]

    weight = findWeight(w.tags)
    weight_list[i] = weight

    nodes = w.get_nodes(resolve_missing=True)
    N = len(nodes)
    lat_list = [0] * N
    lon_list = [0] * N

    for j in range(N):
        n = nodes[j]
        lat_list[j] = n.lat
        lon_list[j] = n.lon

        #print(n.id, n.lon, n.lat,weight)
    poly_geo = Polygon(zip(lon_list, lat_list))
    meta_geo[i] = poly_geo

weight_df = pd.DataFrame(weight_list, columns=["infra_score"])

poly_df = GeoDataFrame(weight_df, crs=crs, geometry=meta_geo)
print("Infrastrukturfaktoren der Umrisse: ")
print(weight_df)

polysjoin = gpd.sjoin(grid, poly_df, how="left", op='contains')
polypivot = pd.pivot_table(polysjoin, index="id", aggfunc={'infra_score': np.sum})
# dfpivot.columns = dfpivot.columns.droplevel()

finalgrid = grid.merge(polypivot, how='left', on="id")
#print(finalgrid)
finalgrid["infra_score"] = np.array(finalgrid["infra_score"]) + np.array(dfgridnew["infra_score"])
#print("Gebäudeumrisse mit Infrastrukturwert: ")
#print(poly)
finalgrid.to_file(driver="ESRI Shapefile", filename='infra_scores.shp')
print("Finales Grid")
print(finalgrid)
print("Finished")
