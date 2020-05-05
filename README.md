# unnecessary-activity-tracker

A tool to identify and visualize necessary and unnecessary movement of people using positional data(e. g. anonymized cell phone data)

This project was made in 48h hours _from 20.3 to 22.3.20_ during the [#WirVsVirus](http://www.wirvsvirushackathon.org) -Hackathon, to contribute in finding solutions to the ongoing corona-crisis.
Also see the corresponding [devpost-page](https://devpost.com/software/0045_haustiere_handydaten)(in German).
Feel free to pull request and help continue the development of this tool!

### Dependencies
You'll need Python3.x and the following modules:
_overpy, numpy, pandas, geopandas, rtree_
You can install them using pip or anaconda!

## Generating the Infrastructure.
A Shape File of the grid(the _.shp, .dbf, .shx_ files), and a Table with tagweights.
(You could also use a single _.geojson_ files in all commands instead of _.shp_ files)

With the files in the same folder you may generate the infrastructure score for all grids by running
```
python3 grid_infra_score_gen.py GRID.shp TAG_WEIGHTS.csv INFRA_SCORES.shp
```
Keep in mind that this generates 5 filetypes: _.cpg, .dbf, .prj, .shp, .shx_
If you set this Tool up on a server, you might consider running this very rarely, because the infra_score won't change much over time.

(You can run this command with optional arguments `blur_radius`(between 0 and 1), `grid_width` and `grid height`(each >= 1, numbers of cells in each direction) at the end. This will blur the out coming infrastructure values to more realistically simulate the movement of persons between tiles. Keep in mind that this feature is work in progress!)


## Generating the Shape File
If you got your files(especially the _.shp and .shx_) you can generate the final scores by executing
```
python3 main.py INFRA_SCORES.shp POINTS.shp POINTS_BACKGROUND.shp SCORES_FINAL.shp
```
Where `POINTS.shp` is the Shape File for the Location-Data, and `POINTS_BACKGROUND.shp` is the background noise of given Point distribution.(e. g. at night)
`SCORES_FINAL.shp` is the output Shape File.
This generates 5 filetypes(_.cpg, .dbf, .prj, .shp, .shx_) for the Shape file with values
* **geometry** The Grids geometries to display in any compatible program
* **val** The calculated value of unnecessary-movement calculated with _val(pop, infra)= pop^infra + pop_
* **pop** The calculated value of movement
* **infra** The generated values for the surrounding infra structure
By specifying a different file extension for the output file to _.geojson_, a .geojson with the given data ready to display with web-frontend APIs(e. g. Leaflet, OpenLayers)


## No files? We have example data!
We ran the analysis for Heidelberg, Germany.
You can use the grid and point files found in `data-samples`.

Feel also free to use the hand-crafted `tagweights.csv`.
You can try it out yourself! Just run
```
python3 grid_infra_score_gen.py data-samples/grid/grid.shp tagweights.csv infra_score.shp
```
and after that
```
python3 main.py infra_score.shp data-samples/points/points.shp data-samples/t0_points/t0_points.shp final_score.shp
```
You should now be able to open the generated Shape File in any compatible editor.(e. g [QGis](https://www.qgis.org/de/site/index.html))
Alternatively you can see the result for Heidelberg [here](http://wirvsvirus.lpk-server.de/).


## Built With

* [Overpass](https://github.com/drolbr/Overpass-API) - The API used to access OSM
* [overpy](https://github.com/DinoTools/python-overpy) - Python Wrapper to access the Overpass-API
* [geopandas](https://geopandas.org/) - Mapping Python Library build on pandas to create and import Shape Files
* [QGis](https://www.qgis.org/de/site/index.html) - Used to display Shape Files
* [qgis2web](https://github.com/tomchadwin/qgis2web) - Used to generate the [example page](http://wirvsvirus.lp-kb.eu/)
* [OpenLayers](https://openlayers.org/) - Used by qgis2web to display the map
