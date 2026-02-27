# SUMMARY:      rowcolmap.py (QGIS version)
# USAGE:        Creates a polygon representation of a raster grid in QGIS.
# ORG:          Pacific Northwest National Laboratory & Duke University
# AUTHOR:       Adapted for QGIS by Y.Song
# E-MAIL:       yiyun.song@duke.edu
# ORIG-DATE:    Apr-2017
# LAST CHANGE:  2025-09-09
# DESCRIPTION:  Converts ArcGIS `arcpy` script to work with QGIS `QgsRasterLayer`

from qgis.core import (
    QgsRasterLayer, QgsVectorLayer, QgsField, QgsFeature, QgsGeometry, QgsPointXY, QgsProcessing
)
from PyQt5.QtCore import QVariant
import processing

def rowcolmapfun(elev_path, output_polygon):
    """
    Converts DEM raster into a polygon representation where each polygon represents a DEM cell.
    Each polygon is assigned row and column numbers.

    Parameters:
    elev_path (str): Path to the DEM raster file (.tif)
    output_polygon (str): Output shapefile path for polygon representation
    """

    # Load DEM raster
    dem_raster = QgsRasterLayer(elev_path, "dem_raster")

    if not dem_raster.isValid():
        print(f"Error: Unable to load DEM file {elev_path}")
        return

    # Convert raster to polygon
    params = {
        'INPUT': elev_path,
        'BAND': 1,
        'FIELD': "VALUE",
        'OUTPUT': output_polygon
    }
    processing.run("gdal:polygonize", params)
    print(f"Raster to polygon conversion completed: {output_polygon}")

    # Load polygon layer
    polygon_layer = QgsVectorLayer(output_polygon, "rowcolpoly", "ogr")
    
    if not polygon_layer.isValid():
        print(f"Error: Unable to load polygon layer {output_polygon}")
        return
    
    # Add "ROW" and "COL" fields
    polygon_layer.startEditing()
    polygon_layer.dataProvider().addAttributes([
        QgsField("ROW", QVariant.Int),
        QgsField("COL", QVariant.Int)
    ])
    polygon_layer.updateFields()

    # Get raster resolution (cell size)
    transform = dem_raster.dataProvider().transform()
    cell_size_x = abs(transform.scaleX())  # X resolution
    cell_size_y = abs(transform.scaleY())  # Y resolution

    # Assign row and column values
    for feature in polygon_layer.getFeatures():
        centroid = feature.geometry().centroid().asPoint()
        row = int((centroid.y() - dem_raster.extent().yMinimum()) / cell_size_y)
        col = int((centroid.x() - dem_raster.extent().xMinimum()) / cell_size_x)

        feature.setAttribute("ROW", row)
        feature.setAttribute("COL", col)
        polygon_layer.updateFeature(feature)

    # Save changes
    polygon_layer.commitChanges()
    print(f"Row and Column mapping completed: {output_polygon}")

