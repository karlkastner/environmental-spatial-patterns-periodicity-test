#!/usr/bin/python3
# Mon 25 Oct 15:08:44 CEST 2021
# Karl Kastner
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 
## fetch vegetation patterns from the google satellite server

import matplotlib.pyplot
import matplotlib.image
from   qgis.core import *
from   qgis.core import QgsVectorLayer,QgsProject,QgsApplication,QgsRasterLayer
from   PyQt5.QtCore import *
from   PyQt5.QtGui import *
import qgis.utils
import time
import os
import math
import scipy.io
import requests 
import tempfile
import sys

show = False

# 0 is script name
pattern_type = sys.argv[1];
dx_str = sys.argv[2];
dx     = float(dx_str);

dmin = 2e-5;

dir_str    = './'

pattern_type_C = ["anisotropic", "isotropic"];

# initialize QGIS
qgs = QgsApplication([], False)
qgs.initQgis()
	
# create a project
project = QgsProject.instance() 

dy      = dx
	
# manually delineated polygons enclosing the patterns
input_str       = dir_str + '/mat/vegetation-patterns-' + pattern_type + '-sampling-interval.shp'
shp_reprojected = 'mat/vegetation-patterns-polygons-' + pattern_type + '-reprojected.shp'

# output qgis file
project_str  = dir_str + '/mat/vegetation-patterns-' + pattern_type + '.qgz'

# output directory names for fetched images
img_dir  = dir_str + '/img/' + pattern_type + '/' + dx_str + "/"
# output filename base for cropped tiles
img_base = 'google-satellite_'
# output file name for sampled transects in matlab format
#mat_str      = dir_str + '/mat/vegetation-patterns-' + pattern_type + '-sampled.mat'

# google_url of google maps tile server
google_url = "mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}" 
google_url = requests.utils.quote(google_url)
google_url = 'type=xyz&url=https://'+google_url;

# create output folders
folder_a = ["/mat", "/dat", "/img/" + pattern_type + '/', '/img/' + pattern_type + '/' + dx_str]
for folder_str in folder_a:
	if (not os.path.exists(dir_str + "/" +folder_str)):
		os.mkdir(dir_str + "/" + folder_str)
	# end if
# end for d in a

# lat-long wgs4
shp_crs_str = "EPSG:4326"
# pseudo-mercator
google_crs_str = "EPSG:900913"

shp_crs  = QgsCoordinateReferenceSystem(shp_crs_str)
google_crs = QgsCoordinateReferenceSystem(google_crs_str)
	
cmd = 'ogr2ogr ' + ' -s_srs ' + shp_crs_str + ' -t_srs ' + google_crs_str + ' ' + shp_reprojected + ' ' + input_str
os.system(cmd);

# add google satellite tile server as source layer
google_layer = QgsRasterLayer(google_url, 'google-satellite', 'wms')  
if google_layer.isValid():
	QgsProject.instance().addMapLayer(google_layer)
else:
	raise Exception('Connection to google-tile server failed');
# end if
# alternatively, load it from the project
# google_layer   = QgsProject.instance().mapLayersByName("google-satellite")[0]

# open the input shapefile
transect_layer = QgsVectorLayer(input_str, 'borders', 'ogr')
if not transect_layer.isValid():
	raise Exception('Layer is invalid')

# read input shape coordinates
feature_a   = transect_layer.getFeatures()
bb_a        = []
dx_sample_a = []
xy          = []
xyc         = []
id_a        = [];
region_a    = [];
length_a    = [];
continent_a = [];
name_a      = [];
i = 0;
print('Reading polygon coordinates');
print(input_str)
for feature in feature_a:
	# todo call by name
	# continent_a.append(feature.attributes()[2]);
	bb = feature.geometry().boundingBox()
	xyc.append(feature.geometry().centroid().asPoint())

	#field_names = [field.name() for field in transect_layer.fields()]
	#print(field_names)

	dx_sample_a.append(feature['dx_sample']);
	#%wait = input("Press Enter to continue.")

	name_str = (img_base + "%+09.5f" % xyc[i][1]) + "_" + ("%+010.5f" % xyc[i][0]);
	name_a.append(name_str);

	# transform to the coordinate system used by the tile server
	transform  = QgsCoordinateTransform(shp_crs, google_crs, QgsProject.instance())
	bb         = transform.transformBoundingBox(bb)
	bb_a.append(bb)
	i = i+1
# end for feature

print('Fetching google-satellite');
options = QgsMapSettings()
options.setLayers([google_layer])
options.setBackgroundColor(QColor(0, 255, 0))
for i in range(len(bb_a)):
	bb = bb_a[i]
	oname = (img_dir + '/' + name_a[i] + '.png')
	print(str(i) + ' ' + dx_str + ' ' + str(dx_sample_a[i]) + ' ' + oname)
	if ((dx_sample_a[i] <= dx) and (dx_sample_a[i]>0) and (not os.path.exists(oname))):
		print(str(i) + " " + name_a[i])

		lx = round((bb.xMaximum() - bb.xMinimum())/dx)
		ly = round((bb.yMaximum() - bb.yMinimum())/dy)
		size    = QSize(lx,ly)

		options.setDestinationCrs(QgsCoordinateReferenceSystem("EPSG:900913"))
		options.setOutputSize(size)
		options.setExtent(bb)
		
		# render image
		img = QImage(size, QImage.Format_ARGB32_Premultiplied)
		painter = QPainter()
		painter.begin(img)
		#p.setRenderHint(QPainter.Antialiasing)
		render = QgsMapRendererCustomPainterJob(options, painter)
		
		render.start()
		render.waitForFinished()
		painter.end()
	
		# save image
		img.save(oname)
		pgw = QgsMapSettingsUtils.worldFileContent(options)
		pgw_str = (img_dir + '/' + name_a[i] + '.pgw')
		#pgw_str = (img_str + id_str + '.pgw')
		with open(pgw_str, "w") as f:
		    f.write(pgw)
	# end if not exists oname

	mskfile = img_dir + '/' + name_a[i] + '_mask.tif'
	if ((dx_sample_a[i] <= dx) and (dx_sample_a[i]>0) and (not os.path.exists(mskfile))):
		# create a mask file
		cmd= ('gdalwarp -s_srs EPSG:3857 -of Gtiff ' + oname + '  ' + mskfile)
		os.system(cmd)
		# zero mask file everywhere
		cmd= ('gdal_calc.py --overwrite --co COMPRESS=LZW -A ' + mskfile + ' --outfile=' + mskfile + ' --calc=0')
		os.system(cmd)
		# burn pattern mask
		cmd=('gdal_rasterize -burn 1 ' + "-where \"id\"='" + str(i) + "' " + shp_reprojected + ' ' + mskfile)
		os.system(cmd)
	# end if not exists mskfile
# end for i

# save project
project.write(project_str);

# there is a segfault when closing Qgis for an unknown reason,
# this does not affect the analysis
qgs.exitQgis();

