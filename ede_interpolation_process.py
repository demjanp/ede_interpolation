
from PyQt5 import (QtCore, QtWidgets, QtGui)

from qgis.core import QgsFieldProxyModel, QgsMapLayerProxyModel, QgsRasterLayer, QgsProject, QgsColorRampShader, QgsRasterShader, QgsSingleBandPseudoColorRenderer, QgsColorRamp
from qgis.utils import iface

from .ede_function import *

import os
import numpy as np
from osgeo import gdal, osr
from collections import defaultdict

def load_curve(fcalib):
	# load calibration curve
	# data from: fcalib 14c file
	# returns: [[CalBP, ConvBP, CalSigma], ...], sorted by CalBP
	
	with open(fcalib, "r", encoding = "ansi") as f:
		data = f.read()
	data = data.split("\n")
	cal_curve = []
	for line in data:
		if line.startswith("#"):
			continue
		cal_curve.append([float(value) for value in line.split(",")])
	cal_curve = np.array(cal_curve)
	cal_curve = cal_curve[np.argsort(cal_curve[:,0])]
	
	return cal_curve

def calibrate(age, uncert, curve_conv_age, curve_uncert):
	# calibrate a 14C measurement
	# calibration formula as defined by Bronk Ramsey 2008, doi: 10.1111/j.1475-4754.2008.00394.x
	# age: uncalibrated 14C age BP
	# uncert: 1 sigma uncertainty
	
	sigma_sum = uncert**2 + curve_uncert**2
	return (np.exp(-(age - curve_conv_age)**2 / (2 * sigma_sum)) / np.sqrt(sigma_sum))

def bp_to_ce(t):
	
	t = int(t) - 1950
	if t > 0:
		return t, "BCE"
	return -t, "CE"

class EDEInterpolationProcess(QtWidgets.QProgressDialog):
	
	def __init__(self, data, s_duration, s_diameter, time_step, time_from, time_to, cell_size, approximate, path_layers, path_summed, crs):
		
		self.running = True
		self.s_duration = s_duration
		self.s_diameter = s_diameter
		self.time_step = int(time_step)
		self.time_from = time_from
		self.time_to = time_to
		self.cell_size = int(cell_size)
		self.approximate = approximate
		self.path_layers = path_layers
		self.path_summed = path_summed
		self.crs = crs
		
		lookup_f_s = {}
		
		QtWidgets.QProgressDialog.__init__(self, iface.mainWindow())
		
		self.setWindowTitle("EDE Interpolation")
		self.setMaximum(100)
		self.setValue(0)
		self.canceled.connect(self.on_canceled)
		
		self.show()
		
		self.process(data)
		
		self.close()
	
	def on_canceled(self):
		
		self.running = False
	
	def save_summed(self, ts, summed):
		
		with open(self.path_summed, "w") as f:
			f.write("Date_BP,Date_CE,EDE_summed\n")
			for ti in range(ts.shape[0]):
				t_ce, cebce = bp_to_ce(ts[ti])
				f.write("%d,%d %s,%f\n" % (int(ts[ti]), t_ce, cebce, summed[ti]))
	
	def save_raster(self, grid, x0, y0, path):
		
		driver = gdal.GetDriverByName("GTiff")
		out_raster = driver.Create(
			path,
			grid.shape[1],
			grid.shape[0],
			1,
			gdal.GDT_Float32,
		)
		out_raster.SetGeoTransform((
			x0,
			self.cell_size, 
			0,
			y0,
			0,
			self.cell_size,
		))
		out_raster.SetProjection(self.crs.toWkt())
		out_raster.GetRasterBand(1).WriteArray(grid)
		out_raster = None
	
	def process(self, data):
		
		# data = {upd/npd: {dating = [calendar years BP, ...], uncert = [calendar years, ...], coords = [[x, y], ...], accur = [accuracy, ...]}}
		
		UPD_t_ds = np.round(data["upd"]["dating"]).astype(int)  # mean datings of archaeological components (calendar years BP)
		UPD_uncert_ds = np.round(data["upd"]["uncert"]).astype(int)  # uncertainties of the datings (calendar years)
		UPD_As = np.round(data["upd"]["coords"]).astype(int)  # spatial coordinates of archaeological components (metres)
		UPD_accurs = np.round(data["upd"]["accur"]).astype(int)  # accuracies of spatial coordinates of archaeological components (+-metres)
		
		NPD_t_ds = np.round(data["npd"]["dating"]).astype(int)  # measured radiocarbon ages of archaeological components (radiocarbon years)
		NPD_uncert_ds = np.round(data["npd"]["uncert"]).astype(int)  # 1-sigma uncertainties of the measured radiocarbon ages (radiocarbon years)
		NPD_As = np.round(data["npd"]["coords"]).astype(int)  # spatial coordinates of archaeological components (metres)
		NPD_accurs = np.round(data["npd"]["accur"]).astype(int)  # accuracies of spatial coordinates of archaeological components (+-metres)
		
		if (not UPD_t_ds.size) and (not NPD_t_ds.size):
			return
		
		s_halflife = self.s_duration / 2  # expected half-life of a settlement in years
		s_radius = self.s_diameter / 2 # expected radius of a settlement in metres
		
		# temporal extent
		t_min = np.inf
		t_max = -np.inf
		if UPD_t_ds.size:
			t_min = min(t_min, (UPD_t_ds - UPD_uncert_ds).min() - self.s_duration)
			t_max = max(t_max, (UPD_t_ds + UPD_uncert_ds).max() + self.s_duration)
		if NPD_t_ds.size:
			t_min = min(t_min, (NPD_t_ds - 2*NPD_uncert_ds).min() - self.s_duration)
			t_max = max(t_max, (NPD_t_ds + 2*NPD_uncert_ds).max() + self.s_duration)
		t_min, t_max = [int(round(value / 10) * 10) for value in [t_min, t_max]]
		
		if self.time_from is not None:
			t_max = min(t_max, self.time_from)
		if self.time_to is not None:
			t_min = max(t_min, self.time_to)
		
		ts_slices = np.arange(t_max, t_min - 2*self.time_step, -self.time_step).tolist()  # times of time slices
		
		# prepare lookup for probability distributions of 14C datings
		self.setLabelText("Calibrating radiocarbon dates")
		cal_curve = load_curve(os.path.join(os.path.dirname(__file__), "intcal13.14c")) # [[CalBP, ConvBP, CalSigma], ...], sorted by CalBP
		
		# filter calibration curve to include only time-step dates
		cal_curve = cal_curve[(cal_curve[:,0] >= t_min) & (cal_curve[:,0] <= t_max)][::-1]
		ts = cal_curve[:,0]
		curve_conv_age = cal_curve[:,1]
		curve_uncert = cal_curve[:,2]
		if ts[-1] < ts_slices[-1]:
			ts_slices.append(ts[-1])
		
		# calculate probability distributions for all combinations of 14c age and uncertainty
		unique_dates = set()  # ((age, uncert), ...)
		for idx in range(NPD_t_ds.shape[0]):
			unique_dates.add((NPD_t_ds[idx], NPD_uncert_ds[idx]))
		lookup_14c = defaultdict(dict)  # {age: {uncert: D, ...}, ...}; D[ti] = p; where ti = index in ts, p = probability
		cmax = len(unique_dates)
		cnt = 0
		for age, uncert in unique_dates:
			QtWidgets.QApplication.processEvents()
			if not self.running:
				return
			self.setValue((cnt / cmax) * 100)
			cnt += 1
			lookup_14c[age][uncert] = calibrate(age, uncert, curve_conv_age, curve_uncert)
		
		# prepare lookup of spatial probability distribution around evidence points
		self.setLabelText("Calculating spatial probability distribution")
		self.setValue(0)
		accurs = set()
		accurs.update(UPD_accurs.tolist())
		accurs.update(NPD_accurs.tolist())
		lookup_f_s = {} # {accur: M, ...}; M[n, n] = f_s(d, accur, s_radius); where center is point A and n is 2 * [maximum distance from A in raster units] + 1; where f_s > 0
		cnt = 0
		cmax = len(accurs)
		for accur in accurs:
			QtWidgets.QApplication.processEvents()
			if not self.running:
				return
			self.setValue((cnt / cmax) * 100)
			cnt += 1
			r = int(round((accur + 2*s_radius) / self.cell_size))
			n = 2 * r + 1
			lookup_f_s[accur] = np.zeros((n, n), dtype = float)
			rcs = np.argwhere(np.ones((n,n), dtype = bool))
			mask = (rcs > r).all(axis = 1)
			for row, col in rcs[mask]:
				d = (((row - r)**2+(col - r)**2)**0.5)*self.cell_size
				if self.approximate:
					p = f_s_approx(d, accur, s_radius)
				else:
					p = f_S_lens(d, accur, s_radius) / f_S(accur, s_radius)
				if (p == np.inf) or np.isnan(p):
					p = 0
				lookup_f_s[accur][row,col] = p
				lookup_f_s[accur][n - row,col] = p
				lookup_f_s[accur][row, n - col] = p
				lookup_f_s[accur][n - row, n - col] = p
			lookup_f_s[accur][0,0] = lookup_f_s[accur][0,1]
		
		# spatial extent
		row_min, col_min = np.inf, np.inf
		row_max, col_max = -np.inf, -np.inf
		for As, accurs in [[UPD_As, UPD_accurs], [NPD_As, NPD_accurs]]:
			for idx in range(As.shape[0]):
				A = As[idx]
				accur = accurs[idx]
				r = int(lookup_f_s[accur].shape[0] / 2)
				col, row = np.round(A / self.cell_size).astype(int)
				row_min = min(row_min, row - r - 1)
				col_min = min(col_min, col - r - 1)
				row_max = max(row_max, row + r)
				col_max = max(col_max, col + r)
		width, height = (col_max - col_min), (row_max - row_min)
		x0, y0 = col_min * self.cell_size, row_min * self.cell_size
		
		# calculate time-slices
		self.setLabelText("Generating time-slices")
		paths = []
		summed = []
		val_max = -np.inf
		grid_summed = np.zeros((height, width), dtype = float)
		t_slice = ts_slices.pop(0)
		n_slice = 1
		for ti in range(ts.shape[0]):
			QtWidgets.QApplication.processEvents()
			if not self.running:
				return
			self.setValue((ti / ts.shape[0]) * 100)
			
			grid = np.ones((height, width), dtype = float)
			
			for idx in range(UPD_t_ds.shape[0]):
				t_d = UPD_t_ds[idx]
				uncert_d = UPD_uncert_ds[idx]
				A = UPD_As[idx]
				accur = UPD_accurs[idx]
				M = 1 - lookup_f_s[accur] * f_t_UPD(ts[ti], t_d, uncert_d, s_halflife)
				r = int((M.shape[0] - 1) / 2)
				col0, row0 = np.round((A - [x0, y0]) / self.cell_size - r - 1).astype(int)
				grid[row0:row0 + M.shape[0],col0:col0 + M.shape[0]] *= M
			
			for idx in range(NPD_t_ds.shape[0]):
				t_d = NPD_t_ds[idx]
				uncert_d = NPD_uncert_ds[idx]
				A = NPD_As[idx]
				accur = NPD_accurs[idx]
				M = 1 - lookup_f_s[accur] * f_t_NPD(ts[ti], s_halflife, lookup_14c[t_d][uncert_d], ts)
				r = int((M.shape[0] - 1) / 2)
				col0, row0 = np.round((A - [x0, y0]) / self.cell_size - r - 1).astype(int)
				grid[row0:row0 + M.shape[0],col0:col0 + M.shape[0]] *= M
			
			grid = 1 - grid
			grid[np.isnan(grid)] = 0
			grid[grid == np.inf] = 0
			
			summed.append(grid.sum())
			
			if ts[ti] <= t_slice:
				val_max = max(val_max, grid_summed.max())
				t_ce, cebce = bp_to_ce(t_slice)
				t_ce2, cebce2 = bp_to_ce(ts_slices[0])
				datestr = "%03d_%d_%s_-_%d_%s" % (n_slice, t_ce, cebce, t_ce2, cebce2)
				paths.append([datestr, os.path.join(self.path_layers, "ede_%s.tif" % (datestr))])
				self.save_raster(grid_summed, x0, y0, paths[-1][1])
				t_slice = ts_slices.pop(0)
				n_slice += 1
				grid_summed[:] = grid
			else:
				grid_summed += grid
		
		if self.path_summed:
			self.save_summed(ts, summed)
		
		project = QgsProject.instance()
		val_max = val_max*0.9
		for datestr, path in paths:
			layer = QgsRasterLayer(path, "EDE_%s" % (datestr))
			layer.setCrs(self.crs)
			
			s = QgsRasterShader()
			c = QgsColorRampShader()
			c.setColorRampType(QgsColorRampShader.Interpolated)
			i = [] 
			i.append(QgsColorRampShader.ColorRampItem(0, QtGui.QColor('#ffffff')))
			i.append(QgsColorRampShader.ColorRampItem(val_max, QtGui.QColor('#000000')))
			c.setColorRampItemList(i)
			s.setRasterShaderFunction(c)
			ps = QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 1, s)
			ps.setClassificationMin(0)
			ps.setClassificationMax(val_max)
			layer.setRenderer(ps)
			
			project.addMapLayer(layer)
