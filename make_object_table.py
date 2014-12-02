import os
import glob
from astropy.table import Table
from astropy.io import fits
import numpy as np
from scipy.stats import mode

import pdb

global REPO_PATH

def get_repo_path():
	repo_path = os.getcwd()
	return repo_path

REPO_PATH = get_repo_path()

def read_detection_table(sn_name):
	detections_filename = sn_name + '_detections_table.fits'
	tbdata = fits.getdata(detections_filename, 1)
	return tbdata

def sort_tbdata_by_ra(tbdata):
	sorted_table = tbdata[np.argsort(tbdata['ra'])]
	return Table(sorted_table)

def match_objects(sorted_tbdata, object_threshold = 0.1):
	'''
	object_threshold set in arcseconds
	'''

	#Convert RA and DEC to arcseconds
	ra_to_arcsec = 86400.  #24h*60m*60s
	dec_to_arcsec = 3600.

	sorted_tbdata.add_column(Table.Column(data = sorted_tbdata['ra']*ra_to_arcsec, name = 'ra_arcsec'))
	sorted_tbdata.add_column(Table.Column(data = sorted_tbdata['dec']*dec_to_arcsec, name = 'dec_arcsec'))

	object_id = np.int_(tbdata['ra']*0)-1 #This array will hold the object ID number, initialize with all 0s
	object_number = 0
	for indx, object_ra in enumerate(sorted_tbdata['ra_arcsec']):
		if object_id[indx] == -1: #object has not previously been associated with another object
			object_number += 1
			object_id[indx] = object_number
			#Better to search item by item or using array np.where?
			search_indx = np.where(((sorted_tbdata['ra_arcsec'] - object_ra) > 0) &
						  ((sorted_tbdata['ra_arcsec'] - object_ra) < 2))[0]
			for nearby_obj_indx in search_indx:
				if (object_ra - sorted_tbdata['ra_arcsec'][nearby_obj_indx])**2 + \
					(sorted_tbdata['dec_arcsec'][indx] - sorted_tbdata['dec_arcsec'][nearby_obj_indx])**2 < object_threshold:
					object_id[nearby_obj_indx] = object_number
	sorted_tbdata.add_column(Table.Column(data = object_id, name = 'id'))
	return sorted_tbdata



def build_object_table(sorted_tbdata):
	id_of_most_observed_obj, max_array_size = mode(sorted_tbdata['id'])
	max_array_size = int(max_array_size[0]) #mode returns a single element array, make this a number


	object_table = Table(names = ['id', 'date', 'counts', 'error', 'filter', 'zeropnt', 'snr'],
				dtype = ['i8', '{},f8'.format(max_array_size),'{},f8'.format(max_array_size),
						'{},f8'.format(max_array_size), '{},S3'.format(max_array_size),
						'{},f8'.format(max_array_size),'{},f8'.format(max_array_size)])

	#It looks like set automatically sorts
	for unique_obj in set(sorted_tbdata['id']):
		dates = np.ones((max_array_size,))*999E8  #This number should be larger than any realistic date
		counts = np.ones((max_array_size,))*-999
		error = np.ones((max_array_size,))*-999
		filter = ['N/A']*max_array_size
		zeropt = np.ones((max_array_size,))*-999
		snr = np.ones((max_array_size,))*-999
		object_detection_indx = np.where(sorted_tbdata['id'] == unique_obj)[0]
		for det_num, detection in enumerate(object_detection_indx):
			dates[det_num] = sorted_tbdata['mjd-obs'][detection]
			counts[det_num] = sorted_tbdata['counts'][detection]
			error[det_num] = sorted_tbdata['error'][detection]
			filter[det_num] = sorted_tbdata['filter'][detection]
			zeropt[det_num] = sorted_tbdata['zeropnt'][detection]
			snr[det_num] = sorted_tbdata['snr'][detection]

		filter = np.array(filter)
		#Put in time order
		time_order_indx = np.argsort(dates)
		dates[dates == 999E8] = -999
		#Convert back to lists to build row
		irow = [unique_obj, dates[time_order_indx], counts[time_order_indx],
				error[time_order_indx], filter[time_order_indx],
				zeropt[time_order_indx], snr[time_order_indx]]
		object_table.add_row(irow)
	return object_table

def write_table(sn_name, table, filename = '_object_table.fits'):
	'''
	Remove table file if it exists and write a new one
	table - astropy table object
	filename - name of file to be written
	'''
	filename = sn_name+filename
	try:
		os.remove(filename)
	except OSError:
		pass
	table.write(filename, format = 'fits')


if __name__ == "__main__":
	tbdata = read_detection_table('SN2013fn')
	sorted_table = sort_tbdata_by_ra(tbdata)
	sorted_tbdata = match_objects(sorted_table)
	write_table('SN2013fn', sorted_tbdata, filename = 'object_intermediate_table.fits')
	object_table = build_object_table(sorted_tbdata)
	write_table('SN2013fn', object_table, filename = 'object_table.fits')
