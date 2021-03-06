import os
import glob
from astropy.table import Table
from astropy.io import fits
import numpy as np

import pdb

global REPO_PATH

#Detections table
#Objects table
def get_repo_path():
	repo_path = os.getcwd()
	return repo_path

REPO_PATH = get_repo_path()

def get_table_filenames(sn_name):
	'''
	Get list of tables of sources found in each file
	'''
	flist_H = glob.glob(os.path.join(REPO_PATH,
						'DM',
						'{}'.format(sn_name),
						'{}_A_H*'.format(sn_name),
						'{}_A_H*'.format(sn_name),
						'src.fits'))

	flist_J = glob.glob(os.path.join(REPO_PATH,
						'DM',
						'{}'.format(sn_name),
						'{}_A_J*'.format(sn_name),
						'{}_A_J*'.format(sn_name),
						'src.fits'))
	flist = flist_H + flist_J

	return flist

def build_table1_row(filename, source_tbdata, hdr, row_num):
	'''
	Create a row to be put into table 1 from info in the original DM output table
	'''
	x, y = source_tbdata['centroid_sdss'][row_num]
	if not (np.isfinite(x) and np.isfinite(y)):
		print source_tbdata['centroid_sdss'][row_num]
		return None
	else:
		irow = [filename,
				source_tbdata['id'][row_num],
				hdr['mjd-obs'],
				source_tbdata['coord'][row_num][0],
				source_tbdata['coord'][row_num][1],
				source_tbdata['centroid_sdss'][row_num][0],
				source_tbdata['centroid_sdss'][row_num][1],
				source_tbdata['flux_psf'][row_num],
				source_tbdata['flux_psf_err'][row_num],
				source_tbdata['flux_psf'][row_num]/source_tbdata['flux_psf_err'][row_num],
				hdr['filter1'], hdr['zeropnt']]
		return irow

def get_image_header(ifile):
	'''
	Get header from image file (based on table file name)
	'''
	image_name = ifile.split('/')[-2]+'.fits'
	hdr = fits.getheader(os.path.join(REPO_PATH, 'images', image_name))
	return hdr

def make_table_of_all_objects(sn_name):
	'''
	Combine information from each image into a table of the form
	im#      obj#   mjd_obs   ra	dec	x	y	det#      counts sigma_counts SNR	filter zeropt
	ID = linear numbers
	det# + im# =
	'''
	table_flist = get_table_filenames(sn_name)
	all_objects_table = Table(names = ['filename', 'det', 'mjd-obs','ra', 'dec', 'x', 'y', 'counts', 'error', 'snr', 'filter', 'zeropnt'],
	dtype = ['S20', '<i8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', 'S1', '<f8'])
	for ifile in table_flist:
		with fits.open(ifile) as ofile:
			source_tbdata = ofile[1].data
			hdr = get_image_header(ifile)
			for row_num, source_row in enumerate(source_tbdata):
				irow = build_table1_row(ifile.split('/')[-2], source_tbdata, hdr, row_num)
				if irow:
					all_objects_table.add_row(irow)
	return all_objects_table

def write_table(sn_name, table, filename = '_detections_table.fits'):
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
	all_objects_table = make_table_of_all_objects('SN2013fn')
	write_table('SN2013fn', all_objects_table)
