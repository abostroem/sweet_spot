import os
import glob
from astropy.table import Table
from astropy.io import fits
import numpy as np

import pdb

global REPO_PATH


def get_repo_path():
	repo_path = os.getcwd()
	return repo_path

REPO_PATH = get_repo_path()

def get_table_filenames(sn_name):
	'''
	Get list of tables of sources found in each file
	'''
	flist = glob.glob(os.path.join(REPO_PATH,
						'DM',
						'{}'.format(sn_name),
						'{}_A_H*'.format(sn_name),
						'{}_A_H*'.format(sn_name),
						'src.fits'))
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
				source_tbdata['coord'][row_num][0],
				source_tbdata['coord'][row_num][1],
				source_tbdata['centroid_sdss'][row_num][0],
				source_tbdata['centroid_sdss'][row_num][1],
				source_tbdata['flux_psf'][row_num],
				source_tbdata['flux_psf_err'][row_num],
				source_tbdata['flux_psf'][row_num]/source_tbdata['flux_psf_err'][row_num],
				hdr['filter1']]
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
	im#      obj#      ra	dec	x	y	det#      counts sigma_counts SNR	filter zeropt
	'''
	table_flist = get_table_filenames(sn_name)
	for ifile in table_flist:
		with fits.open(ifile) as ofile:
			source_tbdata = ofile[1].data
			hdr = get_image_header(ifile)
			all_objects_table = Table(names = ['filename', 'ID', 'ra', 'dec', 'x', 'y', 'counts', 'error', 'snr', 'filter'],
				dtype = ['S20', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', '<f8', 'S1'])
			for row_num, source_row in enumerate(source_tbdata):
				irow = build_table1_row(ifile.split('/')[-2], source_tbdata, hdr, row_num)
				if irow:
					all_objects_table.add_row(irow)
	return all_objects_table

def write_table(table, filename = 'table1.fits'):
	'''
	Remove table file if it exists and write a new one
	table - astropy table object
	filename - name of file to be written
	'''
	os.remove(filename)
	table.write(filename, format = 'fits')

if __name__ == "__main__":
	all_objects_table = make_table_of_all_objects('SN2013fn')
	write_table(all_objects_table)
