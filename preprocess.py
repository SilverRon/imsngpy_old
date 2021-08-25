import os
import ccdproc
from astropy.nddata import CCDData
import matplotlib.pyplot as plt
from ccdproc import ImageFileCollection

def master_zero(
	images,
	fig=False
	):
	"""
	"""
	comment     = '-'*60+'\n' \
				+ 'MAKING MASTER ZERO\n' \
				+ '-'*60+'\n'
	print(comment)
	path_data = images.location
	#	HEADER FOR MASTER FRAME
	ccddata_lst   = []
	# for hdu, fname in images.hdus(imagetyp='Bias', return_fname=True):
	for hdu, fname in images.hdus(imagetyp='BIAS', return_fname=True):
		
		meta = hdu.header
		ccddata_lst.append(ccdproc.CCDData(data=hdu.data, meta=meta, unit="adu"))

	#	HEADER FOR MASTER FRAME
	n = 0
	for hdu, fname in images.hdus(imagetyp='Bias', return_fname=True):	
		n += 1
		newmeta = meta
		newmeta['FILENAME'] = 'zero.fits'
		newmeta['COMB{}'.format(n)] = fname
	print('{} ZERO IMAGES WILL BE COMBINED.'.format(len(ccddata_lst)))
	zeros = ccdproc.Combiner(ccddata_lst)
	mzero = zeros.median_combine()
	mzero.header  = newmeta

	#	SAVE MASTER FRAME
	if '{}/zero.fits'.format(path_data) in glob.glob(path_data+'/*'):
		os.system('rm {}/zero.fits'.format(path_data))
	mzero.write('{}/zero.fits'.format(path_data))
	if fig != False:
		imstats = lambda dat: (dat.min(), dat.max(), dat.mean(), dat.std())
		zero_min, zero_max, zero_mean, zero_std = imstats(np.asarray(mzero))
		plt.figure(figsize=(15, 15))
		plt.imshow(mzero, vmax=zero_mean + 4*zero_std, vmin=zero_mean - 4*zero_std)
		plt.savefig(path_data+'/zero.png')
	return mzero