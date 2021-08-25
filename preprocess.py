import os
import ccdproc
from astropy.nddata import CCDData
import matplotlib.pyplot as plt
from astropy.io import fits

def master_bias(imlist):
	"""
	"""
	comment = f"""{'-'*60}
	#\tBIAS MASTER FRAME (<--{len(imlist)} frames)
	{'-'*60}"""
	print(comment)

	# st = time.time()
	combiner = ccdproc.Combiner([CCDData(fits.getdata(inim), unit="adu", meta=fits.getheader(inim)) for inim in imlist])
	mbias = combiner.median_combine()
	# print(time.time() - st)
	mbias.header['NCOMBINE'] = len(imlist)
	for n, inim in enumerate(imlist):
		if n==0:
			mbias.header = fits.getheader(inim)
		else:
			pass
		mbias.header['COMB{}'.format(n)] = inim
	dateobs = mbias.header['DATE-OBS'].split('T')[0].replace('-', '')
	filename = f'{os.path.dirname(inim)}/{dateobs}-zero.fits'
	mbias.header['FILENAME'] = filename
	return mbias

def master_dark(imlist, mbias, exptime):
	"""
	"""
	comment = f"""{'-'*60}
	#\tDARK MASTER FRAME (<--{len(imlist)} frames)
	{'-'*60}"""
	print(comment)

	# st = time.time()
	combiner = ccdproc.Combiner([CCDData(fits.getdata(inim), unit="adu", meta=fits.getheader(inim)) for inim in imlist])
	mdark = combiner.median_combine()
	# print(time.time() - st)
	mdark.header['NCOMBINE'] = len(imlist)
	mdark.header['SUBBIAS'] = mbias.header['FILENAME']
	for n, inim in enumerate(imlist):
		if n==0:
			mdark.header = fits.getheader(inim)
		else:
			pass
		mdark.header['COMB{}'.format(n)] = inim
	dateobs = mdark.header['DATE-OBS'].split('T')[0].replace('-', '')
	filename = f'{os.path.dirname(inim)}/{dateobs}-dark.fits'
	mdark.header['FILENAME'] = filename
	return mdark
