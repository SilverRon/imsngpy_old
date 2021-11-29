import os
import ccdproc
from astropy.nddata import CCDData
import matplotlib.pyplot as plt
from astropy.io import fits

def master_bias(imlist):
	"""
	"""
	comment = f"""{'-'*60}\n#\tBIAS MASTER FRAME (<--{len(imlist)} frames)\n{'-'*60}"""
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
		mbias.header[f'COMB{n}'] = inim
	dateobs = mbias.header['DATE-OBS'].split('T')[0].replace('-', '')
	filename = f'{dateobs}-zero.fits'
	# filename = f'{os.path.dirname(inim)}/{dateobs}-zero.fits'
	mbias.header['FILENAME'] = filename
	mbias.write(f'{os.path.dirname(inim)}/{dateobs}-zero.fits', overwrite=True)
	return mbias

def master_dark(imlist, mbias):
	"""
	"""
	comment = f"""{'-'*60}\n#\tDARK MASTER FRAME (<--{len(imlist)} frames)\n{'-'*60}"""
	print(comment)

	# st = time.time()
	combiner = ccdproc.Combiner([CCDData(fits.getdata(inim), unit="adu", meta=fits.getheader(inim)) for inim in imlist])
	mdark = ccdproc.subtract_bias(combiner.median_combine(), mbias)
	# print(time.time() - st)
	mdark.header['NCOMBINE'] = len(imlist)
	mdark.header['SUBBIAS'] = mbias.header['FILENAME']
	for n, inim in enumerate(imlist):
		if n==0:
			mdark.header = fits.getheader(inim)
		else:
			pass
		mdark.header[f'COMB{n}'] = inim
	exptime = int(mdark.header['EXPTIME'])
	dateobs = mdark.header['DATE-OBS'].split('T')[0].replace('-', '')
	filename = f'{exptime}-{dateobs}-dark.fits'
	# filename = f'{os.path.dirname(inim)}/{exptime}-{dateobs}-dark.fits'
	mdark.header['FILENAME'] = filename
	mdark.write(f'{os.path.dirname(inim)}/{exptime}-{dateobs}-dark.fits', overwrite=True)
	return mdark

def master_flat(imlist, mbias, mdark):
	"""
	"""
	comment = f"""{'-'*60}\n#\tFLAT MASTER FRAME (<--{len(imlist)} frames)\n{'-'*60}"""
	print(comment)

	# st = time.time()
	combiner = ccdproc.Combiner([CCDData(fits.getdata(inim), unit="adu", meta=fits.getheader(inim)) for inim in imlist])
	mdark = ccdproc.subtract_bias(combiner.median_combine(), mbias)
	# print(time.time() - st)
	mdark.header['NCOMBINE'] = len(imlist)
	mdark.header['SUBBIAS'] = mbias.header['FILENAME']
	for n, inim in enumerate(imlist):
		if n==0:
			mdark.header = fits.getheader(inim)
		else:
			pass
		mdark.header[f'COMB{n}'] = inim
	exptime = int(mdark.header['EXPTIME'])
	dateobs = mdark.header['DATE-OBS'].split('T')[0].replace('-', '')
	filename = f'{exptime}-{dateobs}-dark.fits'
	# filename = f'{os.path.dirname(inim)}/{exptime}-{dateobs}-dark.fits'
	mdark.header['FILENAME'] = filename
	mdark.write(f'{os.path.dirname(inim)}/{exptime}-{dateobs}-dark.fits', overwrite=True)
	return mdark
