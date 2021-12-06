import os
import ccdproc
from astropy.nddata import CCDData
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
import numpy as np
from astropy.io import ascii
from astroscrappy import detect_cosmics
import time
#	Bottleneck function for faster process
def bn_median(masked_array, axis=None):
    """
    Perform fast median on masked array

    Parameters

    masked_array : `numpy.ma.masked_array`
        Array of which to find the median.

    axis : int, optional
        Axis along which to perform the median. Default is to find the median of
        the flattened array.
    """
    import numpy as np
    import bottleneck as bn
    data = masked_array.filled(fill_value=np.NaN)
    med = bn.nanmedian(data, axis=axis)
    # construct a masked array result, setting the mask from any NaN entries
    return np.ma.array(med, mask=np.isnan(med))
#============================================================
#	Fundamental processing
#------------------------------------------------------------
def master_bias(imlist):
	"""
	"""
	comment = f"""{'-'*60}\n#\tBIAS MASTER FRAME (<--{len(imlist)} frames)\n{'-'*60}"""
	print(comment)
	#	Combine
	combiner = ccdproc.Combiner([CCDData(fits.getdata(inim), unit="adu", meta=fits.getheader(inim)) for inim in imlist])
	mbias = combiner.median_combine(median_func=bn_median)
	#	Header
	mbias.header = fits.getheader(imlist[0])
	mbias.header['NCOMBINE'] = len(imlist)
	for n, inim in enumerate(imlist): mbias.header[f'COMB{n}'] = inim
	#	Save
	dateobs = mbias.header['DATE-OBS'].split('T')[0].replace('-', '')
	filename = f'{dateobs}-zero.fits'
	mbias.header['FILENAME'] = filename
	mbias.write(f'{os.path.dirname(inim)}/{dateobs}-zero.fits', overwrite=True)
	return mbias
#------------------------------------------------------------
def master_dark(imlist, mbias):
	"""
	"""
	comment = f"""{'-'*60}\n#\tDARK MASTER FRAME (<--{len(imlist)} frames)\n{'-'*60}"""
	print(comment)
	#	Combine
	combiner = ccdproc.Combiner([CCDData(fits.getdata(inim), unit="adu", meta=fits.getheader(inim)) for inim in imlist])
	mdark = ccdproc.subtract_bias(combiner.median_combine(median_func=bn_median), mbias)
	#	Header
	mdark.header = fits.getheader(imlist[0])
	mdark.header['NCOMBINE'] = len(imlist)
	mdark.header['SUBBIAS'] = mbias.header['FILENAME']
	for n, inim in enumerate(imlist): mdark.header[f'COMB{n}'] = inim
	#	Save
	exptime = int(mdark.header['EXPTIME'])
	dateobs = mdark.header['DATE-OBS'].split('T')[0].replace('-', '')
	filename = f'{exptime}-{dateobs}-dark.fits'
	mdark.header['FILENAME'] = filename
	mdark.write(f'{os.path.dirname(inim)}/{exptime}-{dateobs}-dark.fits', overwrite=True)
	return mdark
#------------------------------------------------------------
def master_flat(imlist, mbias, mdark, filte=''):
	"""
	"""
	comment = f"""{'-'*60}\n#\t{filte} FLAT MASTER FRAME (<--{len(imlist)} frames)\n{'-'*60}"""
	print(comment)

	def subtract_bias_dark(inim, mbias, mdark):
		flat = CCDData(fits.getdata(inim), unit="adu", meta=fits.getheader(inim))
		bflat = ccdproc.subtract_bias(
			flat,
			mbias,
			)
		dbflat = ccdproc.subtract_dark(
			ccd=bflat,
			master=mdark,
			exposure_time='EXPTIME',
			exposure_unit=u.second,
			scale=True,
		)
		return dbflat
	#	Combine
	combiner = ccdproc.Combiner([subtract_bias_dark(inim, mbias, mdark) for inim in imlist])
	combiner.minmax_clipping()
	def scaling_func(arr): return 1/np.ma.median(arr)
	combiner.scaling = scaling_func
	nmflat = combiner.median_combine(median_func=bn_median)	
	#	Header
	nmflat.header = fits.getheader(imlist[0])
	nmflat.header['NCOMBINE'] = len(imlist)
	nmflat.header['SUBBIAS'] = mbias.header['FILENAME']
	nmflat.header['SUBDARK'] = mdark.header['FILENAME']
	for n, inim in enumerate(imlist): nmflat.header[f'COMB{n}'] = inim
	#	Save
	dateobs = nmflat.header['DATE-OBS'].split('T')[0].replace('-', '')
	filename = f'{dateobs}-n{filte}.fits'
	nmflat.header['FILENAME'] = filename
	nmflat.write(f'{os.path.dirname(inim)}/{filename}', overwrite=True)
	return nmflat
#------------------------------------------------------------
def obj_process(inim, gain, readnoise, mbias, mdark, mflat,):
	'''
	Single image
	'''
	data, hdr = fits.getdata(inim, header=True)
	ccd = CCDData(fits.getdata(inim), unit=u.adu, meta=hdr)
	nccd = ccdproc.ccd_process(
		ccd=ccd,
		# oscan='[201:232,1:100]',
		# trim='[1:200, 1:100]',
		# error=True,
		# gain=gain,
		readnoise=readnoise,
		master_bias=mbias,
		dark_frame=mdark,
		exposure_key='EXPTIME',
		exposure_unit=u.second,
		dark_scale=True,
		master_flat=mflat,
	)
	nccd.header['SUBBIAS'] = mbias.header['FILENAME']
	nccd.header['SUBDARK'] = mdark.header['FILENAME']
	nccd.header['DIVFLAT'] = mflat.header['FILENAME']

	return nccd
#------------------------------------------------------------
def defringe(inim, dfim, outim, dfdat, size=5):
	'''
	inim : image to remove fringe
	dfim : master fringe image
	dfdat : master fringe data
	size : pixel
	'''
	#	Image to process
	data, hdr = fits.getdata(inim, header=True)
	#	Master fringe
	dataf, hdrf	= fits.getdata(dfim, header=True)

	master_fri = fringe_calc(dfim, dfdat, size=size)
	image_fri = fringe_calc(inim, dfdat, size=size)
	fscale = np.median(image_fri/master_fri)
	# print(fscale)
	print(f"{os.path.basename(inim)}/{round(fscale, 3)} (<-- {os.path.basename(dfim)}) ==> df{os.path.basename(inim)}")

	fri_scaled = dataf*fscale
	dfdata = data-fri_scaled
	fits.writeto(outim, dfdata, hdr, overwrite=True)
	return outim
#------------------------------------------------------------
def fringe_calc(dfim, dfdat, size=5):
	data, hdr	= fits.getdata(dfim, header=True)

	dfr_list	= []
	dftbl = ascii.read(dfdat)
	#	Bright position
	xb1, xb2 = dftbl['xb']-size, dftbl['xb']+size
	yb1, yb2 = dftbl['yb']-size, dftbl['yb']+size
	#	Faint position
	xf1, xf2 = dftbl['xf']-size, dftbl['xf']+size
	yf1, yf2 = dftbl['yf']-size, dftbl['yf']+size

	for n in range(len(dftbl)):
		fringe_b= np.median(data[yb1[n]:yb2[n], xb1[n]:xb2[n]])
		fringe_f= np.median(data[yf1[n]:yf2[n], xf1[n]:xf2[n]])		
		dfringe	= fringe_b-fringe_f
		dfr_list.append(dfringe)
	return np.array(dfr_list)
#------------------------------------------------------------
def astrometry(inim, pixscale=None, frac=None, ra=None, dec=None, radius=None, cpulimit=60):
	'''
	ra : hh:mm:ss
	dec : dd:mm:ss
	radius [deg]
	'''
	com = f'solve-field {inim} '
	outim = f'{os.path.dirname(inim)}/a{os.path.basename(inim)}'
	if (pixscale != None) & (type(pixscale) == float):
		if frac == None:
			frac = 0.10 # 10% interval of pixel scale as default
		upscl = pixscale+pixscale*frac
		loscl = pixscale-pixscale*frac
		com = f'{com} --scale-unit arcsecperpix --scale-low {loscl} --scale-high {upscl} '
	if (ra != None) & (dec != None):
		if radius == None:
			radius = 1 # 1 deg as default
		com = f'{com} --ra {ra} --dec {dec} --radius {radius}'
	com = f'{com} --no-plots --new-fits {outim} --overwrite --use-sextractor --cpulimit {cpulimit}'

#------------------------------------------------------------
def cosmic_ray_removal(inim, outim, gain, rdnoise):
	'''
	inim 
	obs = 'LOAO'
	gain = 2.68
	rdnoise = 4.84
	'''

	data, hdr = fits.getdata(inim, header=True)
	param_cr = dict(
					indat=data,
					sigclip=4.5,
					sigfrac=0.3,
					objlim=5.0,
					gain=gain, readnoise=rdnoise, 
					pssl=0.0,
					niter=4,
					sepmed=True,
					cleantype='meanmask',
					fsmode='median',
					psfmodel='gauss',
					psffwhm=hdr['seeing'],
					#	Don't touch
					inmask=None,
					satlevel=65536.0,
					psfsize=7, psfk=None, psfbeta=4.765,
					verbose=False
					)
	time_st = time.time()
	mcrdata, crdata = detect_cosmics(**param_cr)
	fits.writeto(outim, crdata, hdr, overwrite=True)
	time_delta = time.time() - time_st
	indx_cr = np.where(mcrdata == True)
	fits.setval(outim, 'CR', value=len(indx_cr[0]), comment='# of removed cosmic ray by Astroscrappy')
	#------------------------------------------------------------
	print(f"{inim} : Remove {len(indx_cr[0])} cosmic-ray pixels [{round(time_delta, 1)} sec]")

