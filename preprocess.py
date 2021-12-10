import os
import ccdproc
from astropy.nddata import CCDData
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
from astropy.io import fits
from astropy import units as u
import numpy as np
from astropy.io import ascii
from astroscrappy import detect_cosmics
import time
from astropy.coordinates import SkyCoord
from astropy.table import Table
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

	hdr['HISTORY'] = 'Fringe pattern correction'
	hdr['FRINGEF'] = (round(fscale, 3), 'Flux scaling for defringe [image/master]')

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
def cosmic_ray_removal(inim, outim, gain, rdnoise, seeing=3*u.arcsec, cleantype='medmask'):
	'''
	inim 
	obs = 'LOAO'
	gain = 2.68
	rdnoise = 4.84
	'''
	if gain.unit == u.electron/u.adu:
		gain = gain.value
	if rdnoise.unit == u.electron:
		rdnoise = rdnoise.value
	if seeing.unit == u.arcsec:
		seeing = seeing.value
		
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
					# cleantype='meanmask',
					# cleantype='idw',
					# cleantype='median',
					# cleantype='medmask',
					cleantype=cleantype,
					fsmode='median',
					psfmodel='gauss',
					psffwhm=hdr['seeing'],
					#	Don't touch
					# inbkg=None,
					# inval=None,
					inmask=None,
					satlevel=65536.0,
					psfsize=7, psfk=None, psfbeta=4.765,
					verbose=False
					)
	time_st = time.time()
	mcrdata, crdata = detect_cosmics(**param_cr)
	ncr = len(crdata[mcrdata])
	hdr['HISTORY'] = 'Cosmic-ray and bad pixels correction with astroscrappy'
	hdr['CRNUMB'] = (ncr, '# of removed cosmic ray by Astroscrappy')
	hdr['CRTYPE'] = (cleantype, 'Clean type for cosmic-ray removal')

	fits.writeto(outim, crdata, hdr, overwrite=True)
	#	Mask array
	moutim = ''.join([os.path.splitext(outim)[0], '.mask.fits'])
	fits.writeto(moutim, mcrdata+0, hdr, overwrite=True)
	time_delta = time.time() - time_st
	#------------------------------------------------------------
	print(f"Remove {ncr} cosmic-ray pixels [{round(time_delta, 1)} sec]")
	print(f'\t{os.path.basename(inim)} --> {os.path.basename(outim)}')
	print(f'\tMask : {os.path.basename(moutim)}')

#------------------------------------------------------------
def astrometry(inim, outim, pixscale=None, frac=None, ra=None, dec=None, radius=None, cpulimit=60):
	'''
	ra : hh:mm:ss
	dec : dd:mm:ss
	radius [deg]
	'''
	if pixscale.unit == u.arcsec/u.pixel:
		pixscale = pixscale.value
	if (radius.unit == u.deg) | (radius.unit == u.arcmin) | (radius.unit == u.arcsec):
		radius = radius.to(u.deg).value
	com = f'solve-field {inim} '
	if (pixscale != None) & (type(pixscale) == float):
		if frac == None:
			frac = 0.10 # 10% interval of pixel scale as default
		upscl = pixscale*(1+frac)
		loscl = pixscale*(1-frac)
		com = f'{com} --scale-unit arcsecperpix --scale-low {loscl} --scale-high {upscl} '
	if (ra != None) & (dec != None):
		if radius == None:
			radius = 1 # 1 deg as default
		com = f'{com} --ra {ra} --dec {dec} --radius {radius}'
	com = f'{com} --no-plots --new-fits {outim} --overwrite --use-sextractor --cpulimit {cpulimit}'
	print(f'{os.path.basename(inim)} --> {os.path.basename(outim)}')
	os.system(com)
#------------------------------------------------------------
def astrometry_analysis(inim, incor, outpng, outdat):
	'''
	inim = 'acrdffdzobj.NGC4108.20210622.0158.fits'
	incor = 'crdffdzobj.NGC4108.20210622.0158.corr'
	outpng = 'crdffdzobj.NGC4108.20210622.0158.astrm.png'
	outdat = 'crdffdzobj.NGC4108.20210622.0158.astrm.dat'
	'''
	data, hdr = fits.getdata(inim, header=True)
	cortbl = Table(fits.getdata(incor))

	#	INDEX astrometry.net (ref. catalog)
	c_indx = SkyCoord(cortbl['index_ra'], cortbl['index_dec'], unit='deg')
	#	Measured position
	c_field = SkyCoord(cortbl['field_ra'], cortbl['field_dec'], unit='deg')

	#	Results
	ra_offset, dec_offset = c_field.spherical_offsets_to(c_indx)
	sep = c_field.separation(c_indx)
	ra_rms = np.sqrt(np.mean(ra_offset.arcsec**2))
	dec_rms = np.sqrt(np.mean(dec_offset.arcsec**2))
	sep_rms = np.sqrt(np.mean(sep.arcsec**2))

	fits.setval(inim, keyword='N_ASTRM', value=len(sep), comment='# of sources for astrometry')
	fits.setval(inim, keyword='SEPRMS', value=round(sep_rms, 3), comment='Astrometry sepearation rms [arcsec]')
	fits.setval(inim, keyword='RARMS', value=round(ra_rms, 3), comment='Astrometry RA offset rms [arcsec]')
	fits.setval(inim, keyword='DERMS', value=round(dec_rms, 3), comment='Astrometry Dec offset rms [arcsec]')

	#	Plot
	#https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/scatter_hist.html
	plt.close('all')
	plt.figure(figsize=(8, 8))
	# fig, ax = plt.subplots(figsize=(5, 5))
	#	Size of Axes
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.65
	# spacing = 0.005
	spacing = 0.

	rect_scatter = [left, bottom, width, height]
	rect_histx = [left, bottom + height + spacing, width, 0.2]
	rect_histy = [left + width + spacing, bottom, 0.2, height]

	#	Main axe
	ax = plt.axes(rect_scatter)
	ax.tick_params(direction='in', top=True, right=True)
	ax.tick_params(axis='both', labelsize=14)
	#	Right and Left axes (ax1, ax2 for eact)
	ax1 = plt.axes(rect_histx)
	ax1.tick_params(direction='in', labelbottom=False)
	ax2 = plt.axes(rect_histy)
	ax2.tick_params(direction='in', labelleft=False)

	#	Main plot
	ax.plot(ra_offset.arcsec, dec_offset.arcsec, c='grey', alpha=0.75, ls='none', marker='o', mfc='none', label=f'{len(sep)} sources')
	ax.errorbar(0, 0, xerr=ra_rms, yerr=dec_rms, capsize=5, c='tomato', label=f"RA offset:{round(ra_rms, 2)}\nDec offset:{round(dec_rms, 2)}\nsep:{round(sep_rms, 2)}")
	if 'SEEING' in hdr.keys():
		ax.add_patch(plt.Circle((0, 0), hdr['SEEING']/2, color='silver', alpha=0.25, label=f"SEEING:{round(hdr['SEEING'], 2)}"))
	ax.set_xlim([-3.5, +3.5])
	ax.set_ylim([-3.5, +3.5])
	ax.set_xlabel('RA offset [arcsec]', fontsize=20)
	ax.set_ylabel('Dec offset [arcsec]', fontsize=20)
	# ax.legend(fontsize=14, loc='lower right', framealpha=0.0, edgecolor='k')
	ax.legend(fontsize=14, loc='lower right', facecolor='none', edgecolor='k')
	ax.grid('both', c='silver', ls='--', alpha=0.5)

	#	x, y ranges
	binwidth = 0.25
	lim = np.ceil(np.abs([ra_offset.arcsec, dec_offset.arcsec]).max() / binwidth) * binwidth
	ax.set_xlim((-lim, lim))
	ax.set_ylim((-lim, lim))

	#	Sub plots
	bins = np.arange(-lim, lim + binwidth, binwidth)
	ax1.hist(ra_offset.arcsec, bins=bins, histtype='step', color='k')
	ax1.axvline(x=np.median(ra_offset.arcsec), ls='--', color='k', alpha=0.5, label=f'Median:{round(np.median(ra_offset.arcsec), 2)}')
	ax1.set_xlim(ax.get_xlim())
	ax1.legend(fontsize=12, loc='upper right', framealpha=0.0)
	ax1.grid('both', c='silver', ls='--', alpha=0.5)
	ax1.set_title(os.path.basename(inim), fontsize=14)

	ax2.hist(dec_offset.arcsec, bins=bins, histtype='step', orientation='horizontal', color='k')
	ax2.axhline(y=np.median(dec_offset.arcsec), ls='--', color='k', alpha=0.5, label=f'{round(np.median(dec_offset.arcsec), 2)}')
	ax2.set_ylim(ax.get_ylim())
	ax2.legend(fontsize=12, loc='upper right', framealpha=0.0)
	ax2.grid('both', c='silver', ls='--', alpha=0.5)

	plt.savefig(outpng, overwrite=True)