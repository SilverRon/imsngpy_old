import numpy as np
import matplotlib.pyplot as plt
import os
from astropy import units as u
from astropy.visualization import (ZScaleInterval, MinMaxInterval, SqrtStretch,ImageNormalize)
from astropy.wcs import WCS
from astropy.io import fits

def sexcom(inim, param_input, dualmode=False):
	'''
	
	'''
	param_sex = dict(
		CONF_NAME = 'default.sex',
		#------------------------------
		#	CATALOG
		#------------------------------
		CATALOG_NAME = 'test.cat',
		CATALOG_TYPE = 'ASCII_HEAD',
		PARAMETERS_NAME = 'default.param',
		#------------------------------
		#	EXTRACTION
		#------------------------------
		DETECT_TYPE = 'CCD',
		DETECT_MINAREA = '5',
		DETECT_MAXAREA = '0',
		DETECT_THRESH = '1.5',
		# ANALYSIS_THRESH = 'RELATIVE',
		ANALYSIS_THRESH = '1.5',						
		FILTER = 'Y',
		FILTER_NAME = 'default.conv',
		DEBLEND_NTHRESH = '64',
		DEBLEND_MINCONT = '0.0001',
		CLEAN = 'Y',
		CLEAN_PARAM = '1.0',
		MASK_TYPE = 'CORRECT',
		#------------------------------
		#	PHOTOMETRY
		#------------------------------
		# PHOT_APERTURES = '3',
		PHOT_AUTOPARAMS = '2.5,3.5',
		PHOT_PETROPARAMS = '2.0,3.5',
		SATUR_LEVEL  = '50000.0',
		SATUR_KEY = 'SQTURATE',
		MAG_ZEROPOINT = '0.0',
		MAG_GAMMA = '4.0',
		GAIN = '1.0',
		GAIN_KEY = 'GAIN',   
		PIXEL_SCALE = '1.0',
		#------------------------------
		#	STAR/GALAXY SEPARATION
		#------------------------------
		SEEING_FWHM = '3.0',
		STARNNW_NAME = 'default.nnw',
		#------------------------------
		#	BACKGROUND
		#------------------------------
		BACK_SIZE = '128',
		BACK_FILTERSIZE = '10',
		BACKPHOTO_TYPE = 'LOCAL',
		#------------------------------
		#	CHECK IMAGE
		#------------------------------
		CHECKIMAGE_TYPE = 'NONE',
		CHECKIMAGE_NAME = 'check.fits',
		#==============================
		#	MEMORY & MISCELLANEOUS
		#==============================
		MEMORY_OBJSTACK = '3000',
		MEMORY_PIXSTACK = '300000',
		MEMORY_BUFSIZE = '1024',
		VERBOSE_TYPE = 'NORMAL',
		HEADER_SUFFIX = '.head',
		WRITE_XML = 'N',
		XML_NAME = 'sex.xml'
	)

	for key in param_input.keys():
		param_sex[key] = param_input[key]

	sexcom_normal = f"sex -c {param_sex['CONF_NAME']} {inim} "
	# sexcom_dual = 'sex -c {} {} '.format(param_sex['CONF_NAME'], inim)
	for key in param_sex.keys():
		if key != 'CONF_NAME':
			sexcom_normal += f'-{key} {param_sex[key]} '

	return sexcom_normal
#------------------------------------------------------------
def generate_best_aperture_with_snrcurve(sctbl, apertures, pixscale,):
	'''
	#------------------------------------------------------------
	#	OPTIMIZED APERTURE FROM SNR CURVE
	#------------------------------------------------------------
	apertures*pixscale --> x?
	'''
	indx_col = np.where('SNR'==np.array(sctbl.keys()))
	x=apertures*pixscale.value	#	[pix] --> [arcsec]
	plt.close('all')
	x_optlist = []
	for raw in range(len(sctbl)):
		# print(raw)
		y = np.array(list(sctbl[raw])[indx_col[0].item():])
		y[np.isnan(y)] = 0.0
		indx_peak = np.where(y==np.max(y))
		if len(y)-1 in indx_peak:
			x_opt=None
		else:
			x_opt=x[indx_peak].item()
			#	Plot consumes too much time
			# plt.plot(x, y, color='silver', alpha=0.125)
			# plt.axvline(x=x_opt, ls='-', linewidth=0.5, color='dodgerblue', alpha=0.125)
			# sctbl['APER_OPT'][raw] = x_opt
			x_optlist.append(x_opt)
	# return np.copy(sctbl['APER_OPT'])
	return np.array(x_optlist)/pixscale.value
#------------------------------------------------------------
def draw_snrcurve(scoutpng, title, seeing, aper_opt, dpi=500):
	'''
	
	'''
	# plt.close('all')
	plt.axvline(x=aper_opt, ls='-', linewidth=2.0, color='tomato', alpha=0.5, label=f'OPT.APERTURE : {round(aper_opt, 3)}\"\n(SEEING*{round(aper_opt/seeing.value, 3)})')
	plt.axvline(x=seeing.value, ls='-', linewidth=2.0, color='gold', alpha=0.5, label=f'SEEING : {round(seeing.value, 3)} \"')
	plt.title(title, fontsize=14)
	plt.grid('both', ls='--', color='silver', alpha=0.5)
	plt.xlabel('Aperture Diameter [arcsec]', fontsize=14)
	plt.ylabel('SNR', fontsize=14)
	plt.legend(fontsize=14, framealpha=0.0, loc='upper right')
	# plt.yscale('log')
	plt.savefig(scoutpng, dpi=dpi, overwrite=True)
#------------------------------------------------------------
def select_point_sources(rawtbl, sep, classstarcut=0.9, flagcut=0, magerrcut=0.05, sepcut=0.5*u.deg):
	'''
	'''
	#	Indivisual index
	#	CLASS_STAR
	indx_cs = np.where(
		(rawtbl['CLASS_STAR']>classstarcut)
	)
	#	FLAGS
	indx_fg = np.where(
		(rawtbl['FLAGS']<=flagcut)
	)
	# #	MAG ERROR CUT
	# indx_mg = np.where(
	# 	(rawtbl[errkey]<=magerrcut)
	# )
	#	POSITION
	indx_sp = np.where(
		(sep<sepcut)
	)
	#	All
	indx_sel = np.where(
		(rawtbl['CLASS_STAR']>classstarcut) &
		(rawtbl['FLAGS']<=flagcut) &
		# (rawtbl[errkey]<=magerrcut) &
		(sep<sepcut)
		)
	# return indx_sel, (indx_cs, indx_fg, indx_mg, indx_sp)
	return indx_sel, (indx_cs, indx_fg, indx_sp)
#------------------------------------------------------------
def draw_space_distribution(outpng, c_cent, rawtbl, seltbl, reftbl, dpi=500):
	#	Space distribution
	plt.close('all')
	plt.figure(figsize=(10, 10))
	#	
	plt.plot(rawtbl['ALPHA_J2000'], rawtbl['DELTA_J2000'], ls='none', marker='o', alpha=0.125, label=f'ALL ({len(rawtbl)})')
	plt.plot(seltbl['ALPHA_J2000'], seltbl['DELTA_J2000'], ls='none', marker='+', alpha=0.75, label=f'SELECTED ({len(seltbl)})')
	plt.plot(reftbl['ra'], reftbl['dec'], ls='none', marker='.', mfc='none', alpha=0.125, label=f'REFERENCE ({len(reftbl)})')
	#	Center
	plt.plot(c_cent.ra.value, c_cent.dec.value, ls='none', marker='x', ms=10, c='tomato', mfc='none', alpha=1)
	#	Setting
	plt.legend(fontsize=14, framealpha=1.0, loc='upper right')
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.xlabel('RA [deg]', fontsize=20)
	plt.ylabel('Dec [deg]', fontsize=20)
	xl, xr = plt.xlim()
	plt.xlim([xr, xl])
	plt.grid('both', ls='--', c='silver', alpha=0.5)
	plt.savefig(outpng, overwrite=True, dpi=dpi)
#------------------------------------------------------------
def select4zp(mtbl, filte, sepcut, refmagerrcut, inmagerrcut, refmaglowercut, refmaguppercut):
	'''
	'''
	#	SEP
	indx_sp = np.where(
		(mtbl['sep']<sepcut)
	)
	#	REF MAG ERROR CUT
	indx_er = np.where(
		(mtbl[f'{filte}err']<refmagerrcut)
	)
	#	MAG < REF MAG LOWER
	indx_0 = np.where(
		(mtbl[filte]<refmaglowercut)
	) 
	#	MAG > REF MAG UPPER 
	indx_1 = np.where(
		(mtbl[filte]>refmaguppercut)
	) 
	#	REF MAG LOWER < MAG < REF MAG UPPER 
	indx_2 = np.where(
		(mtbl[filte]>refmaglowercut) &
		(mtbl[filte]<refmaguppercut)
	)
	#	NaN for REF MAG
	indx_n0 = np.where(
		(~np.isnan(mtbl[filte].mask))
	)
	'''#	NaN for MAG
	indx_n1 = np.where(
		(~np.isnan(mtbl[inmagkey].mask))	
	)'''
	'''#	MAG ERR
	indx_iner = np.where(
		(mtbl[errkey]<inmagerrcut)
	)'''
	#	All
	indx_zp = np.where(
		(mtbl['sep']<sepcut) &
		(mtbl[f'{filte}err']<refmagerrcut) &
		(mtbl[filte]>refmaglowercut) &
		(mtbl[filte]>refmaglowercut) &
		(mtbl[filte]<refmaguppercut) &
		(~np.isnan(mtbl[filte].mask))
		# (~np.isnan(mtbl[inmagkey].mask))
	)
	# return indx_zp, (indx_sp, indx_er, indx_0, indx_1, indx_2, indx_n0, indx_1)
	return indx_zp, (indx_sp, indx_er, indx_0, indx_1, indx_2, indx_n0)

# def draw_zpcal(outpng, zptbl, mtbl, filte, czparr, indxes_zp, inzpkey, refmaglowercut, refmaguppercut, dpi=500):
def draw_zpcal(outpng, zptbl, filte, czparr, indxes_zp, inzpkey, refmaglowercut, refmaguppercut, dpi=500):
	'''
	
	'''
	zp = np.median(czparr.data[~czparr.mask])
	zper = np.std(czparr.data[~czparr.mask])

	plt.close()
	plt.figure(figsize=(9,6))
	#
	plt.errorbar(zptbl[filte][~czparr.mask], zptbl[inzpkey][~czparr.mask], yerr=zptbl[f'{filte}err'][~czparr.mask], mfc='none', marker='s', c='dodgerblue', ls='none', alpha=0.25, label=f'{len(np.where(czparr.mask==False)[0])} ({round(1e2*len(np.where(czparr.mask==False)[0])/len(zptbl), 1)}%)')
	plt.plot(zptbl[filte][czparr.mask], zptbl[inzpkey][czparr.mask], marker='x', c='tomato', ls='none', alpha=0.5, label=f'{len(np.where(czparr.mask==True)[0])} ({round(1e2*len(np.where(czparr.mask==True)[0])/len(zptbl), 1)}%)')
	plt.axhspan(zp-zper, zp+zper, alpha=0.125, color='dodgerblue', label=f'{round(zp, 3)}+/-{round(zper, 3)}')
	plt.axhline(y=zp, ls='--', c='b', alpha=0.5)
	#
	# plt.plot(mtbl[filte][indxes_zp[2]], mtbl[inzpkey][indxes_zp[2]], marker='.', c='grey', ls='none', alpha=0.5,)# label=f'{len(indxes_zp[2][0])}')
	# plt.plot(mtbl[filte][indxes_zp[3]], mtbl[inzpkey][indxes_zp[3]], marker='.', c='grey', ls='none', alpha=0.5,)# label=f'{len(indxes_zp[3][0])}')

	plt.axvline(x=refmaglowercut, ls='--', c='k', alpha=0.5)
	plt.axvline(x=refmaguppercut, ls='--', c='k', alpha=0.5)
	#	Setting
	plt.legend(loc='upper center', framealpha=0, fontsize=20)
	plt.xlim([refmaglowercut-1, refmaguppercut+1])
	plt.ylim([zp-0.5, zp+0.5])
	plt.xlabel('Ref. mag [AB]', fontsize=20)
	plt.ylabel('ZP', fontsize=20)
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)
	plt.grid('both', ls='--', color='silver', alpha=0.5)
	plt.tight_layout()
	plt.savefig(outpng, dpi=dpi, overwrite=True)

def draw_zp_space(inim, outpng, c_cent, inzpkey, zptbl, mask, dpi=200):
	# Create interval object
	interval = ZScaleInterval()
	vmin, vmax = interval.get_limits(fits.getdata(inim))

	# Create an ImageNormalize object using a SqrtStretch object
	norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=SqrtStretch())

	plt.close()
	plt.figure(figsize=(10,10))
	# transform=ax.get_transform('world')
	wcs = WCS(inim)
	ax = plt.subplot(projection=wcs)
	ax.imshow(fits.getdata(inim), origin='lower', cmap='gray', norm=norm)

	vmin, vmax = MinMaxInterval().get_limits(zptbl[inzpkey][~mask])

	plt.scatter(zptbl['ALPHA_J2000'][~mask], zptbl['DELTA_J2000'][~mask], c=zptbl[inzpkey][~mask], transform=ax.get_transform('world'), alpha=0.25, vmin=vmin, vmax=vmax)
	plt.colorbar()

	plt.scatter(zptbl['ALPHA_J2000'][mask], zptbl['DELTA_J2000'][mask], c='tomato', marker='x', alpha=0.5, transform=ax.get_transform('world'),)
	plt.plot(c_cent.ra.value, c_cent.dec.value, ls='none', marker='+', ms=10, c='tomato', mfc='none', alpha=1, transform=ax.get_transform('world'))

	# plt.show()
	xl, xr = plt.xlim()
	plt.xlim([xr, xl])
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.xlabel('RA [deg]', fontsize=24)
	plt.ylabel('Dec [deg]', fontsize=24)
	plt.title(f'{os.path.basename(inim)} ({inzpkey})', fontsize=14)
	plt.tight_layout()

	plt.savefig(outpng, dpi=dpi, overwrite=True)

def calc_depth(N, zp, aperture, skysig):
	'''
	N sigma depth (limiting magnitude)
	aper : diameter [pix]
	'''
	r = float(aperture)/2	#	Diameter --> Radius
	braket = N*skysig*np.sqrt(np.pi*(r**2))
	depth = float(zp)-2.5*np.log10(braket)
	return depth