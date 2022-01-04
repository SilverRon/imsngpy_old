#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#============================================================
#%%
#	Library
#------------------------------------------------------------
import time
import os
import sys
from datetime import date
import numpy as np
import subprocess
#	IMSNGpy modules
sys.path.append('/home/paek/imsngpy')
from tableutil import *
from preprocess import *
from misc import *
from phot import *
from util import *
from query import *
#	Astropy
from astropy.io import ascii
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.coordinates import SkyCoord
from astropy.table import hstack
from astropy.wcs import WCS
from astropy import units as u
#	Snapshot
from astropy.nddata import Cutout2D
from matplotlib.patches import Circle, PathPatch
from astropy.visualization import SqrtStretch, LinearStretch, LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import ZScaleInterval, MinMaxInterval
#	Bottle Neck
import bottleneck as bn
#	Multiprocess tools
from itertools import repeat
import multiprocessing
st_ = time.time()
#============================================================
#%%
#	Path
#------------------------------------------------------------
path_base = '/home/paek/imsngpy'
path_table = f'{path_base}/table'
path_config = f'{path_base}/config'
path_gphot = f'{path_config}/gphot.config'
path_search = f'{path_config}/search.config'
#	Table
ccdtbl = ascii.read(f'{path_table}/ccd.tsv') 
#	gphot Configuration
if os.path.exists(path_gphot) == True:
	gphot_dict = file2dict(path_gphot)
else:
	#	Default gpphot.config
	print('[NOTICE] There is no gphot.config. Use default configuration.')
	gphot_dict = {
		'imkey': './Calib*com.fits',
		'photfraction': '0.75',
		'refcatname': 'PS1',
		'refqueryradius': '1',
		'refmaglower': '14',
		'refmagupper': '18',
		'refmagerupper': '0.05',
		'inmagerupper': '0.05',
		'flagcut': '0',
		'sepfrac': '0.5',
		'DETECT_MINAREA': '5',
		'DETECT_THRESH': '3.0',
		'DEBLEND_NTHRESH': '64',
		'DEBLEND_MINCONT': '0.0001',
		'BACK_SIZE': '64',
		'BACK_FILTERSIZE': '3',
		'BACKPHOTO_TYPE': 'LOCAL',
		'check': 'False'
		}
#	gphot Configuration
if os.path.exists(path_search) == True:
	search_dict = file2dict(path_search)
else:
	#	Default search.config
	print('[NOTICE] There is no search.config. Use default configuration.')
	search_dict = {'frac_sb': '5.0',
		'frac_inv': '1.0',
		'frac': '0.9',
		'inmagkey': "'mag_aper'",
		'inul5key': "'ul5'",
		'cutline_rellip': '5',
		'cutline_flags': '4.0',
		'seeing_up': '3.0',
		'seeing_lo': '0.5',
		'cutline_bkg': '50',
		'key0': "'mag_aper_1'",
		'key1': "'mag_aper_3'",
		'n_delm': '10',
		'frac_dim': '0.3',
		'cutline_null': '0'
		}
#------------------------------------------------------------
#	Configuration files
#------------------------------------------------------------
#	gpphot
prefix_gp = 'gpphot'
path_param = f'{path_config}/{prefix_gp}.param'
path_conv = f'{path_config}/{prefix_gp}.conv'
path_nnw = f'{path_config}/{prefix_gp}.nnw'
path_conf = f'{path_config}/{prefix_gp}.sex'
#	invert
prefix_iv = 'invert'
path_param_iv = f'{path_config}/{prefix_iv}.param'
path_conv_iv = f'{path_config}/{prefix_iv}.conv'
path_nnw_iv = f'{path_config}/{prefix_iv}.nnw'
path_conf_iv = f'{path_config}/{prefix_iv}.sex'
#============================================================
#	Input
#------------------------------------------------------------
tstablename = sys.argv[1]
# tstablename = '/data3/paek/factory/loao/test_fast/transient_search.ecsv'
print(f"{'='*60}\n")
print('INPUT')
print('-'*60)
print(f'Transient Search Table  : {tstablename}')
tstbl = ascii.read(tstablename)
try:
	ncore = int(sys.argv[2])
except:
	ncore = 4
print(f'NCORE  : {ncore}')
print(f"{'='*60}\n")
#------------------------------------------------------------
#	Function
#------------------------------------------------------------
def routine_phot(inim,):
	print(f'PHOTOMETRY START for {os.path.basename(inim)}')
	#
	global path_param
	global path_conv
	global path_nnw
	global path_conf
	#------------------------------------------------------------
	#	Image info
	#------------------------------------------------------------
	hdr = fits.getheader(inim)
	#	IMAGE PHYSICAL CENTER
	a = hdr['naxis1']/2.
	b = hdr['naxis2']/2.
	#	Small squre based on frac
	frac = float(gphot_dict['photfraction'])
	a_ = a*np.sqrt(frac)
	b_ = b*np.sqrt(frac)
	obj = hdr['OBJECT']
	filte = hdr['FILTER']
	#------------------------------------------------------------
	#	CCD INFO
	#------------------------------------------------------------
	if ('OBSERVAT' in hdr.keys()) & ('CCDNAME' in hdr.keys()):
		obs = hdr['OBSERVAT']
		ccd = hdr['CCDNAME']
	else:
		#	observatory from filename
		obs = os.path.basename(inim).split('-')[1]
		if '_' not in obs:
			ccd = ''
		else:
			obs = obs.split('_')[0]
			ccd = obs.split('_')[1]
	indx_ccd = np.where(
		(ccdtbl['obs']==obs) &
		(ccdtbl['ccd']==ccd)
	)
	
	# print(f"""{'-'*60}\n#\tCCD INFO\n{'-'*60}""")
	gain = ccdtbl['gain'][indx_ccd][0]*(u.electron/u.adu)
	rdnoise = ccdtbl['readnoise'][indx_ccd][0]*(u.electron)
	pixscale = ccdtbl['pixelscale'][indx_ccd][0]*(u.arcsec/u.pixel)
	fov = ccdtbl['foveff'][indx_ccd][0]*(u.arcmin)
	# print(f"""GAIN : {gain}\nREAD NOISE : {rdnoise}\nPIXEL SCALE : {pixscale}\nEffective FoV : {fov}""")
	#------------------------------------------------------------
	#	Single
	if ('SEEING' in hdr.keys()) & ('PEEING' in hdr.keys()):
		seeing = hdr['SEEING']*u.arcsec
		peeing = hdr['PEEING']*u.pix
	else:
		print('No seeing information on the header. Calculate ')
		#	If no seeing info, get it with simple conf
		prefix_sp = 'simple'
		path_conf = f'{path_config}/{prefix_sp}.sex'
		path_param = f'{path_config}/{prefix_sp}.param'
		path_nnw = f'{path_config}/{prefix_sp}.nnw'
		path_conv = f'{path_config}/{prefix_sp}.conv'
		seeing, peeing = get_seeing(
			inim,
			gain, 
			pixscale, 
			fov, 
			path_conf, 
			path_param, 
			path_conv, 
			path_nnw, 
			seeing_assume=3*u.arcsec, 
			frac=0.68, 
			n_min_src=5
			)

	#------------------------------------------------------------
	#%%
	#	OPTIMIZED APERTURE FROM SNR CURVE
	#------------------------------------------------------------
	aper_opt = hdr['APER']*pixscale.value
	#------------------------------------------------------------
	#%%
	#	APERTURE DICTIONARY
	#------------------------------------------------------------
	aper_dict = dict(
		MAG_AUTO = dict(
			size=0.0,
			errkey='MAGERR_AUTO',
			suffix='AUTO',
			comment='',
		),
		MAG_APER = dict(
			size=aper_opt/pixscale.value,
			errkey='MAGERR_APER',
			suffix='APER',
			comment='Best aperture diameter based on SNR curve [pix]',
		),
		MAG_APER_1 = dict(
			size=2*0.6731*peeing.value,
			errkey='MAGERR_APER_1',
			suffix='APER_1',
			comment='2*0.6731*seeing (gaussian profile) [pix]',
		),
		MAG_APER_2 = dict(
			size=2*peeing.value,
			errkey='MAGERR_APER_2',
			suffix='APER_2',
			comment='2*seeing aperture diameter [pix]',
		),
		MAG_APER_3 = dict(
			size=3*peeing.value,
			errkey='MAGERR_APER_3',
			suffix='APER_3',
			comment='3*seeing aperture diameter [pix]',
		),
		MAG_APER_4 = dict(
			size=(3*u.arcsec/pixscale).value,
			errkey='MAGERR_APER_4',
			suffix='APER_4',
			comment='Fixed 3\" aperture diameter [pix]',
		),
		MAG_APER_5 = dict(
			size=(5*u.arcsec/pixscale).value,
			errkey='MAGERR_APER_5',
			suffix='APER_5',
			comment='Fixed 5\" aperture diameter [pix]',
		),
	)
	#------------------------------------------------------------
	#%%
	#	SET APERTURES
	#	MAG_AUTO, MAG_APER, MAG_APER_1, MAG_APER_2, MAG_APER_3, MAG_APER_4, MAG_APER_5
	apertures_phot = []
	for magkey in list(aper_dict.keys()):
		if magkey != 'MAG_AUTO':
			apertures_phot.append(aper_dict[magkey]['size'])
	apertures_phot = np.array(apertures_phot)
	apertures_input = ','.join([str(d) for d in apertures_phot])
	#	Run source extractor for photometry
	outcat = f'{os.path.splitext(inim)[0]}.cat'
	param_insex = dict(
		#------------------------------
		#	CATALOG
		#------------------------------
		CATALOG_NAME = outcat,
		#------------------------------
		#	CONFIG FILES
		#------------------------------
		CONF_NAME = path_conf,
		PARAMETERS_NAME = path_param,
		FILTER_NAME = path_conv,    
		STARNNW_NAME = path_nnw,
		#------------------------------
		#	PHOTOMETRY
		#------------------------------
		PHOT_APERTURES = apertures_input,
		GAIN = str(gain.value),
		PIXEL_SCALE = str(pixscale.value),
		#------------------------------
		#	STAR/GALAXY SEPARATION
		#------------------------------
		SEEING_FWHM = str(seeing.value),
		#------------------------------
		#	EXTRACTION
		#------------------------------
		DETECT_MINAREA = gphot_dict['DETECT_MINAREA'],
		DETECT_THRESH = gphot_dict['DETECT_THRESH'],
		DEBLEND_NTHRESH = gphot_dict['DEBLEND_NTHRESH'],
		DEBLEND_MINCONT = gphot_dict['DEBLEND_MINCONT'],
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
	)
	# os.system(sexcom(inim, param_insex))
	# rawtbl = ascii.read(outcat)
	
	sexout = subprocess.getoutput(sexcom(inim, param_insex))
	line = [s for s in sexout.split('\n') if 'RMS' in s]
	skyval, skysig = float(line[0].split('Background:')[1].split('RMS:')[0]), float(line[0].split('RMS:')[1].split('/')[0])
	rawtbl = ascii.read(outcat)
	
	#------------------------------------------------------------
	#	GENERAL INFO. REGARDLESS OF APERTURE
	#------------------------------------------------------------
	#	VALUES
	rsky = 1e2*skysig/skyval
	#	PUT HEADER ON THE IMAGE
	fits.setval(inim, 'AUTHOR', value='Gregory S.H. Paek', comment='PHOTOMETRY AUTHOR')
	fits.setval(inim, 'PHTDATE', value=date.today().isoformat(), comment='PHOTOMETRY DATE [KOR]')
	fits.setval(inim, 'SKYSIG', value=round(skysig, 3), comment='STD of BACKGROUND [pixel]')
	fits.setval(inim, 'SKYVAL', value=round(skyval, 3), comment='MEDIAN of BACKGROUND [pixel]')
	fits.setval(inim, 'RSKY', value=round(rsky, 3), comment='RATIO of SKYSIG/SKYVAL [%]')

	#============================================================
	#%%
	#	ZP CALCULATION
	#------------------------------------------------------------
	# n = 0
	# inmagkey = list(aper_dict.keys())[n]
	for n, inmagkey in enumerate(list(aper_dict.keys())):
		inerrkey = aper_dict[inmagkey]['errkey']
		inzpkey = f"ZP_{aper_dict[inmagkey]['suffix']}"
		inul3key = f"UL3{aper_dict[inmagkey]['suffix'].replace('APER', '')}"
		inul5key = f"UL5{aper_dict[inmagkey]['suffix'].replace('APER', '')}"
		hdrzpkey = inzpkey.replace('_APER', '')
		hdrzperkey = hdrzpkey.replace('ZP', 'ZPER')
		hdrnstdkey = hdrzpkey.replace('ZP', 'NZP')

		#	MERGED TABLE
		rawtbl[inmagkey.lower()] = rawtbl[inmagkey]+hdr[hdrzpkey]
		rawtbl[aper_dict[inmagkey]['errkey'].lower()] = sqsum(rawtbl[inerrkey], hdr[hdrzperkey])

	rawtbl.write(outcat, format='ascii.ecsv', overwrite=True)
	delt = time.time()-st_
	print(f"DONE ({round(delt, 1)} sec)\n")
#------------------------------------------------------------
def mp_phot(imlist, ncore=4):
	if __name__ == '__main__':
		#	Fixed the number of cores (=4)
		with multiprocessing.Pool(processes=ncore) as pool:
			results = pool.starmap(
				routine_phot,
				zip(
					imlist,
					)
				)
#------------------------------------------------------------
def phot(imlist):
	for i, inim in enumerate(imlist):
		print(f"[{i+1}/{len(imlist)}] {os.path.basename(inim)}")
		routine_phot(inim)
#------------------------------------------------------------
def invert_data(inim, outim):
	data, hdr = fits.getdata(inim, header=True)
	fits.writeto(outim, data*-1, header=hdr, overwrite=True)
#------------------------------------------------------------
def routine_search(sciim, refim, subim, search_dict):
	import numpy as np
	print('='*60)
	print(f'Transient Search')
	print('-'*60)
	print(f"SCI :{os.path.basename(sciim)}")
	print(f"REF :{os.path.basename(refim)}")
	print(f"SUB :{os.path.basename(subim)}")
	print('-'*60)
	#	Sci. image
	hdr = fits.getheader(sciim)
	c_cent = SkyCoord(hdr['CRVAL1'], hdr['CRVAL2'], unit=u.deg)
	dateobs = hdr['DATE-OBS']
	epoch = Time(dateobs, format='isot')
	seeing = hdr['SEEING']*u.arcsec
	ccdinfo = get_ccdinfo(sciim, ccdtbl)
	#	Sub. image 
	data_sub, hdr_sub = fits.getdata(subim, header=True)
	peeing = hdr_sub['PEEING']
	skysig = hdr_sub['SKYSIG']
	#------------------------------------------------------------
	#	INVERT
	#------------------------------------------------------------
	isciim = f"{os.path.dirname(sciim)}/inv{os.path.basename(sciim)}"
	isubim = f"{os.path.dirname(subim)}/inv{os.path.basename(subim)}"
	irefim = f"{os.path.dirname(refim)}/inv{os.path.basename(refim)}"

	invert_data(sciim, isciim)
	invert_data(refim, irefim)
	invert_data(subim, isubim)
	#------------------------------------------------------------
	scicat = f'{os.path.splitext(sciim)[0]}.cat'
	#------------------------------------------------------------
	subcat = f'{os.path.splitext(subim)[0]}.cat'
	subtbl = ascii.read(subcat)
	c_sub = SkyCoord(subtbl['ALPHA_J2000'], subtbl['DELTA_J2000'])
	#------------------------------------------------------------
	isubcat = f'{os.path.splitext(isubim)[0]}.cat'
	ccdinfo = get_ccdinfo(sciim, ccdtbl)
	param_insex4isub = dict(
		#------------------------------
		#	CATALOG
		#------------------------------
		CATALOG_NAME = isubcat,
		#------------------------------
		#	CONFIG FILES
		#------------------------------
		CONF_NAME = path_conf_iv,
		PARAMETERS_NAME = path_param_iv,
		FILTER_NAME = path_conv_iv,    
		STARNNW_NAME = path_nnw_iv,
		#------------------------------
		#	PHOTOMETRY
		#------------------------------
		GAIN = str(ccdinfo['gain'].value),
		PIXEL_SCALE = str(ccdinfo['pixscale'].value),
		#------------------------------
		#	STAR/GALAXY SEPARATION
		#------------------------------
		SEEING_FWHM = str(seeing.value),
		#------------------------------
		#	EXTRACTION
		#------------------------------
		DETECT_MINAREA = gphot_dict['DETECT_MINAREA'],
		DETECT_THRESH = gphot_dict['DETECT_THRESH'],
		DEBLEND_NTHRESH = gphot_dict['DEBLEND_NTHRESH'],
		DEBLEND_MINCONT = gphot_dict['DEBLEND_MINCONT'],
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
	)
	os.system(sexcom(isubim, param_insex4isub))
	isubtbl = ascii.read(isubcat)
	c_isub = SkyCoord(isubtbl['ALPHA_J2000'], isubtbl['DELTA_J2000'])
	#------------------------------------------------------------
	irefcat = f'{os.path.splitext(irefim)[0]}.cat'
	param_insex4iref = dict(
		#------------------------------
		#	CATALOG
		#------------------------------
		CATALOG_NAME = irefcat,
		#------------------------------
		#	CONFIG FILES
		#------------------------------
		CONF_NAME = path_conf_iv,
		PARAMETERS_NAME = path_param_iv,
		FILTER_NAME = path_conv_iv,    
		STARNNW_NAME = path_nnw_iv,
		#------------------------------
		#	PHOTOMETRY
		#------------------------------
		GAIN = str(ccdinfo['gain'].value),
		PIXEL_SCALE = str(ccdinfo['pixscale'].value),
		#------------------------------
		#	STAR/GALAXY SEPARATION
		#------------------------------
		SEEING_FWHM = str(seeing.value),
		#------------------------------
		#	EXTRACTION
		#------------------------------
		DETECT_MINAREA = gphot_dict['DETECT_MINAREA'],
		DETECT_THRESH = gphot_dict['DETECT_THRESH'],
		DEBLEND_NTHRESH = gphot_dict['DEBLEND_NTHRESH'],
		DEBLEND_MINCONT = gphot_dict['DEBLEND_MINCONT'],
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
	)
	os.system(sexcom(irefim, param_insex4iref))
	ireftbl = ascii.read(irefcat)
	c_iref = SkyCoord(ireftbl['ALPHA_J2000'], ireftbl['DELTA_J2000'])
	#------------------------------------------------------------
	indx_isub, sep_isub, _ = c_sub.match_to_catalog_sky(c_isub)
	indx_iref, sep_iref, _ = c_sub.match_to_catalog_sky(c_iref)

	subtbl['sep_isub'] = sep_isub
	subtbl['sep_iref'] = sep_iref
	'''
	plt.close('all')
	plt.hist(sep_isub.arcsec, bins=np.arange(0, 10+0.5, 0.5), alpha=0.5, label='sub')
	plt.hist(sep_iref.arcsec, bins=np.arange(0, 10+0.5, 0.5), alpha=0.5, label='ref')
	plt.axvline(x=seeing.value, ls='--', color='grey', alpha=0.5)
	plt.legend(fontsize=20)
	plt.show()
	'''
	#%%
	#============================================================
	#	Transient Filtering (FLAGING)
	#------------------------------------------------------------
	#	FLAG 0-12
	#------------------------------------------------------------
	# import numpy as np
	numbers = np.arange(0, 12+1, 1)
	for num in numbers: subtbl[f'flag_{num}'] = False
	subtbl['flag'] = False
	#------------------------------------------------------------
	#	Criterion
	#------------------------------------------------------------
	#	flag 0		: SkyBoT
	frac_sb = float(search_dict['frac_sb'])
	#	flag 1+2	: Invert source matching
	frac_inv = float(search_dict['frac_inv'])
	#	flag 3		: Sources @edge
	frac = float(search_dict['frac'])
	#	flag 4		: Upperlimit
	inmagkey = str(search_dict['inmagkey'])
	inul5key = str(search_dict['inul5key'])
	#	flag 5		: Ratio of ellipticity
	cutline_rellip = float(search_dict['cutline_rellip'])
	#	flag 6		: FLAGS
	cutline_flags = float(search_dict['cutline_flags'])
	#	flag 7		: seeing
	seeing_up = float(search_dict['seeing_up'])
	seeing_lo = float(search_dict['seeing_lo'])
	#	flag 8		: BKG
	cutline_bkg = float(search_dict['cutline_bkg'])
	#	flag 9		: Mag diff
	key0 = str(search_dict['key0'])
	key1 = str(search_dict['key1'])
	n_delm = float(search_dict['n_delm'])
	#	flag 10		: Dim pixel value
	frac_dim = float(search_dict['frac_dim'])
	#	flag 11		: NULL pixel
	cutline_null = float(search_dict['cutline_null'])
	'''
	#	flag 0		: SkyBoT
	frac_sb = 5.0
	#	flag 1+2	: Invert source matching
	frac_inv = 1.0
	#	flag 3		: Sources @edge
	frac = 0.9
	#	flag 4		: Upperlimit
	inmagkey = 'mag_aper'
	inul5key = 'ul5'
	#	flag 5		: Ratio of ellipticity
	cutline_rellip = 5
	#	flag 6		: FLAGS
	cutline_flags = 4.0
	#	flag 7		: seeing
	seeing_up = 3.0
	seeing_lo = 0.5
	#	flag 8		: BKG
	cutline_bkg = 50
	#	flag 9		: Mag diff
	key0 = 'mag_aper_1'
	key1 = 'mag_aper_3'
	n_delm = 10
	#	flag 10		: Dim pixel value
	frac_dim = 0.3
	#	flag 11		: NULL pixel
	cutline_null = 0
	'''
	#------------------------------------------------------------
	#	flag 0
	#------------------------------------------------------------
	#	Skybot query
	from astroquery.imcce import Skybot
	try:
		sbtbl = Skybot.cone_search(c_cent, ccdinfo['fov'], epoch)
		c_sb = SkyCoord(sbtbl['RA'], sbtbl['DEC'])
		sbtbl['sep'] = c_cent.separation(c_sb).to(u.arcmin)
		
		#	Skybot matching
		indx_sb, sep_sb, _ = c_sub.match_to_catalog_sky(c_sb)
		subtbl['sep_skybot'] = sep_sb
		subtbl['flag_0'][
			(sep_sb<seeing*frac_sb)
			] = True
	except:
		#	RuntimeError: No solar system object was found in the requested FOV
		print(f"No solar system object was found in the requested FOV ({ccdinfo['fov']})")
		pass
	#------------------------------------------------------------
	#	flag 1
	#------------------------------------------------------------
	if len(isubtbl)>0:
		#	Matching with inverted images
		indx_isub, sep_isub, _ = c_sub.match_to_catalog_sky(c_isub)
		subtbl['flag_1'][
			(sep_isub<seeing*frac_inv)
			] = True
	else:
		print('Inverted subtraction image has no source. ==> pass flag1')
		pass
	#------------------------------------------------------------
	#	flag 2
	#------------------------------------------------------------
	if len(ireftbl)>0:
		#	Coordinate
		c_iref = SkyCoord(ireftbl['ALPHA_J2000'], ireftbl['DELTA_J2000'], unit=u.deg)
		#	Matching with inverted images
		indx_iref, sep_iref, _ = c_sub.match_to_catalog_sky(c_iref)
		subtbl['flag_2'][
			(sep_iref<seeing*frac_inv)
			] = True
	else:
		print('Inverted reference image has no source. ==> pass flag2')
		pass
	#------------------------------------------------------------
	#	SEtractor criterion
	#------------------------------------------------------------
	#	flag 3
	#------------------------------------------------------------
	#	IMAGE PHYSICAL CENTER
	a = hdr['naxis1']/2.
	b = hdr['naxis2']/2.
	#	Small squre based on frac
	import numpy as np
	a_ = a*np.sqrt(frac)
	b_ = b*np.sqrt(frac)
	subtbl['flag_3'][
		(
			(subtbl['X_IMAGE']<a-a_) |
			(subtbl['X_IMAGE']>a+a_) |
			(subtbl['Y_IMAGE']<b-b_) |
			(subtbl['Y_IMAGE']>b+b_)
			)
		] = True
	"""
	plt.close('all')
	plt.figure(figsize=(5,5))
	fg3tbl = subtbl[subtbl['flag_3']==True]
	plt.plot(subtbl['ALPHA_J2000'], subtbl['DELTA_J2000'], marker='o', ls='none', alpha=0.5, mfc='none', mec='k')
	plt.plot(fg3tbl['ALPHA_J2000'], fg3tbl['DELTA_J2000'], marker='x', ls='none', alpha=1.0)
	xl, xr = plt.xlim()
	plt.xlim([xr, xl])"""
	#------------------------------------------------------------
	#	flag 4
	#------------------------------------------------------------
	#	More than 5 sigma signal
	subtbl['flag_4'][
		(subtbl[inmagkey]>hdr[inul5key])
		] = True
	"""
	plt.close('all')
	fg4tbl = subtbl[subtbl['flag_4']==True]
	plt.plot(subtbl[inmagkey], subtbl['magerr_aper'], marker='o', ls='none', alpha=0.5, mfc='none', mec='k')
	plt.plot(fg4tbl[inmagkey], fg4tbl['magerr_aper'], marker='x', ls='none', alpha=1.0)
	plt.axvline(hdr['ul5'], ls='--', c='tomato')
	plt.xlim([np.min(subtbl[inmagkey]-2.5), hdr['ul5']+5])
	plt.ylim([-0.25, 2.5])
	plt.grid('both', ls='--', c='silver', alpha=0.5)"""
	#	Empirical criterion
	#------------------------------------------------------------
	#	flag 5
	#------------------------------------------------------------
	subtbl['ratio_ellipticity'] = subtbl['ELLIPTICITY']/hdr['ELLIP']
	subtbl['ratio_elongation'] = subtbl['ELONGATION']/hdr['ELONG']

	subtbl['flag_5'][
		(subtbl['ratio_ellipticity'] > cutline_rellip)
		] = True
	"""
	plt.close('all')
	fg5tbl = subtbl[subtbl['flag_5']==True]
	bins = np.arange(np.min(subtbl['ratio_ellipticity']), np.max(subtbl['ratio_ellipticity'])+0.1, 0.1)
	plt.hist(subtbl['ratio_ellipticity'], bins=bins, alpha=0.5,)
	plt.hist(fg5tbl['ratio_ellipticity'], bins=bins, alpha=0.75)
	plt.axvline(x=5, ls='--', c='tomato')
	plt.grid('both', ls='--', c='silver', alpha=0.5)"""
	#------------------------------------------------------------
	#	flag 6
	#------------------------------------------------------------
	subtbl['flag_6'][
		(subtbl['FLAGS'] > cutline_flags)
		] = True
	"""
	plt.close('all')
	fg6tbl = subtbl[subtbl['flag_6']==True]
	bins = np.arange(np.min(subtbl['FLAGS']), np.max(subtbl['FLAGS'])+0.5, 0.5)
	plt.hist(subtbl['FLAGS'], bins=bins, alpha=0.5,)
	plt.hist(fg6tbl['FLAGS'], bins=bins, alpha=0.75)
	plt.axvline(x=4, ls='--', c='tomato')
	plt.grid('both', ls='--', c='silver', alpha=0.5)"""
	#------------------------------------------------------------
	#	flag 7
	#------------------------------------------------------------
	subtbl['ratio_seeing'] = subtbl['FWHM_WORLD'].to(u.arcsec)/seeing
	subtbl['flag_7'][
		(subtbl['ratio_seeing']>seeing_up) |
		(subtbl['ratio_seeing']<seeing_lo)
		] = True
	"""
	#	Check plot
	plt.close('all')
	fg7tbl = subtbl[subtbl['flag_7']==True]
	bins = np.arange(np.min(subtbl['ratio_seeing']), np.max(subtbl['ratio_seeing'])+0.25, 0.25)
	plt.hist(subtbl['ratio_seeing'], bins=bins, alpha=0.5,)
	plt.hist(fg7tbl['ratio_seeing'], bins=bins, alpha=0.75)

	plt.axvline(x=1, ls='-', c='red')
	plt.axvline(x=seeing_lo, ls='--', c='tomato')
	plt.axvline(x=seeing_up, ls='--', c='tomato')
	plt.grid('both', ls='--', c='silver', alpha=0.5)"""
	#------------------------------------------------------------
	#	flag 8
	#------------------------------------------------------------
	subtbl['flag_8'][
		(subtbl['BACKGROUND']<-cutline_bkg) |
		(subtbl['BACKGROUND']>+cutline_bkg)
		] = True
	"""
	#	Check plot
	plt.close('all')
	fg8tbl = subtbl[subtbl['flag_8']==True]
	bins = np.arange(np.min(subtbl['BACKGROUND']), np.max(subtbl['BACKGROUND'])+1, 1)
	plt.hist(subtbl['BACKGROUND'], bins=bins, alpha=0.5,)
	plt.hist(fg8tbl['BACKGROUND'], bins=bins, alpha=0.75)

	plt.axvline(x=-cutline_bkg, ls='--', c='tomato')
	plt.axvline(x=+cutline_bkg, ls='--', c='tomato')
	plt.grid('both', ls='--', c='silver', alpha=0.5)"""
	#------------------------------------------------------------
	#	flag 9
	#------------------------------------------------------------
	scitbl = ascii.read(scicat)
	scitbl = scitbl[
		(scitbl['FLAGS']==0) &
		(scitbl['CLASS_STAR']>0.9)
	]

	aperdict = {
		'mag_aper':'SNR_curve',
		'mag_aper_1':'Best_Aperture',
		'mag_aper_2':'2seeing',
		'mag_aper_3':'3seeing',
		'mag_aper_4':'3arcsec',
		'mag_aper_5':'5arcsec',	
	}

	#	Sci. sources magnitude diff.
	scidelm = scitbl[key0] - scitbl[key1]

	#	Subt. sources magnitude diff.
	subdelm = subtbl[key0] - subtbl[key1]
	subtbl['del_mag'] = subdelm

	def get_mad(data):
		return bn.median(np.absolute(data - bn.median(data, axis=0)), axis=0)

	#	MED & MAD
	scidelm_med = bn.median(scidelm)
	scidelm_mad = get_mad(scidelm)
	subtbl['del_mag_med'] = scidelm_med
	subtbl['del_mag_mad'] = scidelm_mad
	subtbl['N_del_mag_mad'] = np.abs((subtbl['del_mag']-subtbl['del_mag_med'])/subtbl['del_mag_mad'])

	#	INDEXING
	indx_delm = np.where(
		(subdelm<scidelm_med-scidelm_mad*n_delm) |
		(subdelm>scidelm_med+scidelm_mad*n_delm)
		)
	subtbl['flag_9'][indx_delm] = True
	"""
	plt.close('all')
	bins = np.arange(-2.5, 2.5+0.125, 0.125)
	plt.hist(subdelm, bins=bins, alpha=0.5)
	plt.hist(subdelm[indx_delm], bins=bins, alpha=0.5)
	plt.axvline(x=scidelm_med-scidelm_mad*n_delm, ls='--', c='tomato')
	plt.axvline(x=scidelm_med+scidelm_mad*n_delm, ls='--', c='tomato')"""
	#------------------------------------------------------------
	#	flag 10+11
	#------------------------------------------------------------



	subtbl['n_dim'] = 0
	subtbl['ratio_dim'] = 0.0
	subtbl['n_null'] = 0
	subtbl['ratio_null'] = 0

	for i in range(len(subtbl)):
		tx, ty = subtbl['X_IMAGE'][i], subtbl['Y_IMAGE'][i]
		bkg = subtbl['BACKGROUND'][i]
		badsky = subtbl['BACKGROUND'] - skysig
		crt = badsky[i]
		#	Crop data
		cdata = crop_data(
			data=data_sub,
			tx=tx,
			ty=ty,
			tsize=peeing,
			)
		#	Count dim pixel
		(n_dim, ratio_dim) = count_dim_pixels(
			cdata=cdata,
			crt=crt,
			)
		subtbl['n_dim'][i] = n_dim
		subtbl['ratio_dim'][i] = ratio_dim
		#	Count null pixel
		n_null = len(np.where(cdata == 1e-30)[0])
		ratio_null = n_null/cdata.size
		subtbl['n_null'][i] = n_null
		subtbl['ratio_null'][i] = ratio_dim
	cutline_dim = cdata.size*frac_dim
	#	Flag for Dipole 
	subtbl['flag_10'][subtbl['ratio_dim']>cutline_dim] = True
	#	Flag for Null
	subtbl['flag_11'][subtbl['n_null']>cutline_null] = True
	#------------------------------------------------------------
	#	flag 12
	#------------------------------------------------------------
	x, y = fits.getdata(refim).shape
	w_ref = WCS(refim)
	xim, yim = w_ref.world_to_pixel(c_sub)
	indx_out = np.where(
		(xim < 0) | (xim > x) | (yim < 0) | (yim > y)
	)
	subtbl['flag_12'][indx_out] = True
	'''
	plt.close('all')
	plt.figure(figsize=(5,5))
	fg12tbl = subtbl[subtbl['flag_12']==True]
	plt.plot(subtbl['ALPHA_J2000'], subtbl['DELTA_J2000'], marker='o', ls='none', alpha=0.5, mfc='none', mec='k')
	plt.plot(fg12tbl['ALPHA_J2000'], fg12tbl['DELTA_J2000'], marker='x', ls='none', alpha=1.0)
	xl, xr = plt.xlim()
	plt.xlim([xr, xl])'''
	#------------------------------------------------------------
	#	Final flag
	#------------------------------------------------------------
	flag = subtbl['flag']
	n_all = len(subtbl)
	print(f"{'-'*60}\nFLAG SUMMARY\t\t({len(subtbl)})\n{'-'*60}")
	for n in numbers:
		if n<10:
			interval='\t\t\t'
		else:
			interval='\t\t'
		tmptbl = subtbl[subtbl[f'flag_{n}']==True]
		print(f'flag=={n}{interval}: {len(tmptbl)}\t({int(100*len(tmptbl)/n_all)}%)')
		flag = flag + subtbl[f'flag_{n}']
	subtbl['flag'] = flag
	#	SkyBoT exception
	indx_sb = np.where(subtbl['flag_0']==True)
	subtbl['flag'][indx_sb] = False
	#------------------------------------------------------------
	#	END
	#------------------------------------------------------------
	bgstbl = subtbl[subtbl[f'flag']==True]
	tctbl = subtbl[subtbl[f'flag']==False]
	#	Save catalog
	tccat = f"{os.path.splitext(subcat)[0]}.transients.cat"
	tctbl.write(tccat, format='ascii.ecsv', overwrite=True)
	subtbl.write(subcat, format='ascii.ecsv', overwrite=True)
	print('-'*60)
	print(f'Filtered sources\t: {len(bgstbl)}\t({int(100*len(bgstbl)/n_all)}%)')
	print(f'Transient Candidates\t: {len(tctbl)}\t({int(100*len(tctbl)/n_all)}%)')
	print('='*60)
	'''
	# %%
	plt.figure(figsize=(5,5))
	plt.imshow(data)
	plt.plot(subtbl['X_IMAGE'], subtbl['Y_IMAGE'], marker='o', ls='none', alpha=0.5, mfc='none', mec='k')
	plt.plot(bgstbl['X_IMAGE'], bgstbl['Y_IMAGE'], marker='x', ls='none', alpha=1.0)'''
#------------------------------------------------------------
def mp_search(sciimlist, refimlist, subimlist, search_dict, ncore=4):
	if __name__ == '__main__':
		#	Fixed the number of cores (=4)
		with multiprocessing.Pool(processes=ncore) as pool:
			results = pool.starmap(
				routine_search,
				zip(
					sciimlist,
					refimlist,
					subimlist,
					repeat(search_dict),
					)
				)
#------------------------------------------------------------
def set_normalization(inim, kind, cutsize):
	#------------------------------------------------------------
	#
	#------------------------------------------------------------
	#
	hdu = fits.open(inim)[0]
	wcs = WCS(hdu.header)
	seeing = hdu.header['SEEING']*u.arcsec
	#	size
	size = u.Quantity((cutsize, cutsize))
	size_small = u.Quantity((seeing, seeing))
	#	Make the cutout, including the WCS
	cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs)
	data = 	cutout.data
	scutout = Cutout2D(hdu.data, position=position, size=size_small, wcs=wcs)
	sdata = scutout.data
	#------------------------------------------------------------
	#	Norm
	#------------------------------------------------------------
	if kind == 'ref':
		norm = ImageNormalize(	
			vmin=np.median(data),
			vmax=np.median(data)*1e2,
			stretch=LogStretch(),
			interval=ZScaleInterval(),
			)
	#------------------------------------------------------------
	elif kind == 'sub':
		norm = ImageNormalize(	
			vmin=np.median(data),
			vmax=np.max(sdata),
			stretch=LogStretch(),
			interval=ZScaleInterval(),
			)	
	#------------------------------------------------------------
	elif kind == 'sci':
		norm = ImageNormalize(	
			vmin=np.median(data),
			vmax=np.max(sdata),
			stretch=LogStretch(),
			interval=ZScaleInterval(),
			)
	else:
		#	Something goes wrong
		norm = ImageNormalize(	
			vmin=np.median(data),
			vmax=np.max(sdata),
			stretch=LogStretch(),
			interval=ZScaleInterval(),
			)
	return norm, data
#------------------------------------------------------------
def make_snapshot(data, wcs, norm, peeing, outpng):
	plt.close('all')
	plt.rc('font', family='serif')
	fig = plt.figure(figsize=(1, 1))
	fig.set_size_inches(1. * data.shape[0] / data.shape[1], 1, forward = False)
	x = 720 / fig.dpi
	y = 720 / fig.dpi
	fig.set_figwidth(x)
	fig.set_figheight(y)
	#	No axes
	ax = plt.subplot(projection=wcs)
	ax = plt.Axes(fig, [0., 0., 1., 1.])
	ax.set_axis_off()
	fig.add_axes(ax)

	ax.imshow(
		data, 
		origin='lower',
		norm=norm,
		cmap='gist_gray',
		)
	#	Circle
	circle = Circle(
		(data.shape[0]/2., data.shape[1]/2.),
		2*peeing,
		edgecolor='yellow',
		facecolor=None,
		fill=False
	)
	ax.add_patch(circle)
	#	Option
	plt.gca().invert_xaxis()
	plt.gca().invert_yaxis()
	plt.savefig(outpng, dpi=200, overwrite=True)
#------------------------------------------------------------
def routine_snapshot(sciim, refim, subim, subcat, position, cutsize, n):
	imlist = [sciim, refim, subim]
	kinds = ['sci', 'ref', 'sub']
	for inim, kind in zip(imlist, kinds):
		#	Output name
		outfits = f"{os.path.splitext(subim)[0]}.{kind}.{n}.fits"
		outpng = f"{os.path.splitext(subim)[0]}.{kind}.{n}.png"

		hdu = fits.open(inim)[0]
		wcs = WCS(hdu.header)
		peeing = hdu.header['PEEING']
		#	Crop image and Get Norm
		norm, data = set_normalization(inim, kind, cutsize)
		fits.writeto(outfits, data, header=hdu.header, overwrite=True)
		#	Make Snapshot
		make_snapshot(
			data=data,
			wcs=wcs,
			norm=norm,
			peeing=peeing,
			outpng=outpng,
			)
#------------------------------------------------------------
def mp_snapshot(sciim, refim, subim, subcat, positions, cutsize, numbers, ncore=4):
	if __name__ == '__main__':
		#	Fixed the number of cores (=4)
		with multiprocessing.Pool(processes=ncore) as pool:
			results = pool.starmap(
				routine_snapshot,
				zip(
					repeat(sciim),
					repeat(refim),
					repeat(subim),
					repeat(subcat),
					positions,
					repeat(cutsize),
					numbers,
					)
				)
#------------------------------------------------------------
#	Photometry
#------------------------------------------------------------
mp_phot(tstbl['sub'], ncore=ncore)
#------------------------------------------------------------
#	Transient Filtering
#------------------------------------------------------------
mp_search(
	sciimlist=tstbl['sci'],
	refimlist=tstbl['ref'],
	subimlist=tstbl['sub'],
	search_dict=search_dict,
	ncore=ncore,
	)
#------------------------------------------------------------
#	Snapshot
#------------------------------------------------------------
print(f"\n{'='*60}\nMAKE SNAPSHOTs ({len(tstbl)} IMAGES)\n{'-'*60}")
cutsize = 3.0*u.arcmin
for n in range(len(tstbl)):
	print(f"[{n}]", end=' ')
	sciim = tstbl['sci'][n]
	refim = tstbl['ref'][n]
	subim = tstbl['sub'][n]

	subcat = f'{os.path.splitext(subim)[0]}.cat'
	tccat = f"{os.path.splitext(subcat)[0]}.transients.cat"
	tctbl = ascii.read(tccat)

	positions = []
	for i in range(len(tctbl)):
		tra = tctbl['ALPHA_J2000'][i].item()
		tdec = tctbl['DELTA_J2000'][i].item()
		position = SkyCoord(tra, tdec, frame='icrs', unit='deg')
		positions.append(position)

	mp_snapshot(
		sciim=sciim,
		refim=refim,
		subim=subim,
		subcat=subcat,
		positions=positions,
		cutsize=cutsize,
		numbers=tctbl['NUMBER'],
		ncore=ncore,
		)
	print("DONE")
print('='*60)
#------------------------------------------------------------
print("ALL DONE!")