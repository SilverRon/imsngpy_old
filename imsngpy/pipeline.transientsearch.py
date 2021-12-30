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
#------------------------------------------------------------
#	Configuration files
#------------------------------------------------------------
#	gpphot
prefix_gp = 'gpphot'
path_param = f'{path_config}/{prefix_gp}.param'
path_conv = f'{path_config}/{prefix_gp}.conv'
path_nnw = f'{path_config}/{prefix_gp}.nnw'
path_conf = f'{path_config}/{prefix_gp}.sex'
#============================================================
#	Input
#------------------------------------------------------------
try:
	imkey = sys.argv[1]
	# imkey = '/data3/paek/factory/loao/test_fast/Calib-LOAO-*com.fits'
	print(f'IMKEY  : {imkey}')
except:
	# imkey = '/data3/paek/factory/loao/test/Calib-*0.fits'
	imkey = '/data3/paek/factory/loao/test/Calib-*com.fits'
	print(f'Use default IMKEY : {imkey}')
try:
	ncore = int(sys.argv[2])
except:
	ncore = 4
print(f'NCORE  : {ncore}')
#------------------------------------------------------------
#	Function
#------------------------------------------------------------
def routine_phot(inim,):
	print(f'Photometry is started for {os.path.basename(inim)}')
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
	print(f"""{'-'*60}\n#\tCCD INFO\n{'-'*60}""")
	gain = ccdtbl['gain'][indx_ccd][0]*(u.electron/u.adu)
	rdnoise = ccdtbl['readnoise'][indx_ccd][0]*(u.electron)
	pixscale = ccdtbl['pixelscale'][indx_ccd][0]*(u.arcsec/u.pixel)
	fov = ccdtbl['foveff'][indx_ccd][0]*(u.arcmin)
	print(f"""GAIN : {gain}\nREAD NOISE : {rdnoise}\nPIXEL SCALE : {pixscale}\nEffective FoV : {fov}""")
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
	print(f"{'='*60}")
	print('\nZP CALCULATION for each APERTURES\n')
	print(f"{'-'*60}")
	for n, inmagkey in enumerate(list(aper_dict.keys())):
		print(f"[{n+1}/{len(list(aper_dict.keys()))}] {inmagkey}")
		inerrkey = aper_dict[inmagkey]['errkey']
		inzpkey = f"ZP_{aper_dict[inmagkey]['suffix']}"
		inul3key = f"UL3{aper_dict[inmagkey]['suffix'].replace('APER', '')}"
		inul5key = f"UL5{aper_dict[inmagkey]['suffix'].replace('APER', '')}"
		hdrzpkey = inzpkey.replace('_APER', '')
		hdrzperkey = hdrzpkey.replace('ZP', 'ZPER')
		hdrnstdkey = hdrzpkey.replace('ZP', 'NZP')

		#	FILTER zptbl_ WITH MAGERR
		zp = hdr[hdrzpkey]
		zper = hdr[hdrzperkey]

		#	MERGED TABLE
		rawtbl[inmagkey.lower()] = rawtbl[inmagkey]+zp
		rawtbl[aper_dict[inmagkey]['errkey'].lower()] = sqsum(rawtbl[inerrkey], zper)

	delt = time.time()-st_
	print(f"PHOTOMETRY IS DONE for {os.path.basename(inim)} ({round(delt, 1)} sec)")
#------------------------------------------------------------
def mp_phot(imlist, ncore=4):
	if __name__ == '__main__':
		#	Fixed the number of cores (=4)
		with multiprocessing.Pool(processes=ncore) as pool:
			results = pool.starmap(
				routine,
				zip(
					imlist,
					)
				)
#------------------------------------------------------------
def phot(imlist):
	for i, inim in enumerate(imlist):
		print(f"[{i+1}/{len(imlist)}] {os.path.basename(inim)}")
		routine(inim)
#------------------------------------------------------------
imlist = sorted(glob.glob(imkey))
#	Photometry
mp_phot(imlist, ncore=ncore)