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
#	SNRCURVE
prefix_sc = 'snrcurve'
path_param_sc = f'{path_config}/{prefix_sc}.param'
path_conv_sc = f'{path_config}/{prefix_sc}.conv'
path_nnw_sc = f'{path_config}/{prefix_sc}.nnw'
path_conf_sc = f'{path_config}/{prefix_sc}.sex'
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
	# imkey = gphot_dict['imkey']
	# imkey = './Calib*0.fits'
	# imkey = '/data3/paek/factory/loao/test/Calib-LOAO-NGC3147-20210421-*0.fits'
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
def routine(inim,):
	print(f'Photometry is started for {os.path.basename(inim)}')
	# global gphot_dict
	#
	global path_param_sc
	global path_conv_sc
	global path_nnw_sc
	global path_conf_sc
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
	#	SNR CURVE
	#------------------------------------------------------------
	#	Aperture steps for SNR curve
	apertures = peeing.value*np.linspace(0.25, 5, 32)
	apertures_input = ','.join([str(d) for d in apertures])
	#------------------------------------------------------------
	#	Run source extractor for SNR Curve
	#------------------------------------------------------------
	sccat = f'{os.path.splitext(inim)[0]}.sc.cat'
	param_insex = dict(
		#------------------------------
		#	CATALOG
		#------------------------------
		CATALOG_NAME = sccat,
		#------------------------------
		#	CONFIG FILES
		#------------------------------
		CONF_NAME = path_conf_sc,
		PARAMETERS_NAME = path_param_sc,
		FILTER_NAME = path_conv_sc,    
		STARNNW_NAME = path_nnw_sc,
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
	os.system(sexcom(inim, param_insex))
	scrawtbl = ascii.read(sccat)
	#------------------------------------------------------------
	#	Conditions for growth-curve
	##	Point sources and nearby center
	#------------------------------------------------------------
	sctbl = scrawtbl[
		(scrawtbl['FLAGS']==0) &
		(scrawtbl['CLASS_STAR']>0.9) &
		(scrawtbl['X_IMAGE']>a-a_) & (scrawtbl['X_IMAGE']<a+a_) &
		(scrawtbl['Y_IMAGE']>b-b_) & (scrawtbl['Y_IMAGE']<b+b_)
	]
	#	Generate SNR column
	sctbl['APER_OPT'] = 0.0
	for n in range(len(apertures)):
		if n==0:
			sctbl[f'SNR'] = sctbl[f'FLUX_APER']/sctbl[f'FLUXERR_APER']
		else:
			sctbl[f'SNR_{n}'] = sctbl[f'FLUX_APER_{n}']/sctbl[f'FLUXERR_APER_{n}']
	#%%
	sctbl['APER_OPT'] = generate_best_aperture_with_snrcurve(sctbl, apertures, pixscale)
	aper_opt = np.median(sctbl['APER_OPT'])
	#------------------------------------------------------------
	#%%
	#	OPTIMIZED APERTURE FROM SNR CURVE
	#------------------------------------------------------------
	'''
	draw_snrcurve(
		title = os.path.basename(inim),
		scoutpng = f'{os.path.splitext(inim)[0]}.sc.png',
		seeing=seeing,
		aper_opt=aper_opt,
		dpi=200,
		)'''
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
	#	PRE-SELECTION FOR MATCHING
	##	Select point sources nearby center
	#------------------------------------------------------------
	#	Image center
	w = WCS(inim)
	imcent = w.all_pix2world(a, b, 1)
	raim = imcent[0].item()
	decim = imcent[1].item()
	c_cent = SkyCoord(raim, decim, unit=u.deg)
	#	Catalog center
	c_raw = SkyCoord(rawtbl['ALPHA_J2000'], rawtbl['DELTA_J2000'], unit=u.deg)
	#	Seperation
	sep = c_raw.separation(c_cent)
	sepcut = fov/2*frac
	#	Select 
	indx_sel, indxes = select_point_sources(
		rawtbl=rawtbl,
		errkey=aper_dict['MAG_APER']['errkey'],
		sep=sep,
		classstarcut=0.5,
		flagcut=float(gphot_dict['flagcut']),
		# magerrcut=float(gphot_dict['inmagerupper']),
		sepcut=sepcut,
		)
	#	Print
	#%%
	print(f"{'='*60}")
	print(f'ALL                  : {len(rawtbl)}')
	print(f"{'-'*60}")
	print(f'CLASS_STAR > {0.5}     : {len(indxes[0][0])} ({round(1e2*len(indxes[0][0])/len(rawtbl), 1)}%)')
	print(f"FLAGS <= {gphot_dict['flagcut']}           : {len(indxes[1][0])} ({round(1e2*len(indxes[1][0])/len(rawtbl), 1)}%)")
	# print(f"{aper_dict['MAG_APER']['errkey']} < {gphot_dict['inmagerupper']} : {len(indxes[2][0])} ({round(1e2*len(indxes[2][0])/len(rawtbl), 1)}%)")
	print(f"SEP < {round(sepcut.to(u.arcmin).value, 1)}'          : {len(indxes[2][0])} ({round(1e2*len(indxes[2][0])/len(rawtbl), 1)}%)")
	print(f"{'='*60}")
	print(f"TOTAL                : {len(indx_sel[0])} ({round(1e2*len(indx_sel[0])/len(rawtbl), 1)}%)")
	print(f"{'='*60}")
	#%%
	#------------------------------------------------------------
	#	QUERY REFERENCE CATALOG
	#------------------------------------------------------------
	seltbl = rawtbl[indx_sel]
	refcatname = gphot_dict['refcatname']
	outrefcat = f'{os.path.dirname(inim)}/ref.{obj}.{refcatname}.ecsv'
	if os.path.exists(outrefcat):
		#	If exist, skip query
		print(f'Found {outrefcat}')
		reftbl = ascii.read(outrefcat)
	else:
		reftbl = querybox(
			refcatname=refcatname,
			racent=c_cent.ra.value,
			decent=c_cent.dec.value,
			# outcat=f'{os.path.dirname(inim)}/ref.{obj}.{refcatname}.cat',
			outcat=outrefcat,
			radius=float(gphot_dict['refqueryradius']),
			#	Broad-band  : ''
			#	Median-band : 'med'
			refmagkey=''
			)
	#%%
	'''
	draw_space_distribution(
		outpng=f'{os.path.splitext(inim)[0]}.space.png',
		c_cent=c_cent,
		rawtbl=rawtbl,
		seltbl=seltbl,
		reftbl=reftbl,
		dpi=200,
		)
	'''
	#------------------------------------------------------------
	#%%
	#	Matching catalog
	#------------------------------------------------------------
	c_ref = SkyCoord(reftbl['ra'], reftbl['dec'], unit=u.deg)
	c = SkyCoord(seltbl['ALPHA_J2000'], seltbl['DELTA_J2000'], unit=u.deg)
	#	REFERENCE CATALOG FILTERING
	#	Magnitude cut and error cut
	indx, sep, _ = c.match_to_catalog_sky(c_ref)
	mtbl = hstack([seltbl, reftbl[indx]])
	mtbl['sep'] = sep
	#	Seperation distribution
	'''
	plt.close('all')
	plt.hist(sep.arcsec*1e2/seeing.value, bins=np.arange(0, 100+1, 1), histtype='step', color='k')
	plt.axvline(x=np.median(sep.arcsec*1e2/seeing.value), ls='--', c='tomato', label=f'{round(np.median(sep.arcsec*1e2/seeing.value), 1)}%')
	plt.xlabel('SEP/SEEING [%]')
	plt.ylabel('#')
	plt.legend(fontsize=20)
	'''
	#%%
	#------------------------------------------------------------
	#	Point sources for zeropoint
	#------------------------------------------------------------
	indx_zp, indxes_zp = select4zp(
		mtbl=mtbl,
		filte=filte,
		sepcut=seeing*float(gphot_dict['sepfrac']),
		refmagerrcut=float(gphot_dict['refmagerupper']),
		inmagerrcut=float(gphot_dict['inmagerupper']),
		refmaglowercut=float(gphot_dict['refmaglower']),
		refmaguppercut=float(gphot_dict['refmagupper']),
		)
	zptbl_ = mtbl[indx_zp]
	#	SELECTION SUMMARY
	print(f"{'='*60}")
	print(f'ALL                  : {len(mtbl)}')
	print(f"{'-'*60}")
	print(f"SEP < SEEING*{float(gphot_dict['sepfrac'])}     : {len(indxes[0][0])} ({round(1e2*len(indxes_zp[0][0])/len(mtbl), 1)}%)")
	print(f"{gphot_dict['refmaglower']} < {filte} < {gphot_dict['refmagupper']}          : {len(indxes_zp[4][0])} ({round(1e2*len(indxes_zp[4][0])/len(mtbl), 1)}%)")
	print(f"{filte}err < {gphot_dict['refmagerupper']}          : {len(indxes[1][0])} ({round(1e2*len(indxes_zp[1][0])/len(mtbl), 1)}%)")
	print(f"No masked {filte}          : {len(indxes_zp[5][0])} ({round(1e2*len(indxes_zp[5][0])/len(mtbl), 1)}%)")
	# print(f"No masked {inmagkey}   : {len(indxes_zp[5][0])} ({round(1e2*len(indxes_zp[5][0])/len(mtbl), 1)}%)")
	print(f"{'-'*60}")
	print(f"*{filte} < {gphot_dict['refmaglower']}              : {len(indxes_zp[2][0])} ({round(1e2*len(indxes_zp[2][0])/len(mtbl), 1)}%)")
	print(f"*{filte} > {gphot_dict['refmagupper']}              : {len(indxes_zp[3][0])} ({round(1e2*len(indxes_zp[3][0])/len(mtbl), 1)}%)")
	print(f"{'='*60}")
	print(f"TOTAL                : {len(indx_zp[0])} ({round(1e2*len(indx_zp[0])/len(mtbl), 1)}%)")
	print(f"{'='*60}")
	#------------------------------------------------------------
	#	GENERAL INFO. REGARDLESS OF APERTURE
	#------------------------------------------------------------
	#	VALUES
	# elong = bn.median(mtbl['ELONGATION'][~czparr.mask])
	# ellip = bn.median(mtbl['ELLIPTICITY'][~czparr.mask])
	# skyval = bn.median(mtbl['BACKGROUND'][~czparr.mask])
	# skysig = bn.nanstd(mtbl['BACKGROUND'][~czparr.mask])
	# rsky = 1e2*skysig/skyval
	elong = bn.median(mtbl['ELONGATION'])
	ellip = bn.median(mtbl['ELLIPTICITY'])
	# skyval = bn.median(mtbl['BACKGROUND'])
	# skysig = bn.nanstd(mtbl['BACKGROUND'])
	rsky = 1e2*skysig/skyval
	#	PUT HEADER ON THE IMAGE
	fits.setval(inim, 'AUTHOR', value='Gregory S.H. Paek', comment='PHOTOMETRY AUTHOR')
	fits.setval(inim, 'PHTDATE', value=date.today().isoformat(), comment='PHOTOMETRY DATE [KOR]')
	fits.setval(inim, 'ELLIP', value=round(ellip, 3), comment='MEDIAN ELLIPTICITY 1-B/A [0-1]')
	fits.setval(inim, 'ELONG', value=round(elong, 3), comment='MEDIAN ELONGATION A/B [1-]')
	fits.setval(inim, 'SKYSIG', value=round(skysig, 3), comment='STD of BACKGROUND [pixel]')
	fits.setval(inim, 'SKYVAL', value=round(skyval, 3), comment='MEDIAN of BACKGROUND [pixel]')
	fits.setval(inim, 'RSKY', value=round(rsky, 3), comment='RATIO of SKYSIG/SKYVAL [%]')
	fits.setval(inim, 'REFCAT', value=gphot_dict['refcatname'], comment='REFERENCE CATALOG NAME')
	fits.setval(inim, 'MAGLOW', value=float(gphot_dict['refmaglower']), comment='REFERENCE MAG BRIGHT LIMIT')
	fits.setval(inim, 'MAGUP', value=float(gphot_dict['refmagupper']), comment='REFERENCE MAG DIM LIMIT')
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
		zptbl_[inzpkey] = zptbl_[filte]-zptbl_[inmagkey]
		zptbl = zptbl_[zptbl_[inerrkey]<float(gphot_dict['inmagerupper'])]
		print(f"{'='*60}")
		print(f"{inerrkey} < {float(gphot_dict['inmagerupper'])}")
		print(f"{'-'*60}")
		print(f'ALL : {len(zptbl_)} --> SELECT : {len(zptbl)} ({round(1e2*len(zptbl)/len(zptbl_), 1)} %)')
		# print(f"{'-'*60}")
		#%%
		#------------------------------------------------------------
		#	Sigma-Clipping
		#------------------------------------------------------------
		# sigma=1.5
		sigma=2.0
		iteration=5
		czparr = sigma_clip(
			np.copy(zptbl[inzpkey]),
			sigma=sigma,
			# maxiters=None,
			# maxiters=3,
			maxiters=iteration,
			# cenfunc=np.median,
			cenfunc=bn.nanmedian,
			copy=False
			)
		zp = np.median(czparr.data[~czparr.mask])
		zper = np.std(czparr.data[~czparr.mask])
		#
		print(f"{'-'*60}")
		print(f"SIGMA CLIPPING {len(czparr)} --> {len(czparr[~czparr.mask])} ({round(1e2*len(czparr[~czparr.mask])/len(czparr), 1)} %)")
		print(f"- SIGMA     : {sigma}")
		print(f"- ITERATION : {iteration}")
		print(f"{'-'*60}")
		print(f"{inzpkey} = {round(zp, 3)} +/- {round(zper, 3)}")
		print(f"{'='*60}\n")
		#	MERGED TABLE
		mtbl[inmagkey.lower()] = mtbl[inmagkey]+zp
		# mtbl[aper_dict[inmagkey]['errkey'].lower()] = sqsum(mtbl[inerrkey], mtbl[f'{filte}err'])
		mtbl[aper_dict[inmagkey]['errkey'].lower()] = sqsum(mtbl[inerrkey], zper)
		#%%
		#------------------------------------------------------------
		#	ZP CALCULATION SUMMARY
		#------------------------------------------------------------
		#	ZP TABLE
		zptbl[inmagkey.lower()] = zptbl[inmagkey]+zp
		zptbl[aper_dict[inmagkey]['errkey'].lower()] = sqsum(zptbl[inerrkey], zptbl[f'{filte}err'])
		draw_zpcal(
			outpng=f'{os.path.splitext(inim)[0]}.{hdrzpkey.lower()}.png',
			zptbl=zptbl,
			# mtbl=mtbl,
			filte=filte,
			czparr=czparr,
			indxes_zp=indxes_zp,
			inzpkey=inzpkey,
			refmaglowercut=float(gphot_dict['refmaglower']),
			refmaguppercut=float(gphot_dict['refmagupper']),
			dpi=200,
			)
		#%%
		'''if inmagkey == 'MAG_AUTO':
			#	DEPTH FOR MAG_AUTO
			plt.plot(mtbl[inmagkey.lower()], mtbl[aper_dict[inmagkey]['errkey'].lower()], ls='none', marker='.', mec='k', mfc='none', alpha=0.25)
			plt.axhline(y=1./5., ls='--', c='tomato')
			plt.ylim([0.0, 0.5])'''
		#%%
		#------------------------------------------------------------
		#	SPATIAL DISTRIBUTION
		#------------------------------------------------------------
		
		draw_zp_space(
			inim,
			outpng=f'{os.path.splitext(inim)[0]}.{inzpkey.lower()}.png',
			c_cent=c_cent,
			inzpkey=inzpkey,
			zptbl=zptbl,
			mask=czparr.mask,
			dpi=100,
			)
		
		#%%
		#	DEPTH / LIMITING MAGNITUDE
		if inmagkey == 'MAG_AUTO':
			ul3, ul5 = 0, 0
		else:
			#	MAG_APER			
			ul3 = calc_depth(
				N=3,
				zp=zp,
				aperture=float(aper_dict[inmagkey]['size']),
				skysig=skysig
				)
			ul5 = calc_depth(
				N=5,
				zp=zp,
				aperture=float(aper_dict[inmagkey]['size']),
				skysig=skysig
				)
		#%%
		#------------------------------------------------------------
		#	HEADER INFO
		#------------------------------------------------------------
		#	For each apertures
		fits.setval(inim, aper_dict[inmagkey]['suffix'], value=round(aper_dict[inmagkey]['size'], 3), comment=aper_dict[inmagkey]['comment'])
		fits.setval(inim, hdrnstdkey, value=len(zptbl[~czparr.mask]), comment=f'NUMBER OF STARS to MEASURE {hdrzpkey}')
		fits.setval(inim, hdrzpkey, value=round(zp, 3), comment=f'ZEROPOINT for {inmagkey}')
		fits.setval(inim, hdrzperkey, value=round(zper, 3), comment=f'ZEROPOINT for {inmagkey}')
		fits.setval(inim, inul3key, value=round(ul3, 3), comment=f'3 SIGMA DEPTH for {inmagkey}')
		fits.setval(inim, inul5key, value=round(ul5, 3), comment=f'5 SIGMA DEPTH for {inmagkey}')

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
objlist = sorted(list(set([fits.getheader(inim)['OBJECT'] for inim in imlist])))
print('='*60)
print(f'{len(imlist)} IMAGE(s)')
print('-'*60)
for inim in imlist: print(os.path.basename(inim))
print('-'*60)
print(f'{len(objlist)} OBJECT(s)', end=' ')
for obj in objlist: print(obj, end=',')
print()
print('-'*60)

if len(objlist) > 1:
	preimlist = sorted(
		[glob.glob(f"{os.path.dirname(inim)}/{os.path.splitext(os.path.basename(inim))[0][0]}*{os.path.splitext(os.path.basename(inim))[0][-1]}{os.path.splitext(inim)[1]}")[0] for obj in objlist]
		)
else:
	preimlist = [imlist[0]]
#	Pre-photometry
##	Prohibit the crush event in querying reference catalog
mp_phot(preimlist, ncore=ncore)
#	Remove overlapped image
for inim in preimlist:
	if inim in imlist:
		imlist.remove(inim)
#	Photometry
mp_phot(imlist, ncore=ncore)

'''
#	Singleprocess
for i, inim in enumerate(imlist):
	print(f"[{i+1}/{len(imlist)}] {os.path.basename(inim)}")
	routine(inim)

#	Multiprocess
if __name__ == '__main__':
	#	Fixed the number of cores (=4)
	with multiprocessing.Pool(processes=ncore) as pool:
		results = pool.starmap(
			routine,
			zip(
				imlist,
				)
			)
'''






#	Header sample
header_to_put = '''
AUTHOR  = 'Gregory S.H. Paek'  / PHOTOMETRY AUTHOR
PHOTIME = '2021-12-11'         / PHTOMETRY TIME [KR]

# SEEING  =                3.626 / SEEING [arcsec]
# PEEING  =                 5.01 / SEEING [pixel]

ELLIP   =                 0.14 / ELLIPTICITY 1-B/A [0-1]
ELONG   =                1.162 / ELONGATION A/B [1-]
APERPIX =                11.22 / APERTURE DIAMETER [pixel]
APER    =                11.22 / BEST APERTURE DIAMETER in SNR curve [pix]
NAPER   =                 2.24 / N = APERTURE/PEEING
SKYSIG  =                5.333 / SKY SIGMA VALUE
SKYVAL  =                23.17 / SKY MEDIAN VALUE
REFCAT  = 'APASS   '           / REFERENCE CATALOG NAME

MAGLOW  =                 12.0 / REF MAG RANGE, BRIGHT LIMIT
MAGUP   =                 20.0 / REF MAG RANGE, DIM LIMIT

AUTO    =                    0 / MAG_AUTO DIAMETER [pix]
STDNUM_0=                  103 / # OF STD STARS for MAG_AUTO
ZP_0    =               24.813 / ZERO POINT for MAG_AUTO
ZPER_0  =                0.033 / ZERO POINT ERROR for MAG_AUTO
UL3_0   =                    0 / 3 sigma limit mag for MAG_AUTO
UL5_0   =                    0 / 5 sigma limit mag for MAG_AUTO

STDNUM_1=                  105 / # OF STD STARS for MAG_APER
ZP_1    =               24.668 / ZERO POINT for MAG_APER
ZPER_1  =                0.038 / ZERO POINT ERROR for MAG_APER
UL3_1   =               19.164 / 3 sigma limit mag for MAG_APER
UL5_1   =               18.609 / 5 sigma limit mag for MAG_APER
APER_1  =                6.744 / BEST GAUSSIAN APERTURE DIAMETER [pix]

STDNUM_2=                  103 / # OF STD STARS for MAG_APER_1
ZP_2    =               24.292 / ZERO POINT for MAG_APER_1
ZPER_2  =                0.051 / ZERO POINT ERROR for MAG_APER_1
UL3_2   =               19.341 / 3 sigma limit mag for MAG_APER_1
UL5_2   =               18.786 / 5 sigma limit mag for MAG_APER_1
APER_2  =                10.02 / 2*SEEING APERTURE DIAMETER [pix]
STDNUM_3=                  106 / # OF STD STARS for MAG_APER_2
ZP_3    =               24.606 / ZERO POINT for MAG_APER_2
ZPER_3  =                0.041 / ZERO POINT ERROR for MAG_APER_2
UL3_3   =               19.224 / 3 sigma limit mag for MAG_APER_2
UL5_3   =                18.67 / 5 sigma limit mag for MAG_APER_2
APER_3  =                15.03 / 3*SEEING APERTURE DIAMETER [pix]
STDNUM_4=                   97 / # OF STD STARS for MAG_APER_3
ZP_4    =               24.775 / ZERO POINT for MAG_APER_3
ZPER_4  =                0.031 / ZERO POINT ERROR for MAG_APER_3
UL3_4   =               18.953 / 3 sigma limit mag for MAG_APER_3
UL5_4   =               18.399 / 5 sigma limit mag for MAG_APER_3
APER_4  =                4.144 / FIXED 3" APERTURE DIAMETER [pix]
STDNUM_5=                  113 / # OF STD STARS for MAG_APER_4
ZP_5    =               23.613 / ZERO POINT for MAG_APER_4
ZPER_5  =                0.082 / ZERO POINT ERROR for MAG_APER_4
UL3_5   =                19.19 / 3 sigma limit mag for MAG_APER_4
UL5_5   =               18.636 / 5 sigma limit mag for MAG_APER_4
APER_5  =                6.906 / FIXED 5" APERTURE DIAMETER [pix]
STDNUM_6=                  103 / # OF STD STARS for MAG_APER_5
ZP_6    =               24.316 / ZERO POINT for MAG_APER_5
ZPER_6  =                0.051 / ZERO POINT ERROR for MAG_APER_5
UL3_6   =               19.339 / 3 sigma limit mag for MAG_APER_5
UL5_6   =               18.784 / 5 sigma limit mag for MAG_APER_5
'''