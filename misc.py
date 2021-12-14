import requests
import os
from astropy.io import fits
import sys
sys.path.append('/home/paek/imsngpy')
from phot import *
import numpy as np
from astropy.io import ascii
from astropy import units as u

def slack_bot(token, channel, text):
	response = requests.post("https://slack.com/api/chat.postMessage",
		headers={"Authorization": "Bearer "+token},
		data={"channel": channel,"text": text}
	)
	print(response)

def add_prefix(imlist, prefix):
	'''
	Add prefix on filename provided by list
	e.g.
	input : image list, prefix
	'''
	return [f'{os.path.dirname(inim)}/{prefix}{os.path.basename(inim)}' for inim in imlist]

def add_suffix(imlist, suffix):
	'''
	Add suffix on filename provided by list
	e.g.
	input : image list, suffix
	'''
	return [f'{os.path.splitext(inim)[0]}.{suffix}{os.path.splitext(inim)[1]}' for inim in imlist]


def get_seeing(inim, gain, pixscale, fov, path_conf, path_param, path_conv, path_nnw, seeing_assume=3, frac=0.68, n_min_src=5):
	'''
	inim = '/data3/paek/factory/loao/2020_1215/afzobj.NGC2207.20201216.0211.fits'
	path_config = '/home/paek/config'
	obs = 'LOAO'
	path_obs = '/home/paek/table/obs.dat'
	seeing_assume = 3 * u.arcsecond
	frac = 0.68
	n_min_src = 5 # Calc. seeing value if has more than 5 stars, or not --> Fix 3 arcsec
	'''
	#------------------------------------------------------------
	#	Input
	#------------------------------------------------------------
	hdr = fits.getheader(inim)
	a = hdr['naxis1']/2.
	b = hdr['naxis2']/2.
	#	Small squre based on frac
	a_ = a*np.sqrt(frac)
	b_ = b*np.sqrt(frac)
	
	#------------------------------------------------------------
	#	OUTPUT NAMES
	outcat = f'{os.path.splitext(inim)[0]}.simple.cat'
	#	SE parameters
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
						GAIN = str(gain.value),
						PIXEL_SCALE = str(pixscale.value),
						#------------------------------
						#	STAR/GALAXY SEPARATION
						#------------------------------
						SEEING_FWHM = str(seeing_assume.value),
						)

	com = sexcom(inim, param_insex)
	os.system(com)
	rawtbl = ascii.read(outcat)
	#	Test plots
	# plt.hist(rawtbl['FLAGS'], bins=np.arange(0, 10+1, 1))
	# plt.hist(rawtbl['CLASS_STAR'], bins=np.arange(0, 1.0+0.1, 0.1))
	# plt.hist(rawtbl['FWHM_WORLD'].arcsec, bins=np.arange(0.5, 10+0.1, 0.1))

	#	Point source selection
	indx_sel = np.where(
						#	Good point source
						(rawtbl['FLAGS'] == 0) &
						(rawtbl['CLASS_STAR']>0.9) &
						(rawtbl['FWHM_WORLD']>0.0) &
						#	Nearby center
						(rawtbl['X_IMAGE']>a-a_) & (rawtbl['X_IMAGE']<a+a_) &
						(rawtbl['Y_IMAGE']>b-b_) & (rawtbl['Y_IMAGE']<b+b_)
						#	Former criterion (version 1)
						# (sqsum((rawtbl['X_IMAGE']-a)/a, (rawtbl['Y_IMAGE']-b)/b) < frac) &
						)
	seltbl = rawtbl[indx_sel]
	#	Seeing in arcsecond/pixel as median value
	if len(seltbl) >= n_min_src:
		seeing = np.median(seltbl['FWHM_WORLD'].to(u.arcsecond))
		peeing = np.median(seltbl['FWHM_IMAGE']) * u.pix
	else:
		print(f'{os.path.basename(inim)} has not enough point sources (n<5). Set seeing {seeing_assume}')
		seeing = seeing_assume*u.arcsec
		peeing = seeing_assume*pixscale
	#	Test plots
	# plt.close()
	# plt.title(f'frac = {frac}')
	# plt.plot(rawtbl['X_IMAGE'], rawtbl['Y_IMAGE'], marker='o', ls='none', label=f'all({len(rawtbl)})')
	# plt.plot(seltbl['X_IMAGE'], seltbl['Y_IMAGE'], marker='.', ls='none', label=f'select({len(seltbl)}):seeing={round(seeing.value, 3)*u.arcsec}')
	# plt.xlabel('X_IMAGE')
	# plt.ylabel('Y_IMAGE')
	# plt.legend()
	#	
	#	Header update
	fits.setval(inim, 'NSEEING', value=len(seltbl), comment='# of point sources to measure the seeing')
	fits.setval(inim, 'SEEING', value=round(seeing.value, 3), comment='SEEING [arcsec]')
	fits.setval(inim, 'PEEING', value=round(peeing.value, 3), comment='PEEING [arcsec]')
	return seeing, peeing
#------------------------------------------------------------
def fnamechange(inim):
	with fits.open(inim) as hdul:
		hdul.verify('fix')
		hdr = hdul[0].header
		obs = hdr['OBSERVAT']
		obj = hdr['OBJECT']
		dateobs = hdr['DATE-OBS']
		datestr = dateobs[0:4]+dateobs[5:7]+dateobs[8:10]
		timestr = dateobs[11:13]+dateobs[14:16]+dateobs[17:19]
		filte = hdr['FILTER']
		exptime = int(hdr['EXPTIME'])

	newim = f'Calib-{obs}-{obj}-{datestr}-{timestr}-{filte}-{exptime}.fits'
	return newim
#------------------------------------------------------------
def identify_ccdinfo(ic0, obs, ccdtbl):
	#	Identify CCD
	print(f"""{'-'*60}\n#\tIDENTIFY CCD\n{'-'*60}""")
	for key, val, suf, ccd in zip((ccdtbl['key'][ccdtbl['obs']==obs]), (ccdtbl['value'][ccdtbl['obs']==obs]), (ccdtbl['suffix'][ccdtbl['obs']==obs]), (ccdtbl['ccd'][ccdtbl['obs']==obs])):
		if (key.lower() in ic0.keywords) & (val == ic0.summary[key.lower()][0]):
			ccdkey = key
			ccdval = val
			ccdtype = ccd
			if suf.mask == True:
				#	No suffix
				suffix = ''
				obsccd = f'{obs}'
			else:
				suffix = suf
				obsccd = f'{obs}_{suffix}'
			print(f'OBSERVAT : {obs}\nCCD KEYWORD : {key}\nCCD HEADER VALUE : {val}\nCCD NAME : {ccdtype}\nSUFFIX : {suffix}\n==> OBS_CCD : {obsccd}')
	#	CCD INFO
	indx_ccd = np.where(
		(ccdtbl['obs']==obs) &
		(ccdtbl['key']==ccdkey) &
		(ccdtbl['value']==ccdval)
	)
	return ccdkey, ccdval, ccdtype, obsccd