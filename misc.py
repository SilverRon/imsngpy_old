import requests
import os

def slack_bot(token, channel, text):
	response = requests.post("https://slack.com/api/chat.postMessage",
		headers={"Authorization": "Bearer "+token},
		data={"channel": channel,"text": text}
	)
	print(response)

def add_prepix(imlist, prefix):
	'''
	Add prefix on filename provided by list
	e.g.
	input : image list, prefix
	'''
	return [f'{os.path.dirname(inim)}/{prefix}{os.path.basename(inim)}' for inim in imlist]

def SE_seeing(inim, gain, pixscale, path_config, seeing_assume, frac=0.68, clean=True):
	'''
	inim = '/data3/paek/factory/loao/2020_1215/afzobj.NGC2207.20201216.0211.fits'
	path_config = '/home/paek/config'
	obs = 'LOAO'
	path_obs = '/home/paek/table/obs.dat'
	seeing_assume = 3 * u.arcsecond
	frac = 0.68
	'''
	#------------------------------------------------------------
	#	Input
	#------------------------------------------------------------
	hdr = fits.getheader(inim)
	a = hdr['naxis1']/2.
	b = hdr['naxis2']/2.
	#------------------------------------------------------------
	#	CCD information
	obsdict = getccdinfo(obs, path_obs)
	gain = obsdict['gain']
	pixscale = obsdict['pixelscale']
	fov = obsdict['fov']
	# rdnoise = obsdict['readoutnoise']
	#------------------------------------------------------------
	#	OUTPUT NAMES
	fmt0 = '.fits'
	fmt1 = '.fit'
	fmt2 = '.fts'

	if fmt0 in inim:
		cat = '{}/{}'.format(os.path.dirname(inim), os.path.basename(inim).replace(fmt0, '.cat'))
	elif fmt1 in inim:
		cat = '{}/{}'.format(os.path.dirname(inim), os.path.basename(inim).replace(fmt1, '.cat'))
	elif fmt2 in inim:
		cat = '{}/{}'.format(os.path.dirname(inim), os.path.basename(inim).replace(fmt2, '.cat'))

	# cat = '{}/{}'.format(os.path.dirname(inim), os.path.basename(inim).replace('.fits', '.cat'))
	#	SE configurations

	

	param = '{}/simple.param'.format(path_config)
	conv = '{}/simple.conv'.format(path_config)
	nnw = '{}/simple.nnw'.format(path_config)
	conf = '{}/simple.sex'.format(path_config)
	#	SE parameters
	param_insex = dict(
						#------------------------------
						#	CATALOG
						#------------------------------
						CATALOG_NAME = cat,
						#------------------------------
						#	CONFIG FILES
						#------------------------------
						CONF_NAME = conf,
						PARAMETERS_NAME = param,
						FILTER_NAME = conv,    
						STARNNW_NAME = nnw,
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

	com = phot.sexcom(inim, param_insex)
	os.system(com)
	rawtbl = ascii.read(cat)
	#	Point source selection
	indx_sel = np.where(
						(rawtbl['FLAGS'] == 0) &
						(sqsum((rawtbl['X_IMAGE']-a)/a, (rawtbl['Y_IMAGE']-b)/b) < frac) &
						(rawtbl['CLASS_STAR']>0.9) &
						(rawtbl['FWHM_WORLD']>0.0)
						)
	seltbl = rawtbl[indx_sel]
	#	Seeing in arcsecond/pixel as median value
	seeing = np.median(seltbl['FWHM_WORLD'].to(u.arcsecond))
	peeing = np.median(seltbl['FWHM_IMAGE']) * u.pix
	#	Header update
	try:
		puthdr(inim, hdrkey='SEEING', hdrval=seeing.value, hdrcomment='SEEING [arcsec]')
		puthdr(inim, hdrkey='PEEING', hdrval=peeing.value, hdrcomment='PEEING [pix]')
	except:
		print('try/except: Too low stars to measure seeing. Use 3.0 arcsecond seeing.')
		puthdr(inim, hdrkey='SEEING', hdrval=3.0, hdrcomment='SEEING [arcsec]')
		puthdr(inim, hdrkey='PEEING', hdrval=(3.0*u.arcsecond*pixscale).value, hdrcomment='PEEING [pix]')
	#	Clean output catalog
	if clean == True:
		rmcom = 'rm {}'.format(cat)
		# print(rmcom)
		os.system(rmcom)
	else:
		pass
	return seeing, peeing