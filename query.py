#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#============================================================
#	https://ps1images.stsci.edu/ps1image.html
#	By Sophia Kim 2019.01.22. based on code PS1 suggests on the link above
#	Pan-STARRS DR1 data query
#	from https://michaelmommert.wordpress.com/2017/02/13/accessing-the-gaia-and-pan-starrs-catalogs-using-python/
#	By CS Choi 
#	REVISED AND ORGANIZED BY GREGORY S.H. PAEK
#	UPDATE : 20.01.03
#============================================================
from __future__ import print_function
from astropy.io import fits
from astropy.io import ascii
import glob
import sys
sys.path.append('/home/paek/imsngpy')
from misc import *
import numpy as np
from astroquery.vizier import Vizier 
# from astroquery.vizier import *
from astropy.coordinates import Angle
import astropy.units as u
from astropy.table import Table
#============================================================
def geturl(ra, dec, size=240, output_size=None, filters="grizy", format="jpg", color=False):
	"""
	Get URL for images in the table
	
	ra, dec = position in degrees
	size = extracted image size in pixels (0.25 arcsec/pixel)
	output_size = output (display) image size in pixels (default = size).
				  output_size has no effect for fits format images.
	filters = string with filters to include
	format = data format (options are "jpg", "png" or "fits")
	color = if True, creates a color image (only for jpg or png format).
			Default is return a list of URLs for single-filter grayscale images.
	Returns a string with the URL
	"""	
	#------------------------------------------------------------
	if color and format == "fits":
		raise ValueError("color images are available only for jpg or png formats")
	if format not in ("jpg","png","fits"):
		raise ValueError("format must be one of jpg, png, fits")

	table	= getimages(ra, dec, size=size, filters=filters)
	url		= (	"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
				"ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
	if output_size:
		url = url + "&output_size={}".format(output_size)
	# sort filters from red to blue
	flist = ["yzirg".find(x) for x in table['filter']]
	table = table[np.argsort(flist)]
	if color:
		if len(table) > 3:
			# pick 3 filters
			table = table[[0,len(table)//2,len(table)-1]]
		for i, param in enumerate(["red","green","blue"]):
			url = url + "&{}={}".format(param,table['filename'][i])
	else:
		urlbase = url + "&red="
		url = []
		for filename in table['filename']:
			url.append(urlbase+filename)
	return url
#------------------------------------------------------------
def getps1image(ra,dec,size=240,filters="grizy"):
	"""
	Query ps1filenames.py service to get a list of images
	ra, dec = position in degrees
	size = image size in pixels (0.25 arcsec/pixel)
	filters = string with filters to include
	Returns a table with the results
	#------------------------------------------------------------
	#	SAMPLE LINES FOR DOWNLOADING SINGLE IMAGE FROM PS1

	param_down	= dict(  name = 'AT2019dko',
						ra = 181.4641667,
						dec = 67.2569444,
						size = 5000,
						output_size = None,
						filters = 'r',
						format = 'fits',
						save_dir='.')
	try:
		query.downimage_routine(**param_down)
	except:
		try:
			query.downimage_routine(**param_down)
		except:
			try:
				query.downimage_routine(**param_down)		
			except:
				query.downimage_routine(**param_down)

	#------------------------------------------------------------
	#	SAMPLE LINES FOR DOWNLOADING MULTIPLE IMAGES FROM PS1
	#	intbl : 'ascii' TABLE
	for i in range(len(intbl)):
		print('['+str(i+1)+'/'+str(len(intbl))+']')
		param_down= dict(	name		= intbl['name'][i],
							ra			= intbl['ra'][i],
							dec			= intbl['dec'][i],
							size		= 5000,
							output_size	= None,
							filters		= 'r',
							format		= "fits")
		try:
			downimage_routine(**param_down)
		except:
			try:
				downimage_routine(**param_down)
			except:
				downimage_routine(**param_down)

	"""	
	service	= "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
	url		= ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
		   "&filters={filters}").format(**locals())
	table	= Table.read(url, format='ascii.tab')
	return table
#------------------------------------------------------------
def downimage_routine(outim, name, ra, dec, size, output_size, filters, format, save_dir='.'):
	filt = filters
	param_geturl= dict(	ra			= ra,
						dec			= dec,
						size		= size,
						output_size	= None,
						filters		= filt,
						format		= "fits")
	url			= geturl(**param_geturl)
	fh		= fits.open(url[0])
	# newname	= 'Ref'+'-PS1-'+name+'-'+filt+'.fits'
	newname = outim
	fh.writeto(save_dir+'/'+newname, overwrite=True)
	pan, panhd	= fits.getdata(save_dir+'/'+newname, header=True)
	pan0		= np.nan_to_num(pan)
	fits.writeto(save_dir+'/'+newname, pan0, panhd, overwrite=True)
#------------------------------------------------------------
def querybox(refcatname, racent, decent, outcat, radius=0.5, refmagkey=''):
	'''
	#------------------------------------------------------------
	#	REF. CATALOG QUERY
	#------------------------------------------------------------
	refcatname = 'PS1'
	reftbl = querybox(**param_query)
	outcat = '/data3/paek/factory/refcat/test.cat'
	'''
	#	PS1
	if refcatname.upper() == 'PS1':
		#	Query
		querytbl = query_ps1(racent, decent, radius=radius)
		#	Convert magnitude system to unified system
		reftbl = convert_ps1toTonry(querytbl)
		# reftbl.write(outcat, format='ascii.tab', overwrite=True)
		reftbl.write(outcat, overwrite=True)
	#	SDSS
	elif refcatname.upper()	== 'SDSS':
		#	Query
		querytbl = query_sdss(racent, decent, radius=radius)
		#	Convert magnitude system to unified system
		reftbl = convert_sdsstoBlaton(querytbl)
		# reftbl.write(outcat, format='ascii.tab', overwrite=True)
		reftbl.write(outcat, overwrite=True)
	#	APASS
	elif refcatname.upper() == 'APASS':
		#	Query
		querytbl = query_apass(racent, decent, radius=radius)
		#	Convert magnitude system to unified system
		reftbl = convert_apasstoBlatonAndLupton(querytbl)
		# reftbl.write(outcat, format='ascii.tab', overwrite=True)
		reftbl.write(outcat, overwrite=True)
	#	2MASS for JHK-bands
	elif refcatname.upper()	== '2MASS':
		#	Query
		querytbl = query_2mass(racent, decent, radius=radius)
		#	Convert magnitude system to unified system
		reftbl = convert_2masstoAB(querytbl)
		# reftbl.write(outcat, format='ascii.tab', overwrite=True)
		reftbl.write(outcat, overwrite=True)
	return reftbl
#-------------------------------------------------------------------------#
def query_sdss(radeg, dedeg, radius=1.0):
	"""
	SDSS QUERY
	INPUT   :   NAME, RA [deg], Dec [deg], radius
	OUTPUT  :   QUERY TABLE
				sdss-[NAME].cat
	"""
	from astroquery.vizier import Vizier 
	queryname = 'SDSS'
	comment = f"""{'='*60}
Query {queryname} catalog from Vizier...
{'-'*60}
RA center    [deg]: {radeg}
Dec center   [deg]: {dedeg}
Query radius [deg]: {radius}"""
	print(comment)

	#	QUERY PART
	from astropy.coordinates import SkyCoord
	Vizier.ROW_LIMIT = -1
	import time
	import astropy.units as u
	st = time.time()
	query = Vizier.query_region(
		SkyCoord(ra=radeg, dec=dedeg, unit=(u.deg, u.deg), frame='icrs'),
		width=f'{radius*60}m', catalog=["SDSS12"]
		)
	catalog = query.keys()[0]
	print(f"{'-'*60}\nGet {catalog} ({round(time.time()-st, 1)} sec)")
	rawtbl = query[catalog]
	indx_ps = np.where(
		(rawtbl['q_mode']=='+') &
		(rawtbl['class']==6) &
		((rawtbl['Q']==2) | (rawtbl['Q']==3))
	)
	querytbl = rawtbl[indx_ps]
	print(f'Filter good point sources : {len(rawtbl)} --> {len(querytbl)} ({round(len(querytbl)*1e2/len(rawtbl), 1)}%)')

	return querytbl
#-------------------------------------------------------------------------#
def query_ps1(radeg, dedeg, radius=1.0):
	"""
	Filter point sources to measure zeropoint
	#	SELECT STARS FROM STARS & GALAXIES (iPSF - iKron <= 0.05)
	https://outerspace.stsci.edu/display/PANSTARRS/How+to+separate+stars+and+galaxies#
	outname = 'ps1-test.cat'
	"""
	from astroquery.vizier import Vizier
	queryname = 'PS1'
	comment = f"""{'='*60}
Query {queryname} catalog from Vizier...
{'-'*60}
RA center    [deg]: {radeg}
Dec center   [deg]: {dedeg}
Query radius [deg]: {radius}"""
	print(comment)
	#	QUERY PART
	from astropy.coordinates import SkyCoord
	Vizier.ROW_LIMIT = -1
	import time
	import astropy.units as u
	st = time.time()
	query = Vizier.query_region(
		SkyCoord(ra=radeg, dec=dedeg, unit=(u.deg, u.deg), frame='icrs'),
		width=f'{radius*60}m',
		catalog=["II/349/ps1"]
		)
	catalog = query.keys()[0]
	print(f"{'-'*60}\nGet {catalog} ({round(time.time()-st, 1)} sec)")
	#	Raw Table
	rawtbl = query[0]

	#	Flag integer --> binary
	import numpy as np
	f_objID_bin = np.array([bin(i)[2:] for i in rawtbl['f_objID']])

	#	SELECT POINT SOURCE & NON-VARIABLE & GOOD QUALITY STARS
	flagcol = [f'{n}' for n in np.flip(np.arange(1, 31+1, 1))]
	data_rows = []
	for i in range(len(rawtbl)):
		row = list(f_objID_bin[i])
		while len(row) < 31:
			row.insert(0, 0)
		data_rows.append(row)
	from astropy.table import Table
	flagtbl = Table(
		rows=data_rows,
		names=flagcol,
		dtype=['i']*len(flagcol),
	)

	indx_ps = np.where(
		(flagtbl['23']!=1) &
		(flagtbl['24']!=1) &
		(flagtbl['2']!=1) &
		(flagtbl['3']!=1) &
		(flagtbl['4']!=1) &
		(flagtbl['5']!=1) &
		(flagtbl['6']!=1) &
		(flagtbl['7']!=1) &
		(flagtbl['8']!=1) &
		#	Point source
		(rawtbl['imag']-rawtbl['iKmag']<=0.05) &
		(rawtbl['Qual']<128) &
		#	Filter saturated sources
		(rawtbl['gmag']>12) &
		(rawtbl['rmag']>12) &
		(rawtbl['imag']>12) &
		(rawtbl['zmag']>12) &
		(rawtbl['ymag']>12) 
	)

	#	SELECT STARS FROM STARS & GALAXIES (iPSF - iKron <= 0.05)
	querytbl = rawtbl[indx_ps]

	#	Plot
	'''
	plt.plot(rawtbl['imag'], rawtbl['imag']-rawtbl['iKmag'], ls='', marker='.', c='k', alpha=0.125)
	plt.plot(querytbl['imag'], querytbl['imag']-querytbl['iKmag'], ls='', marker='.', c='tomato', alpha=0.25)
	plt.axhline(y=0.05, ls='-', c='red')
	'''

	#	CHANGE TO GENERTAL COL. NAMES
	querytbl.rename_columns(
		names=[
			'objID',
			'RAJ2000',
			'DEJ2000',
			'Qual',
			'Nd',
			'Ns',
		],
		new_names=[
			'NUMBER',
			'RA_ICRS',
			'DE_ICRS',
			'Q',
			'N_obs',
			'N_img'
		],
	)

	# querytbl['gmag'].name = 'gmag'
	# querytbl['e_gmag'].name = 'e_gmag'
	# querytbl['rmag'].name = 'rmag'
	# querytbl['e_rmag'].name = 'e_rmag'
	# querytbl['imag'].name = 'imag'
	# querytbl['e_imag'].name = 'e_imag'
	# querytbl['zmag'].name = 'zmag'
	# querytbl['e_zmag'].name = 'e_zmag'
	# querytbl['ymag'].name = 'ymag'
	# querytbl['e_ymag'].name = 'e_ymag'

	# querytbl.write(outname, format='ascii.tab', overwrite=True)
	print(f'Filter good point sources : {len(rawtbl)} --> {len(querytbl)} ({round(len(querytbl)*1e2/len(rawtbl), 1)}%)')
	return querytbl
#-------------------------------------------------------------------------#
def query_apass(radeg, dedeg, radius=1.0):
	"""
	APASS QUERY
	INPUT   :   NAME, RA [deg], Dec [deg], radius
	OUTPUT  :   QUERY TABLE
				apass-[NAME].cat
	"""

	from astroquery.vizier import Vizier 

	queryname = 'APASS'
	comment = f"""{'='*60}
	Query {queryname} catalog to from Vizier...
{'-'*60}
RA center    [deg]: {radeg}
Dec center   [deg]: {dedeg}
Query radius [deg]: {radius}"""
	print(comment)

	#	QUERY PART
	from astropy.coordinates import SkyCoord
	Vizier.ROW_LIMIT = -1
	import time
	st = time.time()
	query = Vizier.query_region(
		SkyCoord(ra=radeg, dec=dedeg, unit=(u.deg, u.deg), frame='icrs'),
		width=f'{radius*60}m', catalog=["APASS9"]
		)
	catalog = query.keys()[0]
	print(f"{'-'*60}\nGet {catalog} ({round(time.time()-st, 1)} sec)")
	rawtbl = query[catalog]

	#   Vega    : B, V
	#   AB      : g, r, i
	#   Vega - AB Magnitude Conversion (Blanton+07)
	#   U       : m_AB - m_Vega = 0.79
	#   B       : m_AB - m_Vega =-0.09
	#   V       : m_AB - m_Vega = 0.02
	#   R       : m_AB - m_Vega = 0.21
	#   I       : m_AB - m_Vega = 0.45
	#   J       : m_AB - m_Vega = 0.91
	#   H       : m_AB - m_Vega = 1.39
	#   K       : m_AB - m_Vega = 1.85

	querytbl = rawtbl
	querytbl['recno'].name = 'NUMBER'
	querytbl['RAJ2000'].name = 'RA_ICRS'
	querytbl['DEJ2000'].name = 'DE_ICRS'
	querytbl['nobs'].name = 'Numb_obs'
	querytbl['mobs'].name = 'Numb_img'

	querytbl['g_mag']    .name = 'gmag'
	querytbl['e_g_mag']  .name = 'e_gmag'
	querytbl['r_mag']    .name = 'rmag'
	querytbl['e_r_mag']  .name = 'e_rmag'
	querytbl['i_mag']    .name = 'imag'
	querytbl['e_i_mag']  .name = 'e_imag'

	# querycat['B-V']     = dum['B-V']    + (-0.09 - 0.02)
	# querycat['e_B-V']   = dum['e_B-V']
	# querycat['Bmag']    = dum['Bmag']   - 0.09  # [Vega] to [AB]
	# querycat['e_Bmag']  = dum['e_Bmag']
	# querycat['Vmag']    = dum['Vmag']   + 0.02  # [Vega] to [AB]
	# querycat['e_Vmag']  = dum['e_Vmag']
	print(f'No filtering : {len(querytbl)}')
	return querytbl
#-------------------------------------------------------------------------#
def query_2mass(radeg, dedeg, radius=1.0):
	"""
	QUERY Point Source Catalog(PSC) PROVIDED BY 2MASS
	REMOVE COMTAMINATED SOURCE BY EXTENDED SOURCE AND MINOR PLANET
	IF GIVE BAND INPUT, 

	INPUT	:	NAME, RA [deg], DEC [deg], BAND, RADIUS
	OUTPUT	:	TABLE, MAGNITUDE [AB]
	"""
	from astroquery.vizier import Vizier 

	queryname = '2MASS'
	comment = f"""{'='*60}
	Query {queryname} catalog to from Vizier...
{'-'*60}
RA center    [deg]: {radeg}
Dec center   [deg]: {dedeg}
Query radius [deg]: {radius}"""
	print(comment)

	#	QUERY PART
	from astropy.coordinates import SkyCoord
	Vizier.ROW_LIMIT = -1
	import time
	st = time.time()
	query = Vizier.query_region(
		SkyCoord(ra=radeg, dec=dedeg, unit=(u.deg, u.deg), frame='icrs'),
		width=f'{radius*60}m', catalog=["II/246"]
		)
	catalog = query.keys()[0]
	print(f"{'-'*60}\nGet {catalog} ({round(time.time()-st, 1)} sec)")
	rawtbl = query[catalog]

	for j in ['Q', 'R', 'C']:
		rawtbl[f'J_{j}flg'] = ['']*len(rawtbl)
		rawtbl[f'H_{j}flg'] = ['']*len(rawtbl)
		rawtbl[f'K_{j}flg'] = ['']*len(rawtbl)
		for i, q in enumerate(rawtbl[f'{j}flg']):
			rawtbl[f'J_{j}flg'][i] = str(q[0])
			rawtbl[f'H_{j}flg'][i] = str(q[1])
			rawtbl[f'K_{j}flg'][i] = str(q[2])
	#	http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=II/246/out
	#	SELECT POINT SOURCE & NON-VARIABLE & GOOD QUALITY STARS
	indx_ps = np.where(
		#	Asteroid
		(rawtbl['Aflg']==0) &
		#	Contamination by extended source
		(rawtbl['Xflg']==0) &
		#	Artifcat contamination & confusion
		(rawtbl['J_Cflg']=='0') &
		(rawtbl['H_Cflg']=='0') &
		(rawtbl['K_Cflg']=='0') &
		#	Quality flag
		#	J
		(
			(rawtbl['J_Qflg']=='A') | 
			(rawtbl['J_Qflg']=='B') | 
			(rawtbl['J_Qflg']=='C') | 
			(rawtbl['J_Qflg']=='D')
			) &
		#	H
		(
			(rawtbl['H_Qflg']=='A') | 
			(rawtbl['H_Qflg']=='B') | 
			(rawtbl['H_Qflg']=='C') | 
			(rawtbl['H_Qflg']=='D')
			) &
		#	K
		(
			(rawtbl['K_Qflg']=='A') | 
			(rawtbl['K_Qflg']=='B') | 
			(rawtbl['K_Qflg']=='C') | 
			(rawtbl['K_Qflg']=='D')
			) &
		#	Read flag
		#	J
		(
			(rawtbl['J_Rflg']=='1') | 
			(rawtbl['J_Rflg']=='2') | 
			(rawtbl['J_Rflg']=='3')
			) &
		#	H
		(
			(rawtbl['H_Rflg']=='1') | 
			(rawtbl['H_Rflg']=='2') | 
			(rawtbl['H_Rflg']=='3')
			) &
		#	K
		(
			(rawtbl['K_Rflg']=='1') | 
			(rawtbl['K_Rflg']=='2') | 
			(rawtbl['K_Rflg']=='3') 
			)
	)
	#	CHANGE TO GENERTAL COL. NAMES
	querytbl = rawtbl[indx_ps]
	querytbl['_2MASS'].name ='name'
	querytbl['RAJ2000'].name ='ra'
	querytbl['DEJ2000'].name ='dec'
	#	AB OFFSET
	querytbl['Jmag'].name ='J'	
	querytbl['e_Jmag'].name ='Jerr'
	querytbl['Hmag'].name ='H'
	querytbl['e_Hmag'].name ='Herr'
	querytbl['Kmag'].name ='K'
	querytbl['e_Kmag'].name ='Kerr'
	# querytbl['Qflg'].name ='Qflg'
	# querytbl['Rflg'].name ='Rflg'
	# querytbl['Bflg'].name ='Bflg'
	# querytbl['Cflg'].name ='Cflg'
	print(f'Filter good point sources : {len(rawtbl)} --> {len(querytbl)} ({round(len(querytbl)*1e2/len(rawtbl), 1)}%)')
	return querytbl
#-------------------------------------------------------------------------#
def convert_2masstoAB(intbl):
	intbl['Jmag'] = intbl['J']+0.91
	intbl['Hmag'] = intbl['H']+1.39
	intbl['Kmag'] = intbl['K']+1.85
	return intbl
#-------------------------------------------------------------------------#
def convert_sdsstoBlaton(intbl):
	"""
	=====================================================================
	CONVERSION SDSS FILTER SYSTEM TO JOHNSON FILTER SYSTEM [AB]
	INPUT   :   QUERIED SDSS CATALOG
	OUTPUT  :   1.  ONLY STAR CLASS (=6)
				2.  (2)ACCEPTABLE & (3)GOOD QUALITY
				3.  NAME, RA, Dec, ugriz, BVRI [AB] (Blaton+07), clean
				sdss-conv.cat
	---------------------------------------------------------------------
	Blaton+07
	CONVERSION TABLE to AB
	---------------------------------------------------------------------
	Equation                                        Color Dispersion
	---------------------------------------------------------------------
	U   = u - 0.0682 - 0.0140[ (u - g) - 1.2638 ]   sigma[u - g] = 0.26
	B   = u - 1.0286 - 0.7981[ (u - g) - 1.2638 ]   sigma[u - g] = 0.26
	*B   = g + 0.2354 + 0.3915[ (g - r) - 0.6102 ]   sigma[g - r] = 0.15
	V   = g - 0.3516 - 0.7585[ (g - r) - 0.6102 ]   sigma[g - r] = 0.15
	R   = r - 0.0576 - 0.3718[ (r - i) - 0.2589 ]   sigma[r - i] = 0.10
	I   = i - 0.0647 - 0.7177[ (i - z) - 0.2083 ]   sigma[i - z] = 0.10
	=====================================================================
	"""	
	# name    = intbl['SDSS12']
	# ra, de  = intbl['RA_ICRS'],    intbl['DE_ICRS']
	# clean   = intbl['q_mode']
	# clas    = intbl['class']
	# Q       = intbl['Q']

	u, uer  = intbl['umag'], intbl['e_umag']
	g, ger  = intbl['gmag'], intbl['e_gmag']
	r, rer  = intbl['rmag'], intbl['e_rmag']
	i, ier  = intbl['imag'], intbl['e_imag']
	z, zer  = intbl['zmag'], intbl['e_zmag']

	uger, grer, rier, izer	= 0.26, 0.15, 0.10, 0.10
	ug, gr, ri, iz			= u-g, g-r, r-i, i-z
	'''
	U		= u - 0.0682 - 0.0140*(ug-1.2638)
	Uer		= np.sqrt( ((uer)**2.) + ((-0.0140*uger)**2.) )
	'''
	B		= g + 0.2354 + 0.3915*(gr-0.6102)
	Ber		= np.sqrt( ((ger)**2.) + ((0.3915*grer)**2.) )
	V		= g - 0.3516 - 0.7585*(gr-0.6102)
	Ver		= np.sqrt( ((ger)**2.) + ((-0.7585*grer)**2.) )
	R		= r - 0.0576 - 0.3718*(ri-0.2589)
	Rer		= np.sqrt( ((rer)**2.) + ((-0.3718*rier)**2.) )
	I		= i - 0.0647 - 0.7177*(iz-0.2083)
	Ier		= np.sqrt( ((ier)**2.) + ((-0.7177*izer)**2.) )

	intbl['RA_ICRS'].name='ra'
	intbl['DE_ICRS'].name='dec'
	# intbl['mode'].name=
	# intbl['q_mode'].name=
	intbl['class'].name='class'
	intbl['SDSS12'].name='name'
	# intbl['m_SDSS12'].name=
	# intbl['ObsDate'].name=
	# intbl['Q'].name=
	intbl['umag'].name='u'
	intbl['e_umag'].name='uerr'
	intbl['gmag'].name='g'
	intbl['e_gmag'].name='gerr'
	intbl['rmag'].name='r'
	intbl['e_rmag'].name='rerr'
	intbl['imag'].name='i'
	intbl['e_imag'].name='ierr'
	intbl['zmag'].name='z'
	intbl['e_zmag'].name='zerr'
	# intbl['zsp'].name=
	# intbl['zph'].name=
	# intbl['e_zph'].name=
	# intbl['__zph_'].name=

	intbl['B'] = B
	intbl['Berr'] = Ber
	intbl['V'] = V
	intbl['Verr'] = Ver
	intbl['R'] = R
	intbl['Rerr'] = Rer
	intbl['I'] = I
	intbl['Ierr'] = Ier
	'''
	outbl	= Table(
		[name, ra, de, u, uer, g, ger, r, rer, i, ier, z, zer, B, Ber, V, Ver, R, Rer, I, Ier, clean, Q, clas],
		names=['name', 'ra', 'dec', 'u', 'uerr', 'g', 'gerr', 'r', 'rerr', 'i', 'ierr', 'z', 'zerr', 'B', 'Berr', 'V', 'Verr', 'R', 'Rerr', 'I', 'Ierr', 'clean', 'Q', 'class']
		)
	'''
	return intbl
#-------------------------------------------------------------------------#
def convert_apasstoBlatonAndLupton(intbl):
	"""
	=====================================================================
	CONVERSION SDSS FILTER SYSTEM TO JOHNSON FILTER SYSTEM [AB]
	INPUT   :   QUERIED SDSS CATALOG
	OUTPUT  :   1.  ONLY STAR CLASS (=6)
				2.  (2)ACCEPTABLE & (3)GOOD QUALITY
				3.  NAME, RA, Dec, ugriz, BVRI [AB] (Blaton+07), clean
				sdss-conv.cat
	---------------------------------------------------------------------
	Blaton+07
	CONVERSION TABLE to AB
	---------------------------------------------------------------------
	Equation                                        Color Dispersion
	---------------------------------------------------------------------
	U   = u - 0.0682 - 0.0140[ (u - g) - 1.2638 ]   sigma[u - g] = 0.26
	B   = u - 1.0286 - 0.7981[ (u - g) - 1.2638 ]   sigma[u - g] = 0.26
	*B   = g + 0.2354 + 0.3915[ (g - r) - 0.6102 ]   sigma[g - r] = 0.15
	V   = g - 0.3516 - 0.7585[ (g - r) - 0.6102 ]   sigma[g - r] = 0.15
	R   = r - 0.0576 - 0.3718[ (r - i) - 0.2589 ]   sigma[r - i] = 0.10
	I   = i - 0.0647 - 0.7177[ (i - z) - 0.2083 ]   sigma[i - z] = 0.10
	=====================================================================
	"""

	reftbl	= intbl
	#	Vega --> AB
	reftbl['B-V'] = reftbl['B-V']+(-0.09 - 0.02)
	reftbl['Bmag'] = reftbl['Bmag']-0.09
	reftbl['Vmag'] = reftbl['Vmag']+0.02
	
	name    = reftbl['NUMBER']
	ra, de  = reftbl['RA_ICRS'],    reftbl['DE_ICRS']
	Numb_obs= reftbl['Numb_obs']
	Numb_img= reftbl['Numb_img']
	B		= reftbl['Bmag']
	Ber		= reftbl['e_Bmag']
	V		= reftbl['Vmag']
	Ver		= reftbl['e_Vmag']	
	BV		= reftbl['B-V']
	e_BV	= reftbl['e_B-V']

	g, ger  = reftbl['gmag'],       reftbl['e_gmag']
	r, rer  = reftbl['rmag'],       reftbl['e_rmag']
	i, ier  = reftbl['imag'],       reftbl['e_imag']

	grer, rier		= 0.15, 0.10
	gr, ri			= g-r, r-i

	R		= r - 0.0576 - 0.3718*(ri-0.2589)
	Rer		= np.sqrt( ((rer)**2.) + ((-0.3718*rier)**2.) )
	'''
	#	Blaton
	I		= i - 0.0647 - 0.7177*(iz-0.2083)
	Ier		= np.sqrt( ((ier)**2.) + ((-0.7177*izer)**2.) )
	'''
	#	Lupton
	Isig  = 0.0078
	I       = r - 1.2444*(r - i) - 0.3820
	Ier1    = ((1.-1.2444)**2)*(ger**2)+((+1.2444)**2)*(rer**2)
	Ier     = np.sqrt(Ier1**2 + Isig**2)
	#	Vega to AB
	I = I + 0.45

	outbl	= Table(
		[name, ra, de, Numb_obs, Numb_img, g, ger, r, rer, i, ier, B, Ber, V, Ver, R, Rer, I, Ier],
		names=['name', 'ra', 'dec', 'numb_obs', 'numb_img', 'g', 'gerr', 'r', 'rerr', 'i', 'ierr', 'B', 'Berr', 'V', 'Verr', 'R', 'Rerr', 'I', 'Ierr']
		)
	return outbl
#-------------------------------------------------------------------------#
def convert_ps1toTonry(intbl):
	'''
	PS1 -> Johnson/COusins [Vega] -> [AB]	(Tonry+12)
	#   Vega - AB Magnitude Conversion (Blanton+07)
	U       : m_AB - m_Vega = 0.79
	B       : m_AB - m_Vega =-0.09
	V       : m_AB - m_Vega = 0.02
	R       : m_AB - m_Vega = 0.21
	I       : m_AB - m_Vega = 0.45
	J       : m_AB - m_Vega = 0.91
	H       : m_AB - m_Vega = 1.39
	K       : m_AB - m_Vega = 1.85
	'''
	#	Table component to variable
	name    = intbl['NUMBER']
	# name    = intbl['objID']
	Q		= intbl['Q']
	ra, dec = intbl['RA_ICRS'], intbl['DE_ICRS']
	g, ger  = intbl['gmag'],	intbl['e_gmag']
	r, rer  = intbl['rmag'],	intbl['e_rmag']
	i, ier  = intbl['imag'],	intbl['e_imag']
	z, zer  = intbl['zmag'],	intbl['e_zmag']
	y, yer  = intbl['ymag'],	intbl['e_ymag']

	#	TRANSF. ERROR FOR B CONST. TERMS
	Bsig, Vsig, Rsig, Isig	= 0.034, 0.012, 0.01, 0.016
	#	COLOR TERM
	gr		= intbl['gmag']-intbl['rmag']
	grer	= sqsum(intbl['e_gmag'], intbl['e_rmag'])
	#	CONVERT TO B
	B0		= 0.213
	B1		= 0.587
	B		= B0 + B1*gr + intbl['gmag'] - 0.09
	Ber		= sqsum( Bsig, sqsum(B1*grer, intbl['e_gmag']) )
	#	CONVERT TO V
	B0		= 0.006
	B1		= 0.474
	V		= B0 + B1*gr + intbl['rmag'] + 0.02
	Ver	= sqsum( Bsig, sqsum(B1*grer, intbl['e_rmag']) )
	#	CONVERT TO R
	B0		=-0.138
	B1		=-0.131
	R		= B0 + B1*gr + intbl['rmag'] + 0.21
	Rer		= sqsum( Rsig, sqsum(B1*grer, intbl['e_rmag']) )
	#	CONVERT TO I
	B0		=-0.367
	B1		=-0.149
	I		= B0 + B1*gr + intbl['imag'] + 0.45
	Ier		= sqsum( Isig, sqsum(B1*grer, intbl['e_imag']) )
	#	REULST
	outbl = Table(
		[name, ra, dec, g, ger, r, rer, i, ier, z, zer, y, yer, B, Ber, V, Ver, R, Rer, I, Ier, Q],
		names=['name', 'ra', 'dec', 'g', 'gerr', 'r', 'rerr', 'i', 'ierr', 'z', 'zerr', 'y', 'yerr', 'B', 'Berr', 'V', 'Verr', 'R', 'Rerr', 'I', 'Ierr', 'Q']
		)
	return outbl
#-------------------------------------------------------------------------#
#	Move to misc?
