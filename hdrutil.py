from astropy.io import fits



def correcthdr_routine(path_data, hdrtbl, obs):
	'''
	1.	NGC337		-> NGC0337
		ngc0337		-> NGC0337
	'''
	from astropy.io import fits
	from astropy.time import Time
	comment = '-'*60+'\n' \
			+ 'CHANGE TO UNIVERSAL HEADER SERIES ...\n' \
			+ '-'*60+'\n'
	print(comment)
	objfilterlist = []
	objexptimelist = []
	flatfilterlist = []
	darkexptimelist = []
	imlist = []
	for ims in ('{}/*.fits'.format(path_data), '{}/*.fit'.format(path_data), '{}/*.fts'.format(path_data)):
		imlist.extend(sorted(glob.glob(ims)))
	# for inim in glob.glob(path_data+'/*.fit*'):
	for inim in imlist:
		print(inim)
		tool_tbd.puthdr(inim, hdrkey='OBS', hdrval=obs, hdrcomment='observation location')
		data, hdr   = fits.getdata(inim, header=True)
		#	CHECK IMAGETYP HDR
		if hdr['IMAGETYP'] == 'Light': hdr['IMAGETYP'] = 'object'
		if hdr['IMAGETYP'] == 'zero' : hdr['IMAGETYP'] = 'Bias'
		if hdr['IMAGETYP'] == 'Dark' : hdr['IMAGETYP'] = 'dark'
		if hdr['IMAGETYP'] == 'Flat Field' : hdr['IMAGETYP'] = 'FLAT'
		#	CHECK FILTER HDR
		if hdr['IMAGETYP'] in ['object', 'Bias', 'dark',]:
			if ((hdr['FILTER'] == 'U101') | (hdr['FILTER'] == 1)): hdr['FILTER'] = 'U'
			if ((hdr['FILTER'] == 'B102') | (hdr['FILTER'] == 2)): hdr['FILTER'] = 'B'
			if ((hdr['FILTER'] == 'V103') | (hdr['FILTER'] == 3)): hdr['FILTER'] = 'V'
			if ((hdr['FILTER'] == 'R104') | (hdr['FILTER'] == 4)): hdr['FILTER'] = 'R'
			if ((hdr['FILTER'] == 'I105') | (hdr['FILTER'] == 5)): hdr['FILTER'] = 'I'
			if (hdr['FILTER'] == '0') | (hdr['FILTER'] == 0):
				hdr['FILTER'] = 'R'
				print('PLEASE CHECK SOAO FILTER HDR (FILTER:0?)')
		#	CHECK DATE-OBS
		if 'T' not in hdr['DATE-OBS']: hdr['DATE-OBS'] = hdr['DATE-OBS']+'T'+hdr['TIME-OBS']
		t = Time(hdr['DATE-OBS'], format='isot')
		hdr['JD'] = t.jd
		hdr['MJD'] = t.mjd
		#	CHECK OBJECT HDR
		if 'ngc' in hdr['OBJECT']:
			while len(hdr['OBJECT'])<7:
				head = hdr['OBJECT'][0:3]
				tail = hdr['OBJECT'][3: ]
				tail = '0'+tail
				hdr['OBJECT'] = head+tail
		hdr['OBJECT'] = hdr['OBJECT'].upper()
		#	CHECK USER SETTING HDR
		for i in range(len(hdrtbl)):
			key = hdrtbl['key'][i]
			val = hdrtbl['val'][i]
			newval = hdrtbl['newval'][i]
			if hdr[key] == val:
				# print(hdr[key], key, newval)
				hdr[key] = newval
		fits.writeto(inim, data, hdr, overwrite=True)
		if hdr['IMAGETYP'] == 'object':
			objfilterlist.append(hdr['FILTER'])
			objexptimelist.append(hdr['EXPTIME'])
		if hdr['IMAGETYP'].upper() == 'FLAT':
			flatfilterlist.append(hdr['FILTER'])
		if hdr['IMAGETYP'].upper() == 'DARK':
			darkexptimelist.append(hdr['EXPTIME'])
	#	OBJECT FILTER, FLAT FILTER
	objfilterlist = list(set(objfilterlist))#;		objfilterlist.sort()
	objexptimelist = list(set(objexptimelist))#;		objexptimelist.sort()
	flatfilterlist = list(set(flatfilterlist))#;		flatfilterlist.sort()
	darkexptimelist = list(set(darkexptimelist))#;	darkexptimelist.sort()
	return objfilterlist, objexptimelist, flatfilterlist, darkexptimelist, t