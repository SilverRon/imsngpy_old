#	Library
import os
import glob
import numpy as np
import ccdproc
import astroalign as aa
#------------------------------------------------------------
from astropy.time import Time
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.io import fits
from astropy.nddata import CCDData
#------------------------------------------------------------
def scaling_func(arr): return 1/np.ma.median(arr)
#------------------------------------------------------------
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
#------------------------------------------------------------
def align_astroalign(srcim, tgtim, zero=False):
	'''
	http://quatrope.github.io/astroalign/
	'''
	print(f'Align {os.path.basename(srcim)} to {os.path.basename(tgtim)}')
	tdata, thdr = fits.getdata(tgtim, header=True)

	#	Registered image
	rdata, footprint = aa.register(
		fits.getdata(srcim),
		fits.getdata(tgtim),
		fill_value=np.NaN,
		)
	if zero==True:
		tzero = np.median(tdata[~np.isnan(tdata)].flatten())
		rzero = np.median(rdata[~np.isnan(rdata)].flatten())
		zero_offset = tzero-rzero
		print(f'\t--> (Aligned image) = {os.path.basename(srcim)} - ({round(zero_offset, 3)})')
		rdata = rdata+zero_offset
	else:
		print(f'\t--> Skip scaling with zero')
	return CCDData(rdata, unit='adu', header=fits.getheader(srcim))
#------------------------------------------------------------
def imcombine_ccddata(aligned_imlist, fluxscale=False, zpkey=None, nref=None, imlist=None,):
	"""
	"""
	comment = f"""{'-'*60}\n#\tCOMBINE {len(aligned_imlist)} OBJECT FRAMES\n{'-'*60}"""
	print(comment)
	#	Info. for combined image
	if imlist==None:
		imlist = [f"{indata.header['PATHPRC']}/{indata.header['IMNAME']}" for indata in aligned_imlist]
	jd = calc_jd_median(imlist)
	comim = combine_name(imlist)
	#	Print
	for n, inim in enumerate(imlist):
		print(f'[{n}] {os.path.basename(inim)}')
	print()

	#	Flux scaling with ZP
	if (fluxscale==True) & (zpkey!=None) & (nref!=None):
		zp_ref = aligned_imlist[nref].meta[zpkey]
		print("Flux scaling")
		print('-'*60)
		print(f"{nref} {aligned_imlist[nref].meta['IMNAME']} (ZP_ref={zp_ref})")
		print('-'*60)
		for n in range(len(aligned_imlist)):
			if n!=nref:
				zp = aligned_imlist[n].meta[zpkey]
				factor = scale_flux_zp(zp, zp_ref)
				aligned_imlist[n] = CCDData(aligned_imlist[n].data*factor, unit='adu', header=aligned_imlist[n].meta)
				print(f"{n} {aligned_imlist[n].meta['IMNAME']}*{round(factor, 3)} (ZP={zp})")
	else:
		print(f'Skip flux scaling with zeropoint (fluxscale={fluxscale}, zpkey={zpkey})')

	#	Combine
	combiner = ccdproc.Combiner(aligned_imlist, dtype=np.float32)
	comdata = combiner.median_combine(median_func=bn_median)	

	#	Header
	comdata.header = aligned_imlist[0].header
	comdata.header['IMNAME'] = os.path.basename(comim)
	comdata.header['NCOMBINE'] = len(aligned_imlist)
	for n, inim in enumerate(imlist):
		comdata.header[f'COMB{n}'] = os.path.basename(inim)
	comdata.header['DATE-OBS'] = (jd.isot, 'YYYY-MM-DDThh:mm:ss observation start, UT')
	comdata.header['JD'] = (jd.value, 'Julian Date at start of exposure')
	comdata.header['MJD'] = (jd.mjd, 'Modified Julian Date at start of exposure')
	comdata.header['EXPTIME'] = (calc_total_exptime(imlist), 'Total Exposure time [sec]')
	#	Save
	# comdata.write(comim, overwrite=True)
	fits.writeto(comim, comdata.data, header=comdata.meta, overwrite=True)
	print(f'==> {os.path.basename(comim)}')
	return comim
#------------------------------------------------------------
#	Alignment
def gregistering(imlist, refim):
	import alipy
	identifications = alipy.ident.run(
		refim,
		imlist,
		visu=False
		)
	# for inim, id in zip(imlist, identifications):
	for id in identifications:
		if id.ok == True: # i.e., if it worked
			# print("%20s : %20s, flux ratio %.2f" % (id.ukn.name, id.trans, id.medfluxratio))
			print(f"{id.ukn.name} : {id.trans}, flux ratio {round(id.medfluxratio, 2)}" % (id.ukn.name, id.trans, id.medfluxratio))
			params_align = dict(
				filepath = id.ukn.filepath,
				uknstarlist = id.uknmatchstars,
				refstarlist = id.refmatchstars,
				shape = alipy.align.shape(refim),
				outdir = os.path.dirname(refim),
				makepng = False
				)
			alipy.align.irafalign(**params_align)
		else:
			print(f"{id.ukn.name} : No transformation found!")
	outputshape = alipy.align.shape(refim)
#------------------------------------------------------------
def imcombine_iraf(imlist):
	'''
	path_data = '/data3/paek/factory/doao/20201209-1m-IMSNG'
	images_to_align = sorted(glob.glob('/data3/paek/factory/doao/20201209-1m-IMSNG/Calib-DOAO*-R-60.fits'))
	ref_image = '/data3/paek/factory/doao/20201209-1m-IMSNG/Calib-DOAO-NGC6946-20201209-094720-R-60.fits'
	'''
	from pyraf import iraf
	import glob, os
	from astropy.io import fits
	from astropy.time import Time
	import numpy as np

	path_data = os.path.dirname(imlist[0])
	make_list_file(
		imlist,
		outname=f'{path_data}/imcombine.tmp',
		)
	jd = calc_jd_median(imlist)

	comim = combine_name(imlist)

	print(f'#\t{len(imlist)} IMAGE IMCOMBINE --> {os.path.basename(comim)}')
	param_imcomb = dict(
						input=f"@{path_data}/imcombine.tmp",
						output=comim,
						combine="median",
						project="no",
						reject="none",
						scale="none",
						zero="mode",
						)
	iraf.imcombine(**param_imcomb)

	for i, inim in enumerate(imlist):
		fits.setval(comim, keyword=f'IMCOMB{i}', value=os.path.basename(inim), comment=f'Combined image {i}')
	fits.setval(comim, keyword='DATE-OBS', value=jd.isot, comment='YYYY-MM-DDThh:mm:ss observation start, UT')
	fits.setval(comim, keyword='JD', value=jd.value, comment='Julian Date at start of exposure')
	fits.setval(comim, keyword='MJD', value=jd.mjd, comment='Modified Julian Date at start of exposure')
	fits.setval(comim, keyword='EXPTIME', value=exptime, comment='Total Exposure time [sec]')

	return outim
#------------------------------------------------------------
def calc_jd_median(imlist):
	return Time(np.median([fits.getheader(inim)['JD'] for inim in imlist]), format='jd')
#------------------------------------------------------------
def calc_total_exptime(imlist):
	return int(np.sum([fits.getheader(inim)['EXPTIME'] for inim in imlist]))
#------------------------------------------------------------
def make_list_file(imlist, outname,):
	f = open(f'{outname}', 'w')
	for i, inim in enumerate(imlist):
		f.write(f'{inim}\n')
	f.close()
#------------------------------------------------------------
def split_dateobs(dateobs):
	utdate = dateobs.split('T')[0].replace('-', '')
	uttime = dateobs.split('T')[1].replace(':', '')[:6]
	return utdate, uttime
#------------------------------------------------------------
def combine_name(imlist):
	"""
	"""
	#	The first image
	hdr = fits.getheader(imlist[0])
	path_data = os.path.dirname(imlist[0])
	#	Component
	jdmed = calc_jd_median(imlist)
	dateobsmed = jdmed.isot
	utdate, uttime = split_dateobs(dateobsmed)
	exptime = calc_total_exptime(imlist)
	filte = hdr['FILTER']
	obj = hdr['OBJECT']
	obs = hdr['OBSERVAT']

	outim = f'{path_data}/Calib-{obs}-{obj}-{utdate}-{uttime}-{filte}-{exptime}.com.fits'
	return outim
#------------------------------------------------------------
def scale_flux_zp(zp, zp_ref):
	fref_f = 10.**(-0.4*(zp-zp_ref))
	f_fref = 1./fref_f
	return f_fref
#------------------------------------------------------------
def hotpants(inim, refim, outim, outconvim, iu=60000, tu=6000000000, tl=-100000):
	'''
	'''
	com = f"hotpants -c t -n i -iu {iu} -tu {tu} -tl {tl} -v 0 -inim {inim} -tmplim {refim} -outim {outim} -oci {outconvim}"
	return com
#------------------------------------------------------------
def subtraction_routine(tgtim, path_ref):
	'''
	'''
	srchdr = fits.getheader(tgtim)
	obj, filte = srchdr['OBJECT'], srchdr['FILTER']
	rimlist = glob.glob(f"{path_ref}/Ref*{obj}*{filte}*.fits")
	if len(rimlist)>0:
		srcim = rimlist[0]
		#	Alignment
		# srcim = f"{os.path.splitext(refim)[0]}_aligned{os.path.splitext(refim)[1]}"
		srcdata = align_astroalign(
			srcim=srcim,
			tgtim=tgtim,
			zero=False
			)
		refim = f"{os.path.dirname(tgtim)}/aa{os.path.basename(srcim)}"
		fits.writeto(refim, srcdata.data, header=srcdata.meta, overwrite=True)
		#
		outim=f"{os.path.dirname(tgtim)}/hd{os.path.basename(tgtim)}"
		outconvim=f"{os.path.dirname(tgtim)}/hc{os.path.basename(tgtim)}"
		subcom = hotpants(
			inim=tgtim,
			refim=refim,
			outim=outim,
			outconvim=outconvim,
			iu=60000,
			tu=6000000000,
			tl=-100000,
			)
		os.system(subcom)