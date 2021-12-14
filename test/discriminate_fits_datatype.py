from astropy.io import fits
# path_data = ''
import glob
imlist = sorted(glob.glob(f'{path_data}/Calib*NGC*0.fits'))

for inim in imlist:
	with fits.open(inim) as hdul:
		print(hdul.info())