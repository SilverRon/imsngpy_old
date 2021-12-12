import time
import fitsio
from astropy.io import fits

path_test = '/data3/paek/factory/master_frames/LOAO/zero'

import glob
imlist = sorted(glob.glob(f'{path_test}/*.fits'))[:100]

t_fitsio = []
t_astro = []

for inim in imlist:
#	fitsio
	st = time.time()
	fitsio.read(inim, header=True)
	t_fitsio.append(time.time()-st)
#	astropy.io
	st = time.time()
	fits.getdata(inim, header=True)
	t_astro.append(time.time()-st)

import numpy as np
t_astro = np.array(t_astro)
t_fitsio = np.array(t_fitsio)

y_astro = np.cumsum(t_astro)
y_fitsio = np.cumsum(t_fitsio)
x = np.arange(1, len(imlist)+1, 1)

plt.close('all')
plt.plot(x, y_astro, mfc='none', color='dodgerblue', label=f'astropy.io.fits ({round(y_astro[-1]/100, 3)} sec/image)',)
plt.plot(x, y_fitsio, mfc='none', color='tomato', label=f'fitsio ({round(y_fitsio[-1]/100, 3)} sec/image)',)
plt.plot(x, y_fitsio-y_astro, mfc='none', color='k', alpha=0.5, label='residual')

plt.title(f'{len(imlist)} images', fontsize=14)
plt.legend(fontsize=14, framealpha=0)
plt.xlabel('# of images', fontsize=20)
plt.ylabel('Time [sec]', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
# plt.yscale('log')
plt.grid('both', ls='--', c='silver', alpha=0.5)
plt.tight_layout()