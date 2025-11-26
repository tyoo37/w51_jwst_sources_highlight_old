import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
from astropy.nddata import Cutout2D
from astropy import units as u
import os

import sys
sys.path.append('/home/t.yoo/Paths')
import Paths.Paths as paths

Path = paths.filepaths()

image_filenames ={
    "f140m": "/orange/adamginsburg/jwst/w51/F140M/pipeline/jw06151-o001_t001_nircam_clear-f140m-merged_i2d.fits",
    "f150w": "/orange/adamginsburg/jwst/w51/F150W/pipeline/jw06151-o001_t001_nircam_clear-f150w-merged_i2d.fits",
    "f162m": "/orange/adamginsburg/jwst/w51/F162M/pipeline/jw06151-o001_t001_nircam_clear-f162m-merged_i2d.fits",
    "f182m": "/orange/adamginsburg/jwst/w51/F182M/pipeline/jw06151-o001_t001_nircam_clear-f182m-merged_i2d.fits",
    "f187n": "/orange/adamginsburg/jwst/w51/F187N/pipeline/jw06151-o001_t001_nircam_clear-f187n-merged_i2d.fits",
    "f210m": "/orange/adamginsburg/jwst/w51/F210M/pipeline/jw06151-o001_t001_nircam_clear-f210m-merged_i2d.fits",
    "f335m": "/orange/adamginsburg/jwst/w51/F335M/pipeline/jw06151-o001_t001_nircam_clear-f335m-merged_i2d.fits",
    "f360m": "/orange/adamginsburg/jwst/w51/F360M/pipeline/jw06151-o001_t001_nircam_clear-f360m-merged_i2d.fits",
    "f405n": "/orange/adamginsburg/jwst/w51/F405N/pipeline/jw06151-o001_t001_nircam_clear-f405n-merged_i2d.fits",
    "f410m": "/orange/adamginsburg/jwst/w51/F410M/pipeline/jw06151-o001_t001_nircam_clear-f410m-merged_i2d.fits", # weird, the filename is different from what is downloaded with the STScI pipeline...
    "f480m": "/orange/adamginsburg/jwst/w51/F480M/pipeline/jw06151-o001_t001_nircam_clear-f480m-merged_i2d.fits",
    "f560w": "/orange/adamginsburg/jwst/w51/F560W/pipeline/jw06151-o002_t001_miri_f560w_i2d.fits",
    "f770w": "/orange/adamginsburg/jwst/w51/F770W/pipeline/jw06151-o002_t001_miri_f770w_i2d.fits",
    "f1000w": "/orange/adamginsburg/jwst/w51/F1000W/pipeline/jw06151-o002_t001_miri_f1000w_i2d.fits",
    "f1280w": "/orange/adamginsburg/jwst/w51/F1280W/pipeline/jw06151-o002_t001_miri_f1280w_i2d.fits",
    "f1500w": "/orange/adamginsburg/jwst/w51/F1500W/pipeline/jw06151-o002_t001_miri_f1500w_i2d.fits",
    "f2100w": "/orange/adamginsburg/jwst/w51/F2100W/pipeline/jw06151-o002_t001_miri_f2100w_i2d.fits",
    
}


reproj_dir = '/orange/adamginsburg/jwst/w51/data_reprojected/'
repr140_image_filenames = {x: y.replace("i2d", "i2d_reprj_f140") for x,y in image_filenames.items()}
repr140_image_filenames = {x: (reproj_dir+os.path.basename(y)) for x,y in repr140_image_filenames.items()}

img1 = fits.getdata(repr140_image_filenames['f187n'])
img2 = fits.getdata(repr140_image_filenames['f360m'])
img3 = fits.getdata(repr140_image_filenames['f405n'])
ref_fits = fits.open(image_filenames['f140m'])
tgt_header = fits.getheader(image_filenames['f140m'], ext=('SCI', 1))

wcs_header = WCS(tgt_header)

CO_moment0_fits = fits.open('/orange/adamginsburg/w51/TaehwaYoo/fitsfile_jwst/moment0.fits')
CO_moment0_wcs = WCS(CO_moment0_fits[0].header, naxis=2)
CO_moment0_image = CO_moment0_fits[0].data
CO_moment0_reproj_image, _ = reproject_interp((CO_moment0_image, CO_moment0_wcs), wcs_header, )

reproj_dir2 = '/orange/adamginsburg/jwst/w51/data_reprojected_for_outflowplot/'
hdu = fits.PrimaryHDU(data=CO_moment0_reproj_image, header=tgt_header)
hdu.writeto(reproj_dir2 + 'COmom0_reprj_f140.fits', overwrite=True)

vla_kuband_fits = fits.open('/orange/adamginsburg/w51/TaehwaYoo/vla/2016paper/W51Ku_C_Aarray_continuum_2048_high_uniform.clean.image.fits')
vla_reproj, _ = reproject_interp((vla_kuband_fits[0].data, WCS(vla_kuband_fits[0].header, naxis=2)), wcs_header)
hdu = fits.PrimaryHDU(data=vla_reproj, header=tgt_header)
hdu.writeto(reproj_dir2 + 'vla_reproj_f140.fits', overwrite=True)

h2k_gtc_image = fits.getdata('/orange/adamginsburg/w51/TaehwaYoo/gtc/adendawson/real_reduction/reduced_images/H2_minus_K.fits')
h2k_gtc_wcs = WCS(fits.getheader('/orange/adamginsburg/w51/TaehwaYoo/gtc/adendawson/real_reduction/reduced_images/H2_minus_K.fits'))
h2k_gtc_reproj_image, _ = reproject_interp((h2k_gtc_image, h2k_gtc_wcs), wcs_header, )
hdu = fits.PrimaryHDU(data=h2k_gtc_reproj_image, header=tgt_header)
hdu.writeto(reproj_dir2 + 'h2k_reproj_f140.fits', overwrite=True)


alma_w51n_b6 = fits.open(Path.w51n_b6_tt0)
alma_w51n_b6_image = alma_w51n_b6[0].data
if len(alma_w51n_b6_image.shape)!=2:
    alma_w51n_b6_image = alma_w51n_b6_image[0, 0, :, :]  # Adjust if the data is 3D (e.g., spectral cube)
alma_w51n_b6_wcs = WCS(alma_w51n_b6[0].header, naxis=2)
if CO_moment0_image.ndim > 2:
    CO_moment0_image = np.squeeze(CO_moment0_image)
CO_moment0_reproj_image_to_alma, _ = reproject_interp(
    (CO_moment0_image, CO_moment0_wcs),
    alma_w51n_b6_wcs,
    shape_out=alma_w51n_b6_image.shape
)

hdu = fits.PrimaryHDU(data=CO_moment0_reproj_image_to_alma, header=alma_w51n_b6[0].header)
hdu.writeto(reproj_dir2 + 'COmom0_reprj_alma_b6.fits', overwrite=True)
