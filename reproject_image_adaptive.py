from reproject import reproject_adaptive, reproject_interp
from astropy.wcs.utils import proj_plane_pixel_scales

from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
import sys
sys.path.append('/home/t.yoo/Paths')
import Paths.Paths as paths
image_filenames ={
    "f140m": "/orange/adamginsburg/jwst/w51/F140M/pipeline/jw06151-o001_t001_nircam_clear-f140m-merged_i2d.fits",
    "f150w": "/orange/adamginsburg/jwst/w51/F150W/pipeline/jw06151-o001_t001_nircam_clear-f150w-merged_i2d.fits",
    "f162m": "/orange/adamginsburg/jwst/w51/F162M/pipeline/jw06151-o001_t001_nircam_clear-f162m-merged_i2d.fits",
    "f182m": "/orange/adamginsburg/jwst/w51/F182M/pipeline/jw06151-o001_t001_nircam_clear-f182m-merged_i2d.fits",
    "f187n": "/orange/adamginsburg/jwst/w51/F187N/pipeline/jw06151-o001_t001_nircam_clear-f187n-merged_i2d.fits",
    "f210m": "/orange/adamginsburg/jwst/w51/F210M/pipeline/jw06151-o001_t001_nircam_clear-f210m-merged_i2d.fits",
    "f335m": "/orange/adamginsburg/jwst/w51/F300M/pipeline/jw06151-o001_t001_nircam_clear-f335m-merged_i2d.fits",
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

image_sub_filenames = {
    "f405n-f410m": "/orange/adamginsburg/jwst/w51/F405_minus_F410cont_pipeline_v0.1.fits",
    "f410m-f405n": "/orange/adamginsburg/jwst/w51/F410_minus_F405_fractional_bandwidth_pipeline_v0.1.fits",
    "f187n-f182m": "/orange/adamginsburg/jwst/w51/F187_minus_F182cont_pipeline_v0.1.fits",

}
Path = paths.filepaths()
alma_w51n_b3_fits = fits.open(Path.w51n_b3_tt0)
alma_w51n_b6_fits = fits.open(Path.w51n_b6_tt0)
alma_w51e_b3_fits = fits.open(Path.w51e_b3_tt0)
alma_w51e_b6_fits = fits.open(Path.w51e_b6_tt0)

f140m_fits = fits.open(image_filenames["f140m"])['SCI']
wcs_f140m = WCS(f140m_fits.header, naxis=2)

alma_w51n_b6_image = alma_w51n_b6_fits[0].data.squeeze()
alma_w51n_b6_wcs = WCS(alma_w51n_b6_fits[0].header, naxis=2)
alma_w51n_b3_image = alma_w51n_b3_fits[0].data.squeeze()
alma_w51n_b3_wcs = WCS(alma_w51n_b3_fits[0].header, naxis=2)

alma_w51e_b6_image = alma_w51e_b6_fits[0].data.squeeze()
alma_w51e_b6_wcs = WCS(alma_w51e_b6_fits[0].header, naxis=2)
alma_w51e_b3_image = alma_w51e_b3_fits[0].data.squeeze()
alma_w51e_b3_wcs = WCS(alma_w51e_b3_fits[0].header, naxis=2)

from reproject.mosaicking import find_optimal_celestial_wcs
science_images = [
    (alma_w51n_b6_image, alma_w51n_b6_wcs),
    (alma_w51n_b3_image, alma_w51n_b3_wcs),
    (alma_w51e_b6_image, alma_w51e_b6_wcs),
    (alma_w51e_b3_image, alma_w51e_b3_wcs),
]
all_pixels = []
for img, wcs in science_images:
    ny, nx = img.shape
    corners = np.array([[0,0],[nx,0],[0,ny],[nx,ny]])
    world = wcs.all_pix2world(corners, 0)
    pix_in_ref = wcs_f140m.all_world2pix(world, 0)
    all_pixels.append(pix_in_ref)
all_pixels = np.vstack(all_pixels)

# 3. Find bounding box in ref_image pixel coordinates
xmin, ymin = np.floor(all_pixels.min(axis=0)).astype(int)
xmax, ymax = np.ceil(all_pixels.max(axis=0)).astype(int)
new_nx = xmax - xmin
new_ny = ymax - ymin

# 4. Shift CRPIX so that the new image covers only the science images
new_wcs = wcs_f140m.deepcopy()
new_wcs.wcs.crpix -= [xmin, ymin]
shape_out = (new_ny, new_nx)

# 5. Reproject all science images to this new WCS/shape
alma_w51n_b6_reprojected, _ = reproject_interp((alma_w51n_b6_image, alma_w51n_b6_wcs), new_wcs, shape_out=shape_out)
alma_w51n_b3_reprojected, _ = reproject_interp((alma_w51n_b3_image, alma_w51n_b3_wcs), new_wcs, shape_out=shape_out)
alma_w51e_b6_reprojected, _ = reproject_interp((alma_w51e_b6_image, alma_w51e_b6_wcs), new_wcs, shape_out=shape_out)
alma_w51e_b3_reprojected, _ = reproject_interp((alma_w51e_b3_image, alma_w51e_b3_wcs), new_wcs, shape_out=shape_out)

# Save as FITS
fits.PrimaryHDU(data=alma_w51n_b6_reprojected, header=new_wcs.to_header()).writeto('/orange/adamginsburg/jwst/w51/reduced_fits_images/alma_w51n_b6_reprojected.fits', overwrite=True)
fits.PrimaryHDU(data=alma_w51n_b3_reprojected, header=new_wcs.to_header()).writeto('/orange/adamginsburg/jwst/w51/reduced_fits_images/alma_w51n_b3_reprojected.fits', overwrite=True)
fits.PrimaryHDU(data=alma_w51e_b6_reprojected, header=new_wcs.to_header()).writeto('/orange/adamginsburg/jwst/w51/reduced_fits_images/alma_w51e_b6_reprojected.fits', overwrite=True)
fits.PrimaryHDU(data=alma_w51e_b3_reprojected, header=new_wcs.to_header()).writeto('/orange/adamginsburg/jwst/w51/reduced_fits_images/alma_w51e_b3_reprojected.fits', overwrite=True)



