"""
Verification stamps for Tier-2 injection test.
Two models:
  A: i1j2k4  Sigma=3e21  M_BE=0.080  T_BE=1.75 K  -- identical to the 18-core test (cross-check)
  B: i2j2k4  Sigma=6e21  M_BE=0.080  T_BE=1.75 K  -- same mass, double cloud column

Stamp format: 73x73 pix, 3"/pix, center at pixel (37,37) 1-based.
Profile: 2D Gaussian, peak = IC_Xbs, FWHM = FWHMX_bs from catalog.
Clipping: zero outside the radius where the Gaussian drops below 1e-3 * peak
(matches the convention used in the existing nc.Xum.bs stamps).
"""
import numpy as np
from astropy.io import fits
from fit_mass import load_catalog

cat = load_catalog('bes_model_params_catalog')
PIX = 3.0; N = 73; CX = CY = N // 2

def gauss_stamp(fwhm_as, peak, n=73, pix=3.0, clip_frac=1e-3):
    sig = (fwhm_as / pix) / 2.3548200450309493
    y, x = np.mgrid[0:n, 0:n] - (n // 2)
    g = peak * np.exp(-(x**2 + y**2) / (2 * sig**2))
    g[g < clip_frac * peak] = 0.0
    return g

def write_stamp(arr, path, wave_um, bunit='MJy/sr'):
    hdu = fits.PrimaryHDU(arr.astype(np.float32))
    h = hdu.header
    h.update({'BZERO':0.,'BSCALE':1.,'BUNIT':bunit,
              'CRPIX1':float(CX+1),'CRPIX2':float(CY+1),
              'CRVAL1':0.,'CRVAL2':0.,
              'CDELT1':-PIX/3600.,'CDELT2':PIX/3600.,
              'CD1_1':-PIX/3600.,'CD1_2':0.,'CD2_1':0.,'CD2_2':PIX/3600.,
              'CTYPE1':'RA---TAN','CTYPE2':'DEC--TAN',
              'CROTA1':0.,'CROTA2':0.,'EQUINOX':2000.,
              'WAVE':float(wave_um),'OBJECT':'RT Model','CREATOR':'make_stamps_v2.py'})
    hdu.writeto(path, overwrite=True)

models = [
    (1,2,4,'3e21_i1j2k4','Sigma=3e21, M_BE=0.080, identical to 18-core test'),
    (2,2,4,'6e21_i2j2k4','Sigma=6e21, M_BE=0.080, same mass higher column'),
]
bands = [('160',160),('250',250),('350',350),('500',500)]

print(f'{"Model":18s} {"band":6s} {"IC(MJy/sr)":12s} {"FWHM(\")":8s} '
      f'{"peak_out":10s} {"nonzero_px":10s} {"sum":10s}')
print('-'*74)

for ii,jj,kk,tag,desc in models:
    m = (cat['i']==ii)&(cat['j']==jj)&(cat['k']==kk)
    R_BE = cat['R_BE_as'][m][0]
    print(f'\n{tag}  {desc}')
    print(f'  R_BE={R_BE:.2f}"  M_SED3bs={cat["M_SED3bs"][m][0]:.4f}  '
          f'f=M_BE/M_SED3bs={cat["M_BE"][m][0]/cat["M_SED3bs"][m][0]:.3f}')
    for b,wave in bands:
        ic   = cat[f'IC{b}bs'][m][0]
        fwhm = cat[f'FWHM{b}bs'][m][0]
        s = gauss_stamp(fwhm, ic)
        path = f'inject_{tag}.{b}um.bs.fits'
        write_stamp(s, path, wave, 'MJy/sr')
        print(f'  {b}um: IC={ic:8.4f}  FWHM={fwhm:6.2f}"  '
              f'peak={s.max():8.4f}  nz={np.count_nonzero(s):4d}  sum={s.sum():9.2f}')
    # surfdens
    ic_sd = cat['ICSDbs'][m][0]; fw_sd = cat['FWHMSDbs'][m][0]
    s_sd  = gauss_stamp(fw_sd, ic_sd)
    write_stamp(s_sd, f'inject_{tag}.surfdens.bs.fits', 0, 'H2/cm^2')
    print(f'  surfdens: IC={ic_sd:.4e}  FWHM={fw_sd:.2f}"')

# verification: compare model A with original nc stamp
print('\n--- Cross-check: model A vs original nc.250um.bs stamp ---')
orig = fits.getdata('nc.250um.bs.r18p2x0.rs3p0as.fits').astype(float)
new  = fits.getdata('inject_3e21_i1j2k4.250um.bs.fits').astype(float)
print(f'  original: peak={orig.max():.5f}  nz={np.count_nonzero(orig)}  sum={orig.sum():.2f}')
print(f'  new:      peak={new.max():.5f}  nz={np.count_nonzero(new)}   sum={new.sum():.2f}')
print(f'  peak diff = {abs(orig.max()-new.max()):.2e} MJy/sr  ({100*abs(orig.max()-new.max())/orig.max():.3f}%)')
print(f'  Note: sum differs because the original nc stamp uses the actual')
print(f'  RT model profile, not a pure Gaussian (wings differ slightly).')
print(f'  The peak and FWHM are what matter for the SED fit and mass.')
