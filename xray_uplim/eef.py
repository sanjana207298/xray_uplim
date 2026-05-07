"""
xray_uplim.eef
--------------
Encircled Energy Fraction (EEF) calculation from the NuSTAR CALDB PSF.

Physics
-------
The EEF at aperture radius r_src for a source at off-axis angle theta is

    EEF(r_src, theta) = sum of PSF pixels within r_src
                        ─────────────────────────────
                        sum of all PSF pixels

where the PSF image has been normalised to unit sum before integration.
The full 2D PSF image is used rather than an azimuthally averaged profile,
preserving the genuine asymmetry of the off-axis PSF.

The NuSTAR CALDB tabulates the PSF at off-axis angles 0'-8.5' in steps of
0.5' (18 angles total; Harrison et al. 2013).  At a non-integer angle theta,
pixel-by-pixel bilinear interpolation is performed between the two bracketing
images:

    PSF_theta = PSF_floor(theta) + (theta - floor(theta))
                * (PSF_ceil(theta) - PSF_floor(theta))          (Eq. 12)

For custom energy bands spanning multiple CALDB files, PSFs are combined
with power-law spectral weights assuming photon index Gamma = 2:

    w_i  proportional to  1/E_i,lo - 1/E_i,hi                  (Eq. 13)

The EEF-corrected total-source count-rate upper limit is then

    CR_ul_total = S_ul / (t_eff * EEF(r_src, theta))            (Eq. 11)

CALDB file layout expected
--------------------------
<caldb_dir>/data/nustar/fpm/bcf/psf/
    nu[A/B]psf<YYYYMMDD>v<NNN>.fits

Each file may contain PSFs for one energy sub-band or for the full NuSTAR
band; the code detects which by reading ENERG_LO/ENERG_HI from the header.

Two internal FITS layouts are supported:

  Layout A — BINTABLE with image column (most common CALDB format):
      Extension 1, column 'THETA' : off-axis angles (arcmin), shape (N,)
      Extension 1, column 'IMAGE' : 2-D PSF images,          shape (N, ny, nx)
      (column may also be named 'PSF' or 'PSFINPUT')

  Layout B — Multiple image extensions (one per off-axis angle):
      Each IMAGE extension stores one 2-D PSF; its off-axis angle is given
      by the header keyword THETA (arcmin).

References
----------
Harrison et al. 2013, ApJ 770, 103
"""

import os
import glob
import warnings
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u


PSF_MAX_OFFAXIS     = 8.5    # arcmin — maximum tabulated off-axis angle in CALDB

# Default pixel scale for NuSTAR PSF CALDB images (arcsec/pixel)
# Used only as a fallback when CDELT1 is absent from the header.
PSF_DEFAULT_PSCALE  = 2.4588                          # arcsec/pix


# =============================================================================
# OFF-AXIS ANGLE
# =============================================================================

def off_axis_angle(src_ra_deg, src_dec_deg, evt_hdr):
    """
    Compute the source off-axis angle from the EVENTS header pointing keywords.

    The telescope pointing direction is read from RA_NOM/DEC_NOM (standard
    NuSTAR pipeline output) or RA_PNT/DEC_PNT (fallback).

    Parameters
    ----------
    src_ra_deg  : float         — source right ascension (degrees)
    src_dec_deg : float         — source declination (degrees)
    evt_hdr     : fits.Header   — EVENTS extension header

    Returns
    -------
    theta_arcmin : float  — off-axis angle in arcminutes
    pointing_ra  : float  — pointing RA used (degrees)
    pointing_dec : float  — pointing Dec used (degrees)

    Raises
    ------
    KeyError if no recognised pointing keyword is present.
    """
    for ra_key, dec_key in [('RA_NOM',  'DEC_NOM'),
                             ('RA_PNT',  'DEC_PNT'),
                             ('RA_OBJ',  'DEC_OBJ')]:
        if ra_key in evt_hdr and dec_key in evt_hdr:
            pointing_ra  = float(evt_hdr[ra_key])
            pointing_dec = float(evt_hdr[dec_key])
            break
    else:
        raise KeyError(
            "Cannot determine telescope pointing direction: none of "
            "RA_NOM/DEC_NOM, RA_PNT/DEC_PNT found in EVENTS header. "
            "Check that the event file is a standard nupipeline product.")

    src   = SkyCoord(src_ra_deg  * u.deg, src_dec_deg  * u.deg, frame='icrs')
    point = SkyCoord(pointing_ra * u.deg, pointing_dec * u.deg, frame='icrs')
    theta = float(src.separation(point).to(u.arcmin).value)

    return theta, pointing_ra, pointing_dec


# =============================================================================
# CALDB FILE LOCATION
# =============================================================================

def find_psf_files(caldb_dir, module, e_lo_kev, e_hi_kev):
    """
    Locate NuSTAR PSF CALDB file(s) for a given FPM and energy range.

    Searches <caldb_dir>/data/nustar/fpm/bcf/psf/ for files matching
    nu[A/B]psf*.fits.  If ENERG_LO/ENERG_HI header keywords are present,
    only files whose energy range overlaps [e_lo_kev, e_hi_kev] are kept.
    If no energy keywords exist, the file is treated as broadband.

    The returned list is sorted by lower energy bound; if there is only
    one broadband file, a single-element list is returned.

    Parameters
    ----------
    caldb_dir : str or None
        CALDB root directory.  If None the $CALDB environment variable
        is used.
    module    : str    — 'A' or 'B'
    e_lo_kev  : float
    e_hi_kev  : float

    Returns
    -------
    psf_files : list of str        — absolute paths, sorted by lower energy
    e_ranges  : list of (lo, hi)   — energy coverage of each file (keV)

    Raises
    ------
    RuntimeError  if caldb_dir is not given and $CALDB is not set.
    FileNotFoundError  if the PSF directory or files cannot be found.
    """
    if caldb_dir is None:
        caldb_dir = os.environ.get('CALDB', '')
    if not caldb_dir:
        raise RuntimeError(
            "CALDB directory not set. Either pass caldb_dir= to run_uplim() "
            "or set the $CALDB environment variable to the CALDB root path.")

    psf_dir = os.path.join(caldb_dir, 'data', 'nustar', 'fpm', 'bcf', 'psf')
    if not os.path.isdir(psf_dir):
        raise FileNotFoundError(
            f"NuSTAR PSF CALDB directory not found:\n  {psf_dir}\n"
            "Verify that caldb_dir points to the CALDB root (the directory "
            "that contains data/nustar/…).")

    M = module.upper()

    # Prefer energy-specific PSF files (nuA2dpsfen[1-6]_*.fits).
    # Fall back to the broadband 2D PSF (nuA2dpsf*.fits, excluding grppsf).
    en_pattern   = os.path.join(psf_dir, f"nu{M}2dpsfen*_*.fits")
    broad_pattern = os.path.join(psf_dir, f"nu{M}2dpsf[0-9]*.fits")

    matches = sorted(glob.glob(en_pattern))
    if not matches:
        matches = sorted(glob.glob(broad_pattern))
    if not matches:
        raise FileNotFoundError(
            f"No NuSTAR PSF CALDB files found in:\n  {psf_dir}\n"
            f"Tried patterns:\n  {en_pattern}\n  {broad_pattern}\n"
            "Check that the NuSTAR CALDB is correctly installed.")

    # Read energy coverage from each candidate file
    candidates = []
    for path in matches:
        try:
            with fits.open(path) as hdul:
                hdr = hdul[0].header
                # Energy keywords may be in primary or first extension header
                for h in (hdr, hdul[1].header if len(hdul) > 1 else {}):
                    if 'ENERG_LO' in h and 'ENERG_HI' in h:
                        elo = float(h['ENERG_LO'])
                        ehi = float(h['ENERG_HI'])
                        break
                else:
                    # No energy keywords — treat as full NuSTAR band
                    elo, ehi = 0.0, 1000.0
            candidates.append((elo, ehi, path))
        except Exception:
            continue   # skip unreadable files

    if not candidates:
        raise FileNotFoundError(
            f"Could not read any PSF CALDB files from:\n  {psf_dir}")

    # Keep files whose energy range overlaps the requested band
    relevant = [
        (elo, ehi, path) for elo, ehi, path in candidates
        if elo < e_hi_kev and ehi > e_lo_kev
    ]

    if not relevant:
        # Fall back to the latest (alphabetically last) file
        elo, ehi, path = sorted(candidates)[-1]
        relevant = [(elo, ehi, path)]
        warnings.warn(
            f"No PSF CALDB file overlaps energy range "
            f"{e_lo_kev:.1f}–{e_hi_kev:.1f} keV; "
            f"using {os.path.basename(path)} as a broadband fallback.",
            UserWarning, stacklevel=2)

    relevant.sort(key=lambda x: x[0])
    psf_files = [path for _, _, path in relevant]
    e_ranges  = [(elo, ehi) for elo, ehi, _ in relevant]
    return psf_files, e_ranges


# =============================================================================
# PSF IMAGE LOADING
# =============================================================================

def _pixel_scale_from_header(hdr):
    """
    Extract the pixel scale (arcsec/pixel) from a FITS header.

    Tries CDELT1 (degrees/pixel, converted to arcsec), PSCALE (arcsec/pixel),
    and CD1_1 (degrees/pixel).  Falls back to PSF_DEFAULT_PSCALE with a
    warning if none are found.
    """
    if 'CDELT1' in hdr:
        return abs(float(hdr['CDELT1'])) * 3600.0
    if 'PSCALE' in hdr:
        return abs(float(hdr['PSCALE']))
    if 'CD1_1' in hdr:
        return abs(float(hdr['CD1_1'])) * 3600.0
    warnings.warn(
        "PSF pixel scale not found in FITS header (tried CDELT1, PSCALE, "
        f"CD1_1); using default {PSF_DEFAULT_PSCALE}\" /pix.",
        UserWarning, stacklevel=3)
    return PSF_DEFAULT_PSCALE


def load_psf_images(psf_file):
    """
    Load PSF images at all tabulated off-axis angles from a CALDB file.

    Handles two FITS layouts (see module docstring for details).

    Parameters
    ----------
    psf_file : str — path to NuSTAR PSF CALDB FITS file

    Returns
    -------
    psf_cube         : (N_angles, ny, nx) float array
                       Each image is normalised so that its pixel sum = 1.
    thetas           : (N_angles,) float array — off-axis angles (arcmin),
                       sorted ascending.
    pix_scale_arcsec : float — pixel scale (arcsec/pixel)

    Raises
    ------
    ValueError if the PSF images cannot be found or normalise to zero.
    """
    with fits.open(psf_file) as hdul:

        # ---- Layout A: BINTABLE with image column ---------------------------
        ext1 = hdul[1]
        if not ext1.is_image and hasattr(ext1, 'columns'):
            col_names = [c.name.upper() for c in ext1.columns]
            if 'THETA' in col_names:
                tbl    = ext1.data
                thetas = tbl['THETA'].astype(float).ravel()

                # Find the PSF image column by trying common names
                img_col = None
                for candidate in ('IMAGE', 'PSF', 'PSFINPUT', 'PSFIMAGE'):
                    if candidate in col_names:
                        img_col = candidate
                        break
                if img_col is None:
                    raise ValueError(
                        f"BINTABLE in {psf_file} has a THETA column but no "
                        f"recognisable image column.  "
                        f"Columns present: {col_names}")

                images = np.array([row for row in tbl[img_col]], dtype=float)
                pix_scale_arcsec = _pixel_scale_from_header(ext1.header)

                # Done — break out of the with block
                psf_cube, thetas_sorted = _normalise_and_sort(
                    images, thetas, psf_file)
                return psf_cube, thetas_sorted, pix_scale_arcsec

        # ---- Layout B: multiple image extensions ----------------------------
        # NuSTAR CALDB format: each extension is named "2DPSF_Xarcmin" where
        # X is the off-axis angle (e.g. "2DPSF_0.5ARCMIN").  The angle is
        # NOT stored as a header keyword — it must be parsed from the name.
        import re
        thetas_list = []
        images_list = []
        pix_scale_arcsec = None

        for ext in hdul[1:]:
            if ext.data is None or ext.data.ndim < 2:
                continue

            # Try header keywords first, then parse from extension name
            theta_val = ext.header.get('THETA',
                        ext.header.get('OFFAXIS',
                        ext.header.get('THETA_NOM', None)))

            if theta_val is None:
                m = re.match(r'2DPSF_([0-9.]+)ARCMIN', ext.name.upper())
                if m:
                    theta_val = float(m.group(1))

            if theta_val is None:
                continue   # cannot determine off-axis angle for this extension

            thetas_list.append(float(theta_val))
            images_list.append(ext.data.astype(float))
            if pix_scale_arcsec is None:
                pix_scale_arcsec = _pixel_scale_from_header(ext.header)

        if not thetas_list:
            raise ValueError(
                f"No PSF images found in {psf_file}. "
                "Expected extensions named '2DPSF_Xarcmin' or a BINTABLE "
                "with a THETA column. Check the CALDB file format.")

        if pix_scale_arcsec is None:
            pix_scale_arcsec = PSF_DEFAULT_PSCALE

        images = np.array(images_list, dtype=float)
        thetas = np.array(thetas_list, dtype=float)
        psf_cube, thetas_sorted = _normalise_and_sort(images, thetas, psf_file)
        return psf_cube, thetas_sorted, pix_scale_arcsec


def _normalise_and_sort(images, thetas, filename):
    """
    Sort PSF images by off-axis angle and normalise each to unit pixel sum.

    Returns
    -------
    images : sorted, normalised (N, ny, nx) array
    thetas : sorted (N,) array
    """
    order  = np.argsort(thetas)
    thetas = thetas[order]
    images = images[order]

    for i in range(len(images)):
        s = images[i].sum()
        if s <= 0:
            raise ValueError(
                f"PSF image at theta = {thetas[i]:.1f}' in {filename} "
                "sums to zero or negative; the CALDB file may be corrupt.")
        images[i] /= s

    return images, thetas


# =============================================================================
# PSF INTERPOLATION
# =============================================================================

def interpolate_psf(theta_arcmin, psf_cube, thetas):
    """
    Return a 2D PSF image at arbitrary off-axis angle theta via pixel-by-pixel
    bilinear interpolation between the two bracketing tabulated images.

        PSF_theta = PSF_floor + frac * (PSF_ceil - PSF_floor)   (Eq. 12)

    where frac = theta - floor(theta).

    If theta > max(thetas) (beyond the 8.5' CALDB limit), two PSFs are
    returned: one evaluated at the maximum tabulated angle (capped), and
    one extrapolated linearly from the last two tabulated images.
    A UserWarning is emitted in this case.

    Parameters
    ----------
    theta_arcmin : float
    psf_cube     : (N, ny, nx) array — normalised PSF images
    thetas       : (N,) array        — tabulated off-axis angles (arcmin)

    Returns
    -------
    psf_primary  : (ny, nx) float array  — primary PSF (interpolated, or
                                           capped at max(thetas) if beyond limit)
    psf_extrap   : (ny, nx) or None      — linearly extrapolated PSF; only set
                                           when theta > max(thetas), else None
    extrapolated : bool
    """
    theta_max = float(thetas[-1])

    # On-axis or below minimum tabulated angle
    if theta_arcmin <= float(thetas[0]):
        return psf_cube[0].copy(), None, False

    # Within the tabulated range — interpolate
    if theta_arcmin < theta_max:
        idx_hi = int(np.searchsorted(thetas, theta_arcmin, side='right'))
        idx_lo = idx_hi - 1
        frac   = (theta_arcmin - thetas[idx_lo]) / (thetas[idx_hi] - thetas[idx_lo])
        psf    = psf_cube[idx_lo] + frac * (psf_cube[idx_hi] - psf_cube[idx_lo])
        psf    = _renorm(psf)
        return psf, None, False

    # Exactly at the maximum tabulated angle
    if theta_arcmin == theta_max:
        return psf_cube[-1].copy(), None, False

    # Beyond the CALDB limit — cap + extrapolate
    warnings.warn(
        f"Source off-axis angle {theta_arcmin:.2f}' exceeds the NuSTAR "
        f"CalDB PSF limit ({theta_max:.0f}'). "
        f"EEF will be reported at {theta_max:.0f}' (capped) and via linear "
        f"extrapolation from the {thetas[-2]:.0f}'–{theta_max:.0f}' slope. "
        "Both values are returned; the capped value is used as the primary.",
        UserWarning, stacklevel=2)

    psf_capped = psf_cube[-1].copy()

    # Linear extrapolation: slope is (PSF[-1] - PSF[-2]) per arcmin
    d_theta_tab   = float(thetas[-1] - thetas[-2])
    slope         = (psf_cube[-1] - psf_cube[-2]) / d_theta_tab
    delta         = theta_arcmin - theta_max
    psf_extrap    = psf_capped + delta * slope
    psf_extrap    = np.clip(psf_extrap, 0.0, None)
    psf_extrap    = _renorm(psf_extrap)

    return psf_capped, psf_extrap, True


def _renorm(arr):
    """Renormalise a 2-D array to unit sum; clip negatives first."""
    arr = np.clip(arr, 0.0, None)
    s   = arr.sum()
    return arr / s if s > 0 else arr


# =============================================================================
# SPECTRAL WEIGHTING (Eq. 13)
# =============================================================================

def _spectral_weight(e_lo_kev, e_hi_kev, gamma=2.0):
    """
    Power-law spectral weight for combining PSFs across energy sub-bands.

        w  proportional to  integral_{e_lo}^{e_hi} E^{-Gamma} dE

    For Gamma = 2:  w = 1/e_lo - 1/e_hi   (Eq. 13)
    For Gamma = 1:  w = ln(e_hi/e_lo)
    For other Gamma: general formula.

    Parameters
    ----------
    e_lo_kev, e_hi_kev : float
    gamma              : float — photon index (default 2)

    Returns
    -------
    w : float — unnormalised weight (always >= 0)
    """
    if e_lo_kev <= 0 or e_hi_kev <= e_lo_kev:
        return 0.0
    if abs(gamma - 1.0) < 1e-9:
        return float(np.log(e_hi_kev / e_lo_kev))
    return float(
        (e_lo_kev**(1.0 - gamma) - e_hi_kev**(1.0 - gamma)) / (gamma - 1.0))


# =============================================================================
# EEF INTEGRATION
# =============================================================================

def integrate_eef(psf_image, r_src_arcsec, pix_scale_arcsec):
    """
    Numerically integrate the EEF within a circular aperture from a 2D PSF.

    EEF = sum of PSF pixel values within r_src_arcsec of the image centre
          divided by the total pixel sum.  Since psf_image is pre-normalised
          to unit sum, this simplifies to:

        EEF = sum(psf_image[circle])

    Parameters
    ----------
    psf_image        : (ny, nx) float array  — normalised (sums to 1)
    r_src_arcsec     : float — source aperture radius (arcsec)
    pix_scale_arcsec : float — PSF pixel scale (arcsec/pixel)

    Returns
    -------
    eef : float in [0, 1]
    """
    ny, nx    = psf_image.shape
    cy, cx    = (ny - 1) / 2.0, (nx - 1) / 2.0
    r_pix     = r_src_arcsec / pix_scale_arcsec

    iy, ix = np.ogrid[:ny, :nx]
    dist   = np.sqrt((ix - cx)**2 + (iy - cy)**2)
    mask   = dist <= r_pix

    eef = float(psf_image[mask].sum())
    return float(np.clip(eef, 0.0, 1.0))


# =============================================================================
# MULTI-FILE COMBINATION
# =============================================================================

def _combine_psf_files(psf_files, e_ranges, theta_arcmin, e_lo_kev, e_hi_kev,
                       gamma=2.0):
    """
    Build a single effective PSF by combining multiple energy-band CALDB files
    using power-law spectral weights (Gamma = gamma, Eq. 13).

    The weight for file i is proportional to the integral of E^{-gamma} over
    the overlap of that file's energy range with the requested band [e_lo, e_hi].

    Parameters
    ----------
    psf_files    : list of str        — PSF CALDB paths
    e_ranges     : list of (lo, hi)   — energy range of each file (keV)
    theta_arcmin : float
    e_lo_kev     : float
    e_hi_kev     : float
    gamma        : float              — photon index for spectral weighting (default 2)

    Returns
    -------
    psf_combined     : (ny, nx) float array — normalised weighted PSF
    pix_scale_arcsec : float
    thetas           : (N,) array of off-axis angles from the first file
    files_used       : list of file paths that contributed non-zero weight
    """
    weighted_sum     = None
    total_weight     = 0.0
    pix_scale_arcsec = None
    thetas_ref       = None
    files_used       = []

    for path, (elo, ehi) in zip(psf_files, e_ranges):
        # Clamp to the overlap with the requested band
        overlap_lo = max(elo, e_lo_kev)
        overlap_hi = min(ehi, e_hi_kev)
        if overlap_lo >= overlap_hi:
            continue

        w = _spectral_weight(overlap_lo, overlap_hi, gamma=gamma)
        if w <= 0:
            continue

        psf_cube, thetas, pscale = load_psf_images(path)

        if pix_scale_arcsec is None:
            pix_scale_arcsec = pscale
            thetas_ref       = thetas
        elif abs(pscale - pix_scale_arcsec) / pix_scale_arcsec > 0.01:
            warnings.warn(
                f"PSF pixel scales differ across CALDB files "
                f"({pix_scale_arcsec:.4f} vs {pscale:.4f} arcsec/pix). "
                "Using pixel scale from the first file.",
                UserWarning, stacklevel=2)

        psf_at_theta, _, _ = interpolate_psf(theta_arcmin, psf_cube, thetas)

        if weighted_sum is None:
            weighted_sum = w * psf_at_theta
        else:
            if weighted_sum.shape != psf_at_theta.shape:
                raise ValueError(
                    "PSF image shapes differ across CALDB files "
                    f"({weighted_sum.shape} vs {psf_at_theta.shape}). "
                    "Cannot combine PSFs of different sizes.")
            weighted_sum += w * psf_at_theta

        total_weight += w
        files_used.append(path)

    if weighted_sum is None or total_weight == 0:
        raise RuntimeError(
            "No PSF CALDB files contributed to the spectral combination. "
            f"Requested: {e_lo_kev}–{e_hi_kev} keV.  "
            f"Available ranges: {e_ranges}")

    psf_combined = _renorm(weighted_sum / total_weight)
    return psf_combined, pix_scale_arcsec, thetas_ref, files_used


# =============================================================================
# PUBLIC API
# =============================================================================

def compute_eef(src_ra_deg, src_dec_deg, evt_hdr,
                r_src_arcsec, e_lo_kev, e_hi_kev,
                module, caldb_dir=None, gamma=2.0):
    """
    Compute the EEF for a NuSTAR source aperture from the CALDB PSF.

    Steps
    -----
    1. Compute the source off-axis angle theta from EVENTS header pointing.
    2. Locate NuSTAR PSF CALDB file(s) for this FPM and energy range.
    3. Load PSF images at tabulated off-axis angles (0'–8.5').
    4. Interpolate the PSF at theta (bilinear between bracketing images).
       If multiple CALDB files cover the band, combine with Gamma=2 weights.
    5. Integrate EEF within r_src_arcsec of the PSF image centre.

    Parameters
    ----------
    src_ra_deg   : float         — source RA (degrees)
    src_dec_deg  : float         — source Dec (degrees)
    evt_hdr      : fits.Header   — EVENTS extension header
    r_src_arcsec : float         — source aperture radius (arcsec)
    e_lo_kev     : float         — lower energy bound (keV)
    e_hi_kev     : float         — upper energy bound (keV)
    module       : str           — 'A' or 'B'
    caldb_dir    : str or None   — CALDB root directory (default: $CALDB)
    gamma        : float         — photon index for spectral weighting (default 2.0)
                                   Only relevant when multiple CALDB files are
                                   combined (i.e. energy_band = 'full' or a
                                   custom band spanning multiple sub-bands).

    Returns
    -------
    dict with keys:
        eef              float        EEF at r_src_arcsec (primary value)
        theta_arcmin     float        source off-axis angle (arcmin)
        pointing_ra      float        pointing RA used (degrees)
        pointing_dec     float        pointing Dec used (degrees)
        psf_files        list[str]    CALDB file(s) used
        pix_scale_arcsec float        PSF pixel scale (arcsec/pix)
        extrapolated     bool         True if theta > 8.5' (CalDB limit)
        eef_capped       float|None   EEF at theta = 8.5' (only if extrapolated)
        eef_extrap       float|None   EEF from linear extrap (only if extrapolated)

    Notes
    -----
    When extrapolated=True, eef holds the capped value (theta = 8.5') as the
    conservative primary estimate; eef_extrap is the linearly extrapolated
    value.  Both should be reported in the output with an explicit warning.
    """
    # -- Step 1: off-axis angle -----------------------------------------------
    theta, pt_ra, pt_dec = off_axis_angle(
        src_ra_deg, src_dec_deg, evt_hdr)

    # -- Step 2: CALDB file location ------------------------------------------
    psf_files, e_ranges = find_psf_files(caldb_dir, module, e_lo_kev, e_hi_kev)

    # -- Steps 3 & 4: load, interpolate, (optionally combine) -----------------
    extrapolated = False
    psf_extrap   = None

    if len(psf_files) == 1:
        psf_cube, thetas, pscale = load_psf_images(psf_files[0])

        psf_primary, psf_extrap_img, extrapolated = interpolate_psf(
            theta, psf_cube, thetas)
        files_used = psf_files

    else:
        # Multiple files: spectral combination
        psf_primary, pscale, thetas, files_used = _combine_psf_files(
            psf_files, e_ranges, theta, e_lo_kev, e_hi_kev, gamma=gamma)
        psf_extrap_img = None
        # Extrapolation flag: check theta against the standard CalDB limit
        extrapolated = theta > PSF_MAX_OFFAXIS

    # -- Step 5: EEF integration ----------------------------------------------
    eef = integrate_eef(psf_primary, r_src_arcsec, pscale)

    eef_capped = None
    eef_extrap = None
    if extrapolated and psf_extrap_img is not None:
        eef_capped = eef                                  # primary = capped
        eef_extrap = integrate_eef(psf_extrap_img, r_src_arcsec, pscale)

    return {
        'eef':              eef,
        'theta_arcmin':     theta,
        'pointing_ra':      pt_ra,
        'pointing_dec':     pt_dec,
        'psf_files':        files_used,
        'pix_scale_arcsec': pscale,
        'extrapolated':     extrapolated,
        'eef_capped':       eef_capped,
        'eef_extrap':       eef_extrap,
    }
