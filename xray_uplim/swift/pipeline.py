"""
xray_uplim.swift.pipeline
--------------------------
Full upper-limit pipeline for Swift XRT.

Readout mode (PC / WT) is auto-detected from the event files present on
disk.  Results are reported for the single XRT instrument.

Public API
----------
run_uplim(**kwargs)          — entry point; builds SwiftConfig and runs
process_observation(cfg)     — full pipeline for one Swift observation
"""

import csv
import os
import warnings
import numpy as np

from .config   import SwiftConfig
from .io       import locate_files, load_events, load_expmap
from .aperture import extract_src_bkg_counts, extract_exposure
from .eef      import compute_swift_eef
from ..coords  import parse_coord, sky_to_evt_pixel, sky_to_img_pixel
from ..statistics import net_count_rate, kraft_upper_limit, gehrels_upper_limit


# =============================================================================
# RESULTS TABLE
# =============================================================================

def _print_results_table(N_src, B_scaled, t_eff, N_bkg_raw, area_ratio,
                          confidence_levels, eef=None):
    """Compute and print upper limits at every confidence level."""
    CR_net, CR_sigma = net_count_rate(N_src, B_scaled, t_eff,
                                       N_bkg_raw, area_ratio)

    print(f"\n  Point estimate  (N_src - B) / t_eff  [NOT an upper limit]")
    print(f"    = ({N_src} - {B_scaled:.1f}) / {t_eff:.1f} s")
    print(f"    = {CR_net:+.4e} cts/s  ±  {CR_sigma:.4e}  (1-sigma Poisson)")
    if CR_net < 0:
        print(f"    (Negative — source aperture below expected background: "
              f"clean non-detection)")

    print(f"\n  Upper limits:")
    if eef is not None:
        header = (
            f"  {'CL':>8}  {'Net CR':>13}  "
            f"{'Kraft S_ul':>10}  {'Kraft CR_ap':>13}  {'Kraft CR_tot':>13}  "
            f"{'Geh S_ul':>10}  {'Geh CR_ap':>13}  {'Geh CR_tot':>13}"
        )
    else:
        header = (
            f"  {'CL':>8}  {'Net CR':>13}  "
            f"{'Kraft S_ul':>10}  {'Kraft CR_ap':>13}  "
            f"{'Geh S_ul':>10}  {'Geh CR_ap':>13}"
        )
    divider = "  " + "-" * (len(header) - 2)
    print(header)
    print(divider)

    results = []
    for cl in confidence_levels:
        S_k  = kraft_upper_limit(N_src, B_scaled, cl)
        S_g  = gehrels_upper_limit(N_src, B_scaled, cl)
        CR_k_ap = S_k / t_eff
        CR_g_ap = S_g / t_eff

        CR_k_tot = S_k / (t_eff * eef) if (eef is not None and eef > 0) else None
        CR_g_tot = S_g / (t_eff * eef) if (eef is not None and eef > 0) else None

        results.append({
            'cl':                  cl,
            'CR_net':              CR_net,
            'CR_sigma':            CR_sigma,
            'S_kraft':             S_k,
            'CR_kraft_aperture':   CR_k_ap,
            'CR_kraft_total':      CR_k_tot,
            'S_gehrels':           S_g,
            'CR_gehrels_aperture': CR_g_ap,
            'CR_gehrels_total':    CR_g_tot,
        })

        if eef is not None:
            print(
                f"  {cl:8.4f}  {CR_net:+13.4e}  "
                f"{S_k:10.3f}  {CR_k_ap:13.4e}  {CR_k_tot:13.4e}  "
                f"{S_g:10.3f}  {CR_g_ap:13.4e}  {CR_g_tot:13.4e}"
            )
        else:
            print(
                f"  {cl:8.4f}  {CR_net:+13.4e}  "
                f"{S_k:10.3f}  {CR_k_ap:13.4e}  "
                f"{S_g:10.3f}  {CR_g_ap:13.4e}"
            )

    print(divider)
    if eef is not None:
        print(f"  CR_ap  = aperture count-rate upper limit = S_ul / t_eff.")
        print(f"  CR_tot = EEF-corrected total source rate = S_ul / (t_eff × EEF).")
        print(f"  EEF used: {eef:.4f}")
    else:
        print(f"  CR_ap is the aperture count-rate upper limit.")
        print(f"  EEF correction skipped (psfconst_xrt.fits not found).")

    return results


# =============================================================================
# CSV OUTPUT
# =============================================================================

def _build_csv_rows(mode, e_lo, e_hi, N_src, N_bkg_raw, B_scaled,
                    area_ratio, t_eff, ul_results, eef_info, obsid):
    """Build a list of CSV row dicts (one per confidence level)."""
    rows = []
    for r in ul_results:
        row = {
            'obsid':               obsid,
            'mode':                mode,
            'energy_lo_kev':       e_lo,
            'energy_hi_kev':       e_hi,
            'N_src':               N_src,
            'N_bkg_raw':           N_bkg_raw,
            'B_scaled':            f"{B_scaled:.4f}",
            'area_ratio':          f"{area_ratio:.6f}",
            't_eff_s':             f"{t_eff:.2f}",
            'confidence_level':    r['cl'],
            'CR_net':              f"{r['CR_net']:.6e}",
            'CR_sigma':            f"{r['CR_sigma']:.6e}",
            'S_kraft':             f"{r['S_kraft']:.4f}",
            'CR_kraft_aperture':   f"{r['CR_kraft_aperture']:.6e}",
            'S_gehrels':           f"{r['S_gehrels']:.4f}",
            'CR_gehrels_aperture': f"{r['CR_gehrels_aperture']:.6e}",
            # EEF fields (empty when EEF skipped)
            'theta_arcmin':        '',
            'eef':                 '',
            'energy_kev':          '',
            'psf_file':            '',
            'eef_extrapolated':    '',
            'eef_capped':          '',
            'CR_kraft_total':      '',
            'CR_gehrels_total':    '',
        }
        if eef_info is not None:
            row['theta_arcmin']     = f"{eef_info['theta_arcmin']:.4f}"
            row['eef']              = f"{eef_info['eef']:.6f}"
            row['energy_kev']       = f"{eef_info['energy_kev']:.3f}"
            row['psf_file']         = os.path.basename(eef_info['psf_file'])
            row['eef_extrapolated'] = str(eef_info['extrapolated'])
            row['eef_capped']       = (f"{eef_info['eef_capped']:.6f}"
                                       if eef_info['eef_capped'] is not None
                                       else '')
            if r['CR_kraft_total'] is not None:
                row['CR_kraft_total']   = f"{r['CR_kraft_total']:.6e}"
                row['CR_gehrels_total'] = f"{r['CR_gehrels_total']:.6e}"
        rows.append(row)
    return rows


def write_results_csv(rows, out_dir, obsid):
    """Write upper-limit results to CSV."""
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, f"swift_uplim_{obsid}.csv")

    fieldnames = [
        'obsid', 'mode', 'energy_lo_kev', 'energy_hi_kev',
        'N_src', 'N_bkg_raw', 'B_scaled', 'area_ratio',
        't_eff_s',
        'theta_arcmin', 'eef', 'energy_kev', 'psf_file',
        'eef_extrapolated', 'eef_capped',
        'confidence_level',
        'CR_net', 'CR_sigma',
        'S_kraft',   'CR_kraft_aperture',   'CR_kraft_total',
        'S_gehrels', 'CR_gehrels_aperture', 'CR_gehrels_total',
    ]

    with open(csv_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames,
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n  Results written to: {csv_path}")
    return csv_path


# =============================================================================
# INTERNAL HELPERS
# =============================================================================

def _evt_pixel_to_sky(cx, cy, evt_hdr):
    """Invert sky_to_evt_pixel: event-file pixel (cx, cy) → (ra_deg, dec_deg)."""
    x_col = y_col = None
    for i in range(1, 300):
        if f'TTYPE{i}' not in evt_hdr:
            break
        name = evt_hdr[f'TTYPE{i}'].strip().upper()
        if name == 'X':
            x_col = i
        elif name == 'Y':
            y_col = i

    crpx_x = float(evt_hdr[f'TCRPX{x_col}'])
    crvl_x = float(evt_hdr[f'TCRVL{x_col}'])
    cdlt_x = float(evt_hdr[f'TCDLT{x_col}'])
    crpx_y = float(evt_hdr[f'TCRPX{y_col}'])
    crvl_y = float(evt_hdr[f'TCRVL{y_col}'])
    cdlt_y = float(evt_hdr[f'TCDLT{y_col}'])

    cos_dec = np.cos(np.radians(crvl_y))
    ra_deg  = crvl_x + (cx - crpx_x) * cdlt_x / cos_dec
    dec_deg = crvl_y + (cy - crpx_y) * cdlt_y
    return ra_deg, dec_deg


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def process_observation(cfg: SwiftConfig):
    """
    Full extraction and upper-limit calculation for one Swift XRT observation.

    Steps
    -----
    1.  Locate event file and exposure map (auto-detect PC / WT mode).
    2.  Load and filter events (PI, grade).
    3.  Load exposure map.
    4.  Convert source RA/Dec to event-file and exposure-map pixel coords.
    5.  (Optional) Open interactive region selector GUI.
    6.  Extract source and background counts.
    7.  Compute effective exposure from exposure map.
    8.  Compute EEF from psfconst_xrt.fits (King+Gaussian PSF model).
    9.  Print results table.
    10. Save diagnostic plots.
    11. Write CSV.

    Parameters
    ----------
    cfg : SwiftConfig (validated before calling)

    Returns
    -------
    dict with keys:
        mode, N_src, N_bkg_raw, B_scaled, area_ratio,
        net_counts, t_eff_s, exp_stats, ul, energy, eef_info, csv_rows
    """
    e_lo, e_hi = cfg.resolve_energy_band()
    out_dir    = os.path.join(cfg.data_dir, "ul_products")
    os.makedirs(out_dir, exist_ok=True)

    print(f"\n{'='*70}")
    print(f"  Swift XRT")
    print(f"{'='*70}")

    # -- Step 1: locate files -------------------------------------------------
    evt_file, exp_file, mode = locate_files(cfg)

    # -- Step 2: load events --------------------------------------------------
    events, evt_hdr, pi_lo, pi_hi = load_events(cfg, evt_file, mode)

    # -- Step 3: load exposure map --------------------------------------------
    if exp_file is not None:
        exp_data, exp_hdr = load_expmap(exp_file)
    else:
        exp_data = exp_hdr = None
        warnings.warn(
            "No exposure map found — effective exposure will be estimated "
            "from the event file ONTIME header keyword.",
            UserWarning, stacklevel=2)

    # -- Step 4: source pixel position ----------------------------------------
    src_coord = parse_coord(cfg.ra, cfg.dec)

    cx_evt, cy_evt, pscale_evt = sky_to_evt_pixel(
        src_coord.ra.deg, src_coord.dec.deg, evt_hdr)

    evt_x = np.asarray(events['X'], dtype=float)
    evt_y = np.asarray(events['Y'], dtype=float)

    print(f"\n  Event X range       : [{evt_x.min():.0f}, {evt_x.max():.0f}]")
    print(f"  Event Y range       : [{evt_y.min():.0f}, {evt_y.max():.0f}]")
    print(f"  Source pixel (evt)  : ({cx_evt:.1f}, {cy_evt:.1f})")
    print(f"  Pixel scale (evt)   : {pscale_evt:.4f} \"/pix")

    x_ok = evt_x.min() <= cx_evt <= evt_x.max()
    y_ok = evt_y.min() <= cy_evt <= evt_y.max()
    if not (x_ok and y_ok):
        print(f"  !! WARNING: source pixel is OUTSIDE the event X/Y range — "
              f"check your coordinates!")
    else:
        print(f"  Source position is inside the event image. Good.")

    # -- Step 5: interactive region selector (optional) -----------------------
    label      = f"XRT-{mode}"
    bkg_cx_evt = cx_evt
    bkg_cy_evt = cy_evt

    if cfg.use_gui:
        from ..region_selector import select_regions_interactive
        print(f"\n  Opening interactive region selector for {label}...")
        sel = select_regions_interactive(
            evt_x, evt_y, cx_evt, cy_evt, pscale_evt, cfg, label)

        cx_evt     = sel['cx']
        cy_evt     = sel['cy']
        bkg_cx_evt = sel['bkg_cx']
        bkg_cy_evt = sel['bkg_cy']

        cfg.src_radius_arcsec = sel['src_radius_arcsec']
        cfg.bkg_radius_arcsec = sel['bkg_radius_arcsec']
        cfg.bkg_inner_factor  = sel['bkg_inner_factor']

        bkg_moved = (abs(bkg_cx_evt - cx_evt) > 1.0 or
                     abs(bkg_cy_evt - cy_evt) > 1.0)
        if bkg_moved:
            try:
                bkg_ra, bkg_dec = _evt_pixel_to_sky(bkg_cx_evt, bkg_cy_evt,
                                                     evt_hdr)
                cfg.bkg_mode = 'manual'
                cfg.bkg_ra   = str(float(bkg_ra))
                cfg.bkg_dec  = str(float(bkg_dec))
                print(f"  [GUI] Background → manual mode: "
                      f"RA={bkg_ra:.5f}  Dec={bkg_dec:.5f}")
            except Exception as exc:
                warnings.warn(
                    f"Could not convert background pixel to RA/Dec ({exc}). "
                    "Falling back to annulus mode.",
                    RuntimeWarning, stacklevel=2)
                cfg.bkg_mode = 'annulus'
                bkg_cx_evt   = cx_evt
                bkg_cy_evt   = cy_evt

    print(f"\n  Src aperture : {cfg.src_radius_arcsec:.1f}\"")
    if cfg.bkg_mode == 'annulus':
        r_in = cfg.src_radius_arcsec * cfg.bkg_inner_factor
        print(f"  Bkg annulus  : {r_in:.1f}\" — {cfg.bkg_radius_arcsec:.1f}\"")
    else:
        print(f"  Bkg circle   : r={cfg.bkg_radius_arcsec:.1f}\"  (manual centre)")

    # -- Step 6: source and background counts ---------------------------------
    print()
    N_src, N_bkg_raw, area_ratio, cx_evt, cy_evt, pscale_evt = \
        extract_src_bkg_counts(events, evt_hdr, cfg, mode,
                               bkg_cx_evt=bkg_cx_evt, bkg_cy_evt=bkg_cy_evt)
    B_scaled = N_bkg_raw * area_ratio

    print(f"  Area ratio   (src / bkg) : {area_ratio:.5f}")
    print(f"  Scaled bkg   B           : {B_scaled:.3f} cts")
    print(f"  Net counts   (N_src - B) : {N_src - B_scaled:.3f} cts")

    # -- Step 7: effective exposure -------------------------------------------
    print()
    if exp_data is not None:
        exp_stats, exp_meta, cx_exp, cy_exp = extract_exposure(
            exp_data, exp_hdr, cfg)

        print(f"\n  -- Exposure statistics ------------------------------------------")
        for key, lbl in [('median',       'Median        [RECOMMENDED]        '),
                         ('mean',         'Mean          [diagnostic]         '),
                         ('psf_weighted', 'PSF-wtd mean  [on-axis diag. only] ')]:
            tag = ' <-- PRIMARY' if key == cfg.exp_stat else ''
            print(f"    {lbl} : {exp_stats[key]/1e3:7.3f} ks{tag}")

        t_eff = exp_stats[cfg.exp_stat]
    else:
        # Fallback: ONTIME from event file header
        ontime = float(evt_hdr.get('ONTIME', evt_hdr.get('EXPOSURE', 0.0)))
        t_eff  = ontime
        exp_stats = {'median': ontime, 'mean': ontime, 'psf_weighted': ontime}
        exp_meta  = None
        cx_exp = cy_exp = None
        warnings.warn(
            f"No exposure map — using ONTIME={ontime:.0f} s from event header. "
            "This ignores bad columns and vignetting. Run xrtexpomap for accuracy.",
            UserWarning, stacklevel=2)

    print(f"\n  Using t_eff = {t_eff/1e3:.3f} ks  ({cfg.exp_stat})")

    # -- Step 8: EEF from psfconst_xrt.fits ----------------------------------
    eef_info = None
    try:
        eef_info = compute_swift_eef(
            cfg, evt_hdr, cfg.src_radius_arcsec, e_lo, e_hi)

        print(f"\n  -- EEF (Encircled Energy Fraction) ----------------------------")
        print(f"    Off-axis angle   : {eef_info['theta_arcmin']:.3f} arcmin")
        print(f"    Pointing         : RA={eef_info['pointing_ra']:.5f}  "
              f"Dec={eef_info['pointing_dec']:.5f}")
        print(f"    PSF file         : {os.path.basename(eef_info['psf_file'])}")
        print(f"    Band-centre E    : {eef_info['energy_kev']:.3f} keV")
        print(f"    EEF at {cfg.src_radius_arcsec:.0f}\"       : {eef_info['eef']:.4f}")
        if eef_info['extrapolated']:
            print(f"    !! Off-axis angle exceeds XRT FOV ({eef_info['theta_arcmin']:.1f}'). "
                  f"EEF capped at 12' = {eef_info['eef_capped']:.4f}")

    except (RuntimeError, FileNotFoundError, KeyError) as exc:
        warnings.warn(
            f"EEF computation skipped: {exc}\n"
            "Place psfconst_xrt.fits in xray_uplim/data/swift/psf/ or set "
            "caldb_dir= to enable EEF-corrected upper limits.",
            UserWarning, stacklevel=2)

    # -- Step 9: results table ------------------------------------------------
    eef_val    = eef_info['eef'] if eef_info is not None else None
    ul_results = _print_results_table(
        N_src, B_scaled, t_eff, N_bkg_raw, area_ratio,
        cfg.confidence_levels, eef=eef_val)

    # -- Step 10: diagnostic plots --------------------------------------------
    if cfg.save_plots:
        _save_plots(evt_x, evt_y, cx_evt, cy_evt, pscale_evt,
                    exp_meta, exp_stats, label, e_lo, e_hi, cfg, out_dir,
                    src_coord, bkg_cx_evt, bkg_cy_evt)

    # -- Step 11: CSV ---------------------------------------------------------
    csv_rows = _build_csv_rows(
        mode, e_lo, e_hi, N_src, N_bkg_raw, B_scaled,
        area_ratio, t_eff, ul_results, eef_info, cfg.obsid)

    return {
        'mode':       mode,
        'N_src':      N_src,
        'N_bkg_raw':  N_bkg_raw,
        'B_scaled':   B_scaled,
        'area_ratio': area_ratio,
        'net_counts': N_src - B_scaled,
        't_eff_s':    t_eff,
        'exp_stats':  exp_stats,
        'ul':         ul_results,
        'energy':     (e_lo, e_hi),
        'eef_info':   eef_info,
        'csv_rows':   csv_rows,
    }


# =============================================================================
# DIAGNOSTIC PLOTS
# =============================================================================

def _save_plots(evt_x, evt_y, cx_evt, cy_evt, pscale_evt,
                exp_meta, exp_stats, label, e_lo, e_hi, cfg, out_dir,
                src_coord, bkg_cx_evt, bkg_cy_evt):
    try:
        from ..plots import radial_profile, exposure_histogram, region_image
    except ImportError:
        warnings.warn("Diagnostic plots skipped (import failed).",
                      RuntimeWarning, stacklevel=2)
        return

    radial_profile(
        evt_x, evt_y, cx_evt, cy_evt, pscale_evt,
        label, e_lo, e_hi, cfg.obsid, cfg, out_dir)

    if exp_meta is not None:
        exposure_histogram(exp_meta, exp_stats, label, cfg, out_dir)

    region_image(
        evt_x, evt_y, cx_evt, cy_evt, pscale_evt,
        label, e_lo, e_hi, cfg.obsid, cfg, out_dir,
        src_ra_deg  = src_coord.ra.deg,
        src_dec_deg = src_coord.dec.deg,
        bkg_cx_evt  = bkg_cx_evt,
        bkg_cy_evt  = bkg_cy_evt)


# =============================================================================
# ENTRY POINT
# =============================================================================

def run_uplim(data_dir, obsid, ra, dec, **kwargs):
    """
    Run the full Swift XRT upper-limit pipeline.

    Parameters
    ----------
    data_dir : str           — observation root directory (contains xrt/)
    obsid    : str           — Swift observation ID (e.g. '03000397004')
    ra       : str or float  — source RA  ("HH:MM:SS" or decimal degrees)
    dec      : str or float  — source Dec ("±DD:MM:SS" or decimal degrees)
    **kwargs : any SwiftConfig field, e.g.
                   energy_band='soft',
                   src_radius_arcsec=20.0,
                   confidence_levels=[0.9973],
                   caldb_dir='/path/to/caldb',
                   save_plots=True

    Returns
    -------
    dict — result from process_observation()
    """
    cfg = SwiftConfig(data_dir=data_dir, obsid=obsid, ra=ra, dec=dec, **kwargs)
    cfg.validate()

    e_lo, e_hi = cfg.resolve_energy_band()
    src_coord  = parse_coord(cfg.ra, cfg.dec)
    out_dir    = os.path.join(cfg.data_dir, "ul_products")

    print("Swift XRT Non-Detection Upper Limit")
    print("=" * 70)
    print(f"Source      :  RA = {src_coord.ra.deg:.6f} deg  "
          f"Dec = {src_coord.dec.deg:.6f} deg")
    if isinstance(cfg.energy_band, tuple):
        band_label = f"{e_lo:.2f}–{e_hi:.2f} keV (custom)"
    else:
        band_label = f"'{cfg.energy_band}'  ({e_lo:.2f}–{e_hi:.2f} keV)"
    print(f"Energy band :  {band_label}")
    print(f"Exp stat    :  {cfg.exp_stat}  (primary)")
    print(f"Bkg mode    :  {cfg.bkg_mode}")
    print(f"Data dir    :  {cfg.data_dir}")
    print(f"Obs ID      :  {cfg.obsid}")
    if cfg.caldb_dir:
        print(f"CALDB       :  {cfg.caldb_dir}")
    else:
        caldb_env = os.environ.get('CALDB', '')
        if caldb_env:
            print(f"CALDB       :  {caldb_env}  ($CALDB)")
        else:
            print(f"CALDB       :  not set — using bundled psfconst_xrt.fits")
    print()

    result = process_observation(cfg)

    # -- Summary --------------------------------------------------------------
    print(f"\n{'='*70}")
    print("  SUMMARY")
    print(f"{'='*70}")
    eef_str = (f"{result['eef_info']['eef']:.3f}"
               if result['eef_info'] is not None else "  N/A")
    ul_row = next((u for u in result['ul'] if u['cl'] >= 0.997),
                  result['ul'][-1])
    print(f"  {'Mode':<8}  {'N_src':>6}  {'B_scaled':>9}  "
          f"{'t_eff (ks)':>11}  {'EEF':>6}  "
          f"{'Kraft CR_ap (3σ)':>18}")
    print("  " + "-" * 68)
    print(f"  {result['mode']:<8}  {result['N_src']:>6}  "
          f"{result['B_scaled']:>9.2f}  "
          f"{result['t_eff_s']/1e3:>11.3f}  "
          f"{eef_str:>6}  "
          f"{ul_row['CR_kraft_aperture']:>18.4e}")
    print()

    # -- Write CSV ------------------------------------------------------------
    write_results_csv(result['csv_rows'], out_dir, cfg.obsid)

    return result
