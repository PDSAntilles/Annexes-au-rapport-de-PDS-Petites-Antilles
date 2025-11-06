# -*- coding: utf8 -*-
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Scan a continuous waveform stream using one or more templates.

:copyright:
    2021-2025 Claudio Satriano <satriano@ipgp.fr>
:license:
    GNU General Public License v3.0 or later
    (https://www.gnu.org/licenses/gpl-3.0-standalone.html)
"""
import logging
import numpy as np
import os
import sys
from obspy import read
import matplotlib.pyplot as plt
from tqdm import tqdm
from obspy.signal.cross_correlation import correlate_template
from scipy.signal import find_peaks
from ..config import config, rq_exit
from ..families import (
    read_families, read_selected_families,
    FamilyNotFoundError
)
from ..waveforms import (
    get_waveform_from_client, cc_waveform_pair, get_arrivals,
    NoWaveformError
)
from ..catalog import RequakeEvent, generate_evid
from .._version import get_versions
logger = logging.getLogger(__name__.rsplit('.', maxsplit=1)[-1])
trace_cache = {}


def _build_event(tr, template, p_arrival_absolute_time):
    """Build metadata for matched event, using metadata from template."""
    try:
        trace_lat = template.stats.sac.stla
        trace_lon = template.stats.sac.stlo
        ev_lat = template.stats.sac.evla
        ev_lon = template.stats.sac.evlo
        ev_depth = template.stats.sac.evdp
        p_arrival, _s_arrival, _distance, _dist_deg = get_arrivals(
            trace_lat, trace_lon, ev_lat, ev_lon, ev_depth)
        orig_time = p_arrival_absolute_time - p_arrival.time
    except Exception:  # pylint: disable=broad-except
        # Here we catch a broad exception because get_arrivals can fail
        # in many ways, and we don't want to stop the scan
        orig_time = p_arrival_absolute_time
        ev_lon = ev_lat = ev_depth = None
    ev = RequakeEvent()
    ev.orig_time = orig_time
    ev.lon = ev_lon
    ev.lat = ev_lat
    ev.depth = ev_depth
    ev.trace_id = tr.id
    ev.evid = generate_evid(orig_time)
    ev.author = f"requake{get_versions()['version']}"
    return ev


def _cc_detection(tr, template, lag_sec):
    """Compute cross-correlation between detected event and template."""
    # shorter trace is zero-padded on both sides
    #   --aaaa--
    #   bbbbbbbb
    d_len = 0.5 * (len(tr) - len(template)) * tr.stats.delta
    lag_sec += d_len
    p_arrival = lag_sec + template.stats.sac.a
    p_arrival_absolute_time = tr.stats.starttime + p_arrival
    t0 = p_arrival_absolute_time - config.cc_pre_P
    t1 = t0 + config.cc_trace_length
    tr2 = tr.copy().trim(t0, t1)
    _, _, cc_max = cc_waveform_pair(tr2, template)
    return cc_max, p_arrival_absolute_time


def _scan_family_template(template, catalog_file, t0, t1):
    trace_id = template.id
    key = f'{t0}_{trace_id}'
    try:
        tr = trace_cache[key]
    except KeyError:
        try:
            tr = get_waveform_from_client(trace_id, t0, t1)
            trace_cache[key] = tr
        except NoWaveformError as err:
            raise NoWaveformError(f'No data for {trace_id} : {t0} - {t1}') from err

    tr = tr.copy()
    tpl = template.copy()

    F_min = max(20.0, 2.2 * float(config.cc_freq_max))

    def can_decimate(fs):
        return (fs / 2.0) >= F_min

    def maybe_decimate(trace):
        fs0 = trace.stats.sampling_rate
        if can_decimate(fs0):
            trace.decimate(factor=2, strict_length=False, no_filter=False)
        return trace

    tr = maybe_decimate(tr)
    tpl = maybe_decimate(tpl)

    def match_fs(tr_hi, tr_lo):
        Fs_hi = tr_hi.stats.sampling_rate
        Fs_lo = tr_lo.stats.sampling_rate
        if abs(Fs_hi - Fs_lo) < 1e-9:
            return
        if abs(Fs_hi / 2.0 - Fs_lo) < 1e-9 and can_decimate(Fs_hi):
            tr_hi.decimate(factor=2, strict_length=False, no_filter=False)
            return
        tr_lo.interpolate(Fs_hi, method="lanczos", a=20)

    if tr.stats.sampling_rate >= tpl.stats.sampling_rate:
        match_fs(tr, tpl)
    else:
        match_fs(tpl, tr)

    fs = tr.stats.sampling_rate  # = tpl.stats.sampling_rate

    def safe_bandpass(trace):
        fs_loc = trace.stats.sampling_rate
        nyq = 0.5 * fs_loc
        fmax_eff = min(float(config.cc_freq_max), 0.98 * nyq)
        fmin_eff = min(float(config.cc_freq_min), 0.95 * fmax_eff)
        trace.filter("bandpass", freqmin=fmin_eff, freqmax=fmax_eff, corners=4, zerophase=True)

    safe_bandpass(tr)
    safe_bandpass(tpl)

    valid = False
    template_arr = tpl.data.astype(np.float64)
    L = len(template_arr)
    x = tr.slice(starttime=t0, endtime=t1).data.astype(np.float64)

    if len(x) > L:
        ncc = correlate_template(x, template_arr, mode='valid', normalize='full', demean=True)
        ncc = np.clip(ncc, -1.0, 1.0)

        nominal_stop = min(t1, t0 + config.time_chunk)
        n_nom = int(np.floor((nominal_stop - t0) * fs))
        n_valid = max(0, n_nom - L + 1)

        n_all = max(0, len(x) - L + 1)
        single_chunk = (float(t1 - t0) <= float(config.time_chunk) + 1e-6)

        if single_chunk:
            pad_time = (L - 1) / fs
            try:
                tr_ext = get_waveform_from_client(trace_id, t0 - pad_time, t1 + pad_time).copy()
                tr_ext = maybe_decimate(tr_ext)

                Fs_ext = tr_ext.stats.sampling_rate
                if abs(Fs_ext - fs) > 1e-9:
                    if abs(Fs_ext / 2.0 - fs) < 1e-9 and can_decimate(Fs_ext):
                        tr_ext.decimate(factor=2, strict_length=False, no_filter=False)
                    else:
                        tr_ext.interpolate(fs, method="lanczos", a=20)

                safe_bandpass(tr_ext)
                fs_ext = tr_ext.stats.sampling_rate
                x_ext = tr_ext.slice(starttime=t0 - pad_time, endtime=t1 + pad_time).data.astype(np.float64)

                if len(x_ext) > L and abs(fs_ext - fs) < 1e-9:
                    ncc_ext = correlate_template(x_ext, template_arr, mode='valid', normalize='full', demean=True)
                    ncc_ext = np.clip(ncc_ext, -1.0, 1.0)

                    i_start = int(np.ceil(pad_time * fs))
                    i_stop = int(np.floor((pad_time + float(t1 - t0) - L / fs) * fs))
                    if i_stop >= i_start:
                        ncc_valid = ncc_ext[i_start:i_stop + 1]
                        n_valid = len(ncc_valid)
                    else:
                        ncc_valid = np.array([], dtype=float)
                else:
                    ncc_valid = ncc[:n_all]
                    n_valid = len(ncc_valid)
            except Exception:
                ncc_valid = ncc[:n_all]
                n_valid = len(ncc_valid)
        else:
            ncc_valid = ncc[:n_valid]

        if n_valid > 0:
            valid = True
            abs_ncc = np.abs(ncc_valid)
            peaks, props = find_peaks(abs_ncc, height=float(config.cc_min), plateau_size=True)

            if 'left_edges' in props and 'right_edges' in props:
                peak_indices = [int((le + re) // 2) for le, re in zip(props['left_edges'], props['right_edges'])]
            else:
                peak_indices = peaks.tolist()

            candidates = []
            for i in peak_indices:
                val = float(ncc_valid[i])
                if val <= config.cc_min:
                    continue
                tn = t0 + i / fs
                candidates.append((val, tn, i))

            two_days_sec = 2 * 24 * 3600
            candidates.sort(key=lambda t: t[0], reverse=True)
            selected = []
            for val, tn, idx in candidates:
                if all(abs(tn - stn) >= two_days_sec for _, stn, _ in selected):
                    selected.append((val, tn, idx))

            for val, tn, idx in sorted(selected, key=lambda t: t[1]):
                print(f"\n\n|-----> NCC = {val:.4f} | Template n° {template.stats.family_number} | Temps absolu : {tn.isoformat()}\n")
                ev = _build_event(tr, tpl, tn)
                catalog_file.write(f'{ev.fdsn_text()}|{val:.2f}\n')
                catalog_file.flush()

    if not valid:
        print("Aucune NCC calculée (signal trop court ou paramètres incohérents).")




def _read_template_from_file():
    """
    Read a template from a file provided by the user.
    """
    families = read_families()
    try:
        family_number = sorted(f.number for f in families)[-1] + 1
    except IndexError:
        family_number = 0
    templates = []
    try:
        tr = read(config.args.template_file)[0]
        tr.stats.family_number = family_number
        templates.append(tr)
    except (FileNotFoundError, TypeError) as msg:
        logger.warning(msg)
    return templates


def _read_templates():
    """
    Read templates from files in the template directory or from a file
    provided by the user.
    """
    if config.args.template_file is not None:
        return _read_template_from_file()
    families = read_selected_families()
    templates = []
    for family in families:
        trace_id = family[0].trace_id
        template_file = f'template{family.number:02d}.{trace_id}.sac'
        template_file = os.path.join(config.template_dir, template_file)
        try:
            tr = read(template_file)[0]
            tr.stats.family_number = family.number
            templates.append(tr)
        except (FileNotFoundError, TypeError) as msg:
            logger.warning(msg)
    return templates


def _template_catalog_files(templates):
    """
    Create a catalog file for each template.

    Note: the returned dictionary contains file pointers, not file names.
    These file pointers must be closed by the caller.

    :param templates: list of templates
    :type templates: list of obspy.Trace objects
    :return: dictionary of file pointers
    :rtype: dict
    """
    catalog_files = {}
    for template in templates:
        template_catalog_dir = os.path.join(
            config.args.outdir, 'template_catalogs'
        )
        if not os.path.exists(template_catalog_dir):
            os.makedirs(template_catalog_dir)
        template_signature =\
            f'{template.stats.family_number:02d}.{template.id}'
        template_catalog_file_name = os.path.join(
            template_catalog_dir, f'catalog{template_signature}.txt'
        )
        catalog_files[template_signature] = open(
            template_catalog_file_name, 'w', encoding='utf-8')
    return catalog_files


def scan_templates():
    """
    Perform cross-correlation on catalog events. 
    Manages progression bar. 
    Manages data format errors.  
    """
    try:
        templates = _read_templates()
    except (FileNotFoundError, FamilyNotFoundError) as msg:
        logger.error(msg)
        rq_exit(1)

    catalog_files = _template_catalog_files(templates)
    start_time = config.template_start_time
    end_time = config.template_end_time
    time_chunk = config.time_chunk
    overlap = config.time_chunk_overlap

    if end_time > start_time:
        total_span = float(end_time - start_time)  
    else:
        total_span = 1.0
    nchunks = max(1, int(np.ceil(total_span / float(time_chunk))))

    pbar = tqdm(total=nchunks, unit='chunks', unit_scale=False, disable=not sys.stderr.isatty())
    pbar.set_description_str("[scan_templates]")
    pbar.set_postfix_str(f"{0:5.1f}%  {start_time.isoformat()} / {end_time.isoformat()}")

    for k in range(nchunks):
        t0 = start_time + k * time_chunk
        chunk_stop = min(t0 + time_chunk, end_time)
        t1 = min(chunk_stop + overlap, end_time)

        for template in templates:
            template_signature = f'{template.stats.family_number:02d}.{template.id}'
            catalog_file = catalog_files[template_signature]
            try:
                _scan_family_template(template, catalog_file, t0, t1)
            except NoWaveformError as msg:
                logger.warning(msg)
                continue
        
        if total_span > 0:
            pct = 100.0 * (float(chunk_stop - start_time) / total_span)  
        else:
            pct = 100.0
        pbar.update(1)
        pbar.set_postfix_str(f"{pct:5.1f}%  {chunk_stop.isoformat()} / {end_time.isoformat()}")
        trace_cache.clear()

    pbar.close()
    for fp in catalog_files.values():
        fp.close()
