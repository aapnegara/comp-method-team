"""Role 4 (Q6/Q7) post-processing: compare limiters and flux solvers.

Reads every runs/<name>/output/<name>*.vtu (AMRVAC vtuBCCmpi, uniform grid),
plus <name>.log, timing.txt and <name>.stdout, and writes to analysis/:

  timeseries.csv   one row per snapshot per run
  summary.csv      one row per run (cost, stability, final diagnostics)
  yhead_<exp>.png  jet-head and bow-shock position versus time
  axis_<exp>.png   on-axis profiles (rho, v_y, f, p) at a common time
  slices_<exp>.png xy slices at z~0 with IDENTICAL colour ranges per column

Diagnostic definitions (identical for every run):
  f        = trp1 / rho   (trp1 in the output holds C = rho f, see mod_usr.t)
  y_head   = largest cell-centre y with f > 0.5        (jet material front)
  y_bow    = largest cell-centre y with |p/p_amb - 1| > 0.05, p_amb = 0.6
             (leading edge of the bow shock / pressure disturbance)
  V_jet    = volume with f > 0.5
  V_mix    = volume with 0.1 < f < 0.9  -> width of the jet/ambient transition.
             At this resolution this is dominated by NUMERICAL diffusion; do not
             interpret it as physical turbulent mixing.
  axis     = mean over the 2x2 cell columns closest to x = z = 0

Usage:  python analyze.py                 (all runs found in runs/)
        python analyze.py --time 2.0      (profiles/slices at t = 2)
        python analyze.py --ymax 4        (zoom slices on the jet)
"""
from pathlib import Path
import argparse
import csv
import re
import struct
import xml.etree.ElementTree as ET

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P_AMB = 0.6
FIELDS = ('rho', 'v1', 'v2', 'v3', 'p', 'trp1')
EXPERIMENTS = {
    'A_limiters': ['baseline96', 'woodward96', 'vanleer96', 'koren96', 'cada3_96'],
    'B_solvers': ['baseline96', 'hll96', 'hllc96', 'roe_tvd96'],
}


# --------------------------------------------------------------------- input
def read_vtu(path):
    """Return time, cell-centre axes (x, y, z) and 3D field arrays [ix, iy, iz]."""
    data = path.read_bytes()
    marker = re.search(rb'<AppendedData\s+encoding="raw"\s*>\s*_', data)
    if marker is None:
        raise ValueError(f'{path}: expected raw appended AMRVAC VTU data')
    root = ET.fromstring(data[:marker.start()] + b'</VTKFile>')
    if root.get('byte_order') != 'LittleEndian' or root.get('compressor'):
        raise ValueError(f'{path}: expected uncompressed little-endian VTU')
    start = marker.end()

    def array(element, dtype):
        pos = start + int(element.get('offset'))
        size = struct.unpack_from('<I', data, pos)[0]
        return np.frombuffer(data, dtype=dtype, count=size // 4, offset=pos + 4)

    centres, fields = [], {k: [] for k in FIELDS}
    for piece in root.findall('./UnstructuredGrid/Piece'):
        pts = array(piece.find('./Points/DataArray'), '<f4').reshape(-1, 3)
        conn = array(piece.find('./Cells/DataArray[@Name="connectivity"]'), '<i4').reshape(-1, 8)
        centres.append(pts[conn].mean(axis=1))
        for k in FIELDS:
            fields[k].append(array(piece.find(f'./CellData/DataArray[@Name="{k}"]'), '<f4'))
    c = np.concatenate(centres).astype(np.float64)
    tnode = root.find('./UnstructuredGrid/FieldData/DataArray[@Name="TIME"]')
    time = float(tnode.text.strip())

    # Map scattered cell centres onto a uniform 3D array (round to kill float32 noise).
    axes, idx = [], []
    for d in range(3):
        u, inv = np.unique(np.round(c[:, d], 5), return_inverse=True)
        axes.append(u)
        idx.append(inv)
    shape = tuple(len(a) for a in axes)
    if np.prod(shape) != len(c):
        raise ValueError(f'{path}: not a complete uniform grid (AMR output?)')
    out = {}
    for k in FIELDS:
        arr = np.empty(shape)
        arr[idx[0], idx[1], idx[2]] = np.concatenate(fields[k])
        out[k] = arr
    out['f'] = np.divide(out['trp1'], out['rho'], out=np.zeros(shape), where=out['rho'] > 0)
    return time, axes, out


def read_log(path):
    """Steps, final time and mean cell-updates/s/core from the AMRVAC .log file."""
    if not path.is_file():
        return {}
    lines = [l for l in path.read_text().splitlines() if l.strip()]
    header = lines[0].replace('|', ' ').split()
    rows = [l.replace('|', ' ').split() for l in lines[1:]]
    rows = [r for r in rows if len(r) >= 3]
    if not rows:
        return {}
    res = {'nsteps': int(float(rows[-1][0])), 't_reached': float(rows[-1][1])}
    # 'Cell_Updates /second/core' is quoted with a space -> locate it from the right.
    ncol_right = 6  # Xload Xmemory CellUpd ActiveBlocks Wct/CodeTime TimeToFinish
    if "'Cell_Updates" in header:
        cu = [float(r[-ncol_right + 2]) for r in rows if len(r) >= ncol_right]
        cu = [v for v in cu if v > 0]
        if cu:
            res['cell_updates_per_s_core'] = float(np.mean(cu))
    return res


def read_kv(path):
    """key=value pairs from timing.txt; wall_seconds is SUMMED over resumed parts."""
    if not path.is_file():
        return {}
    kv, wall = {}, 0
    for line in path.read_text().splitlines():
        for tok in line.split():
            if '=' in tok:
                k, v = tok.split('=', 1)
                kv[k] = v
                if k == 'wall_seconds':
                    wall += int(v)
    if 'wall_seconds' in kv:
        kv['wall_seconds'] = str(wall)
    return kv


def scheme_of(par):
    if not par.is_file():
        return '?', '?', '?'
    txt = par.read_text()
    get = lambda key: (re.search(rf"^\s*{key}\s*=\s*(?:\d+\*)?'([^']+)'", txt, re.M) or [None, '?'])[1]
    return get('time_stepper'), get('flux_scheme'), get('limiter')


def failure_note(stdout):
    if not stdout.is_file():
        return ''
    txt = stdout.read_text(errors='replace')
    for pat in ('small value', 'Error', 'ERROR', 'mpistop', 'NaN', 'Segmentation'):
        m = re.search(rf'.*{pat}.*', txt)
        if m:
            return m.group(0).strip()[:120]
    return ''


# ---------------------------------------------------------------- diagnostics
def diagnostics(axes, w):
    x, y, z = axes
    dV = np.prod([a[1] - a[0] for a in axes])
    ygrid = np.broadcast_to(y[None, :, None], w['rho'].shape)
    jet = w['f'] > 0.5
    dist = np.abs(w['p'] / P_AMB - 1) > 0.05
    return {
        'y_head': float(ygrid[jet].max()) if jet.any() else 0.0,
        'y_bow': float(ygrid[dist].max()) if dist.any() else 0.0,
        'rho_min': float(w['rho'].min()), 'rho_max': float(w['rho'].max()),
        'p_min': float(w['p'].min()), 'p_max': float(w['p'].max()),
        'vy_max': float(w['v2'].max()),
        'V_jet': float(jet.sum() * dV),
        'V_mix': float(((w['f'] > 0.1) & (w['f'] < 0.9)).sum() * dV),
    }


def axis_profile(axes, w, key):
    x, _, z = axes
    ix = np.argsort(np.abs(x))[:2]
    iz = np.argsort(np.abs(z))[:2]
    return w[key][np.ix_(ix, np.arange(w[key].shape[1]), iz)].mean(axis=(0, 2))


# ------------------------------------------------------------------- analysis
def analyse_run(rundir):
    name = rundir.name
    snaps = sorted((rundir / 'output').glob(f'{name}[0-9][0-9][0-9][0-9].vtu'))
    series = []
    for s in snaps:
        t, axes, w = read_vtu(s)
        series.append({'run': name, 'snapshot': s.name, 't': t, **diagnostics(axes, w)})
    stepper, flux, lim = scheme_of(rundir / f'{name}.par')
    timing = read_kv(rundir / 'timing.txt')
    info = {'run': name, 'time_stepper': stepper, 'flux_scheme': flux, 'limiter': lim,
            'exit_status': timing.get('exit_status', '?'),
            'wall_seconds': timing.get('wall_seconds', '?'),
            'n_snapshots': len(snaps),
            'failure': failure_note(rundir / f'{name}.stdout'),
            **read_log(rundir / 'output' / f'{name}.log')}
    if series:
        info.update({k: v for k, v in series[-1].items() if k not in ('run', 'snapshot')})
    return info, series, snaps


def snapshot_at(snaps, t_target):
    """Snapshot whose time is closest to t_target (reads headers only)."""
    best, bestdt = None, np.inf
    for s in snaps:
        head = s.read_bytes()[:600].decode('latin1')
        t = float(re.search(r'Name="TIME"[^>]*>\s*([-0-9.Ee+]+)', head).group(1))
        if abs(t - t_target) < bestdt:
            best, bestdt = s, abs(t - t_target)
    return best


def plot_experiment(exp, names, results, t_common, outdir, zslice, ymax=None):
    names = [n for n in names if n in results and results[n][2]]
    if not names:
        return
    label = lambda n: f"{n}  ({results[n][0]['flux_scheme']}+{results[n][0]['limiter']})"

    # 1) y_head(t), y_bow(t)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4), layout='constrained')
    for n in names:
        ser = results[n][1]
        t = [r['t'] for r in ser]
        ax[0].plot(t, [r['y_head'] for r in ser], 'o-', ms=3, label=label(n))
        ax[1].plot(t, [r['y_bow'] for r in ser], 'o-', ms=3, label=label(n))
    ax[0].set(xlabel='t  [r_j / c_s,j]', ylabel='y_head  [r_j]', title='Jet head (f > 0.5)')
    ax[1].set(xlabel='t  [r_j / c_s,j]', ylabel='y_bow  [r_j]', title='Bow shock (|p/p_amb-1| > 0.05)')
    ax[0].legend(fontsize=7)
    fig.savefig(outdir / f'yhead_{exp}.png', dpi=150)
    plt.close(fig)

    # Load the common-time snapshot of every run once.
    data = {}
    for n in names:
        s = snapshot_at(results[n][2], t_common)
        data[n] = read_vtu(s)

    # 2) on-axis profiles
    keys = [('rho', 'density', True), ('v2', 'v_y', False), ('f', 'tracer f', False), ('p', 'pressure', True)]
    fig, ax = plt.subplots(1, 4, figsize=(16, 4), layout='constrained')
    for n in names:
        t, axes, w = data[n]
        for a, (k, title, logy) in zip(ax, keys):
            a.plot(axes[1], axis_profile(axes, w, k), label=f'{label(n)}  t={t:.2f}')
            a.set(xlabel='y  [r_j]', title=f'on-axis {title}')
            if logy:
                a.set_yscale('log')
    ax[0].legend(fontsize=7)
    fig.savefig(outdir / f'axis_{exp}.png', dpi=150)
    plt.close(fig)

    # 3) xy slices at z ~ zslice, identical colour range per column
    cols = [('rho', r'$\log_{10}\rho$', 'viridis', np.log10),
            ('f', r'tracer $f$', 'magma', None),
            ('v2', r'$v_y$', 'RdBu_r', None),
            ('p', r'$\log_{10}p$', 'inferno', np.log10)]
    slices = {}
    for n in names:
        t, axes, w = data[n]
        iz = int(np.argmin(np.abs(axes[2] - zslice)))
        slices[n] = {k: (fn(w[k][:, :, iz]) if fn else w[k][:, :, iz]) for k, _, _, fn in cols}
    lims = {}
    for k, _, _, _ in cols:
        allv = np.concatenate([slices[n][k].ravel() for n in names])
        lims[k] = (0.0, 1.0) if k == 'f' else (float(np.nanmin(allv)), float(np.nanmax(allv)))
    vabs = max(abs(lims['v2'][0]), abs(lims['v2'][1]))
    lims['v2'] = (-vabs, vabs)

    fig, ax = plt.subplots(len(names), len(cols), figsize=(3.2 * len(cols), 3.3 * len(names)),
                           layout='constrained', squeeze=False)
    for i, n in enumerate(names):
        t, axes, _ = data[n]
        for j, (k, title, cmap, _) in enumerate(cols):
            im = ax[i, j].pcolormesh(axes[0], axes[1], slices[n][k].T, shading='nearest',
                                     cmap=cmap, vmin=lims[k][0], vmax=lims[k][1])
            ax[i, j].set_aspect('equal')
            if ymax is not None:
                ax[i, j].set_ylim(0, ymax)
            if i == 0:
                ax[i, j].set_title(title)
            if j == 0:
                ax[i, j].set_ylabel(f'{label(n)}\nt={t:.2f}   y [r_j]', fontsize=8)
            if i == len(names) - 1:
                ax[i, j].set_xlabel('x [r_j]')
                fig.colorbar(im, ax=ax[:, j], location='bottom', shrink=0.8)
    fig.suptitle(f'{exp}: z = {zslice:g} slices, common colour scale per column')
    fig.savefig(outdir / f'slices_{exp}.png', dpi=130)
    plt.close(fig)


def write_csv(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader()
        wr.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--root', type=Path, default=Path(__file__).resolve().parent / 'runs')
    ap.add_argument('--out', type=Path, default=Path(__file__).resolve().parent / 'analysis')
    ap.add_argument('--time', type=float, default=None,
                    help='time for profiles/slices (default: latest time reached by ALL runs)')
    ap.add_argument('--z', type=float, default=0.0)
    ap.add_argument('--ymax', type=float, default=None, help='zoom slices to 0 <= y <= ymax')
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    results = {}
    for d in sorted(p for p in args.root.iterdir() if p.is_dir()):
        info, series, snaps = analyse_run(d)
        results[d.name] = (info, series, snaps)
        print(f"{d.name:16s} snapshots={len(snaps):3d}  steps={info.get('nsteps', '?')}  "
              f"t={info.get('t_reached', '?')}  wall={info['wall_seconds']}s  {info['failure']}")
    if not results:
        raise SystemExit(f'No run directories in {args.root}')

    write_csv(args.out / 'summary.csv', [r[0] for r in results.values()])
    write_csv(args.out / 'timeseries.csv', [row for r in results.values() for row in r[1]])

    for exp, names in EXPERIMENTS.items():
        present = [n for n in names if n in results and results[n][1]]
        if not present:
            continue
        t_common = args.time if args.time is not None else min(results[n][1][-1]['t'] for n in present)
        plot_experiment(exp, present, results, t_common, args.out, args.z, args.ymax)
    print(f'Wrote results to {args.out}')


if __name__ == '__main__':
    main()
