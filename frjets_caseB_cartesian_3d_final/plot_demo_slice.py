"""Plot one uniform-grid AMRVAC vtuBCCmpi snapshot; no simulation validation."""
from pathlib import Path
import argparse
import re
import struct
import xml.etree.ElementTree as ET

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_snapshot(path):
    data = path.read_bytes()
    marker = re.search(rb'<AppendedData\s+encoding="raw"\s*>\s*_', data)
    if marker is None:
        raise ValueError('Expected raw appended AMRVAC VTU data.')
    root = ET.fromstring(data[:marker.start()] + b'</VTKFile>')
    if (root.get('byte_order') != 'LittleEndian'
            or root.get('header_type', 'UInt32') != 'UInt32'
            or root.get('compressor')):
        raise ValueError('Expected uncompressed LittleEndian VTU with UInt32 headers.')
    start = marker.end()

    def array(element, vtk_type, dtype):
        if element is None or element.get('type') != vtk_type or element.get('format') != 'appended':
            raise ValueError(f'Expected appended {vtk_type} array; use the supplied vtuBCCmpi configuration.')
        pos = start + int(element.get('offset'))
        size = struct.unpack_from('<I', data, pos)[0]
        if size % 4 or pos + 4 + size > len(data):
            raise ValueError('Invalid appended-array length.')
        return np.frombuffer(data, dtype=dtype, count=size // 4, offset=pos + 4)

    coordinates = []
    fields = {name: [] for name in ('rho', 'v2', 'trp1')}
    for piece in root.findall('./UnstructuredGrid/Piece'):
        points = array(piece.find('./Points/DataArray'), 'Float32', '<f4').reshape(-1, 3)
        connectivity = array(piece.find('./Cells/DataArray[@Name="connectivity"]'), 'Int32', '<i4')
        offsets = array(piece.find('./Cells/DataArray[@Name="offsets"]'), 'Int32', '<i4')
        if not np.array_equal(offsets, np.arange(1, len(offsets) + 1) * 8):
            raise ValueError('This plotter requires eight-vertex 3D cells.')
        connectivity = connectivity.reshape(-1, 8)
        coordinates.append(points[connectivity].mean(axis=1))
        for name in fields:
            values = array(piece.find(f'./CellData/DataArray[@Name="{name}"]'), 'Float32', '<f4')
            if len(values) != len(connectivity):
                raise ValueError(f'Cell count mismatch for {name}.')
            fields[name].append(values)
    if not coordinates:
        raise ValueError('No grid pieces found.')
    time = root.find('./UnstructuredGrid/FieldData/DataArray[@Name="TIME"]')
    actual_time = float(time.text.strip()) if time is not None and time.text else None
    return np.concatenate(coordinates), {name: np.concatenate(values) for name, values in fields.items()}, actual_time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', nargs='?', default='output/demo0002.vtu', type=Path)
    parser.add_argument('--output', default='demo_slice.png', type=Path)
    parser.add_argument('--z', default=0.0, type=float, help='Requested z; choose nearest cell-center plane.')
    args = parser.parse_args()
    if not args.input.is_file():
        parser.error(f'Input not found: {args.input}. Run demo.par first or supply a VTU path.')
    coordinates, fields, time = read_snapshot(args.input)
    z_values = np.unique(coordinates[:, 2])
    distance = np.abs(z_values - args.z)
    z = z_values[np.flatnonzero(distance == distance.min())[-1]]
    selected = coordinates[:, 2] == z
    xx = np.unique(coordinates[selected, 0])
    yy = np.unique(coordinates[selected, 1])
    if len(xx) < 2 or len(yy) < 2:
        raise ValueError('The requested plane must contain at least two cells in x and y.')
    ix = np.searchsorted(xx, coordinates[selected, 0])
    iy = np.searchsorted(yy, coordinates[selected, 1])
    cells = iy * len(xx) + ix
    if (len(cells) != len(xx) * len(yy) or len(np.unique(cells)) != len(cells)
            or not np.allclose(np.diff(xx), xx[1] - xx[0])
            or not np.allclose(np.diff(yy), yy[1] - yy[0])):
        raise ValueError('Expected a complete uniform Cartesian plane; use ParaView for AMR output.')

    fraction = np.divide(fields['trp1'], fields['rho'],
                         out=np.full_like(fields['rho'], np.nan), where=fields['rho'] != 0)
    figure, axes = plt.subplots(1, 3, figsize=(12, 6), layout='constrained')
    panels = ((fields['rho'], 'Density'), (fraction, 'Jet fraction f = trp1/rho'),
              (fields['v2'], 'Axial velocity vy'))
    for ax, (values, title) in zip(axes, panels):
        grid = np.full((len(yy), len(xx)), np.nan)
        grid[iy, ix] = values[selected]
        plot = ax.pcolormesh(xx, yy, np.ma.masked_invalid(grid), shading='nearest', cmap='viridis')
        figure.colorbar(plot, ax=ax, shrink=.7)
        ax.set(xlabel='x / rj', ylabel='y / rj', title=title, aspect='equal')
    label = f't = {time:g}' if time is not None else f'snapshot {args.input.name}'
    figure.suptitle(f'3D Cartesian jet: {label}, z = {z:g}\nCell-center slice; early injection configuration')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=160)
    plt.close(figure)
    print(f'Saved {args.output.resolve()} (actual z = {z:g})')


if __name__ == '__main__':
    main()
