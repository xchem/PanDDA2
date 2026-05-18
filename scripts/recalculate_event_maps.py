import argparse
import pathlib
import re

import numpy as np
import gemmi


def recalculate_event_map(event_map_path,):
    name = event_map_path.name
    bdc = float(re.search(
        r'1-BDC_([^_]+)',
        name,
    )[1])

    dtag_dir = event_map_path.parent
    dtag = dtag_dir.name

    mean_map_path = dtag_dir / f'{dtag}-ground-state-average-map.native.ccp4'
    xmap_path = (dtag_dir / 'xmap.ccp4').resolve()

    mean_map_ccp4 = gemmi.read_ccp4_map(str(mean_map_path), )
    mean_map_ccp4.setup(0.0)
    xmap_ccp4 = gemmi.read_ccp4_map(str(xmap_path), )
    xmap_ccp4.setup(0.0)

    mean_map_grid = mean_map_ccp4.grid
    xmap_grid = xmap_ccp4.grid

    mean_map_array = np.array(mean_map_grid, copy=False)
    xmap_array = np.array(xmap_grid, copy=False)

    event_map_array = (xmap_array - ((1-bdc) * mean_map_array)) / (bdc)

    event_map_grid = gemmi.FloatGrid(*[xmap_grid.nu, xmap_grid.nv, xmap_grid.nw])
    event_map_grid.spacegroup = gemmi.find_spacegroup_by_name("P 1")
    event_map_grid.set_unit_cell(xmap_grid.unit_cell)
    event_map_grid_array = np.array(event_map_grid, copy=False)
    event_map_grid_array[:, :, :] = event_map_array[:, :, :]

    return event_map_grid


def write_map(grid, path):

    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = grid
    ccp4.update_ccp4_header()
    ccp4.write_ccp4_map(str(path))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog='RecalculateEventMap',
                    description='Recaclulate event maps to cover the unit cell',
                    )
    
    parser.add_argument(
        'event_map_path', 
        # help='The event map to recalculate',
    )

    parser.add_argument(
        'output_map_path',
        # help='Path for output',
    )
    

    args = parser.parse_args()

    path = pathlib.Path(args.event_map_path)

    grid = recalculate_event_map(
        path,
    )

    write_map(grid, args.output_map_path)