import gemmi
import fire
import pathlib


def event_map_to_mtz(event_map_path, reference_structure_path, output_mtz_path):
    # Get the resolution of the map
    st = gemmi.read_structure(str(reference_structure_path))
    resolution = st.resolution

    # Get the event map and set nulls to zero
    event_map_ccp4 = gemmi.read_ccp4_map(str(event_map_path), )
    event_map_ccp4.setup(0.0)

    # FFT to reciprocal space
    sf = gemmi.transform_map_to_f_phi(event_map_ccp4.grid, half_l=True)
    data = sf.prepare_asu_data(dmin=resolution)

    # Make and output the mtz
    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = sf.spacegroup
    mtz.set_cell_for_all(sf.unit_cell)
    mtz.add_dataset('unknown')
    mtz.add_column('FWT', 'F')
    mtz.add_column('PHWT', 'P')
    mtz.set_data(data)
    mtz.write_to_file(str(output_mtz_path))


def event_maps_to_mtz(pandda_path, model_building_path):
    print(f'Making mtzs of event maps from {pandda_path} in {model_building_path}')
    pandda_path = pathlib.Path(pandda_path)
    model_building_path = pathlib.Path(model_building_path)
    for dataset_dir in (pandda_path / 'processed_datasets').glob('*'):
        print(f'Processing {dataset_dir}')
        if not dataset_dir.is_dir():
            print(f'\tNot a directory, skipping!')
            continue

        dtag = dataset_dir.name

        pdb_path = dataset_dir / f"{dtag}-pandda-input.pdb"

        if not pdb_path.exists():
            print(f'\tNo reference pdb {pdb_path}, skipping!')
            continue

        # Get the event maps
        for event_map_path in dataset_dir.glob(f'{dtag}-event_*'):
            event_mtz_path = model_building_path / f'{dtag}' / f"{event_map_path.stem}.mtz"
            print(f'\tMaking event map mtz for {event_map_path.name} at {event_mtz_path}')
            event_map_to_mtz(event_map_path, pdb_path, event_mtz_path)

    print(f'Processed!')

if __name__ == "__main__":
    fire.Fire()