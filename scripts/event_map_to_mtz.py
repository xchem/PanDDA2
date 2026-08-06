import gemmi
import fire


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


if __name__ == "__main__":
    fire.Fire(event_map_to_mtz)