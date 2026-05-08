import time
from pathlib import Path

from rich import print as rprint
import gemmi
import numpy as np

from pandda_gemmi.site_model import HeirarchicalSiteModel, HeirarchicalSiteModelAlignedSequences, Site, get_sites
from pandda_gemmi.event_model.event import Event
from pandda_gemmi.fs import PanDDAInput
from pandda_gemmi.dataset import XRayDataset

def test_HeirarchicalSiteModelAlignedSequences():
    data_dirs = Path('data/XX01ZVNS2B')
    existing_events = {}
    existing_sites = {}

    # Get the datasets
    fs = PanDDAInput(data_dirs, data_dirs, pdb_regex, mtz_regex)
    datasets = {
        dataset_dir.dtag: XRayDataset.from_paths(
            dataset_dir.input_pdb_file,
            dataset_dir.input_mtz_file,
            dataset_dir.input_ligands,
            name=dataset_dir.dtag
        )
        for dataset_dir
        in fs.dataset_dirs.values()
    }

    # Make synthetic events
    event_data = {
        'XX01ZVNS2B-x10175': [12.5, 16.0, 10.5],
        'XX01ZVNS2B-x11441': [12.0, 16.0, 12.0],
        'XX01ZVNS2B-x12469': [14.5, 15.0, 12.5],
        'XX01ZVNS2B-x10225': [12.5, 16.0, 10.0],
        'XX01ZVNS2B-x10978': [-11.0, 11.0, 19.5],
        'XX01ZVNS2B-x10754': [-9.5, 6.5, 17.5],
        'XX01ZVNS2B-x10645': [20.5, 27.5, 24.5],
    }
    pandda_events = {
        (dtag, 1): Event(
            None,
            None,
            0,
            np.array(pos),
        )
        for dtag, pos 
        in event_data.items()
    }

    # Get ref
    ref_dataset = datasets[
            min(
                datasets,
                key=lambda _dtag: datasets[_dtag].reflections.resolution()
            )
        ]

    # Get the sites
    sites = get_sites(
        datasets,
        pandda_events,
        ref_dataset,
        HeirarchicalSiteModelAlignedSequences(
            t=10.0, 
            debug=True,
            ),
        existing_events,
        existing_sites
    )