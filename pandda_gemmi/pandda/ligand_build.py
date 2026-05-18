import pathlib
import os
import dataclasses

import yaml
from rich import print as rprint
import gemmi 

from pandda_gemmi import constants 
from pandda_gemmi.interfaces import *
from pandda_gemmi.args import PanDDAArgs
from pandda_gemmi.pandda.pandda import pandda
from pandda_gemmi.autobuild.inbuilt import mask_dmap, get_conformers, score_conformer
from pandda_gemmi.fs.pandda_input import parse_dir_ligands
from pandda_gemmi.pandda.get_scoring_models import get_scoring_models
from pandda_gemmi.dataset.structure import save_structure, Structure
from pandda_gemmi.dataset.reflections import Reflections

@dataclasses.dataclass
class GetScoringModelsArgs:
    debug: bool
    use_ligand_data: bool

def read_pandda_map(xmap_file):
    """PanDDA 2 maps are often truncated, and PanDDA 1 maps can have misasigned spacegroups. 
    This method handles both."""
    dmap_ccp4 = gemmi.read_ccp4_map(str(xmap_file), setup=False)
    dmap_ccp4.grid.spacegroup = gemmi.find_spacegroup_by_name('P1')
    dmap_ccp4.setup(0.0)
    dmap = dmap_ccp4.grid 
    return dmap

def read_mtz(path):
    return Reflections.from_path(path).transform_f_phi_to_map()

def get_ligand_conformers(ligand_files):
    conformers = {}
    for ligand_key in ligand_files:
        ligand_files = ligand_files[ligand_key]
        conformers[ligand_key] = get_conformers(ligand_files)

    return conformers

def get_event_info(dataset_dir):

    # Determine which builds to perform. More than one binder is unlikely and score ranks 
    # well so build the best scoring event of each dataset.     
    event_yaml_path = dataset_dir / 'events.yaml'
    print(f'Getting event info from: {event_yaml_path}')
    with open(event_yaml_path, 'r') as f:
        events_info = yaml.safe_load(f)
    
    event_scores = {}
    for event_id, event_info in events_info.items():
        event_scores[event_id] = event_info['Score']

    best_event_id = max(event_scores, key=lambda _x: event_scores[_x])

    return events_info[best_event_id], best_event_id



def main(dataset_dir):
    rprint('# PanDDA Dir')
    rprint(dataset_dir)
    dataset_dir = pathlib.Path(dataset_dir).resolve()
    dtag = dataset_dir.name

    out_dir = dataset_dir / constants.PANDDA_MODELLED_STRUCTURES_DIR / 'pandda_autobuild'
    out_dir.mkdir(exist_ok=True)
    rprint(out_dir)

    # Get the event info for the dataset
    rprint('# Getting event info...')
    event_info, best_event_id = get_event_info(dataset_dir)
    rprint(event_info)
    rprint(best_event_id)

    # Get dataset ligand info
    rprint('# Getting ligand data...')
    ligand_files = parse_dir_ligands(
        dataset_dir / 'ligand_files', 
        constants.ARGS_LIGAND_SMILES_REGEX_DEFAULT, 
        constants.ARGS_LIGAND_CIF_REGEX_DEFAULT, 
        constants.ARGS_LIGAND_PDB_REGEX_DEFAULT,
        check_input=False
        )
    rprint(ligand_files)

    # Generate conformers of the dataset ligands to score
    rprint('# Getting ligand conformers...')
    conformers = get_ligand_conformers(ligand_files)
    rprint(conformers)

    # Get the event map grid
    rprint('# Getting event map grid...')
    event_map_grid = read_pandda_map(
        dataset_dir / constants.PANDDA_EVENT_MAP_FILE.format(
            dtag=dtag, 
            event_idx=best_event_id, 
            bdc=round(1-event_info['BDC'], 2),
            )
            )
    rprint(event_map_grid)

    # Get the z grid
    rprint('# Getting z map grid...')
    z_grid = read_pandda_map(dataset_dir / constants.PANDDA_Z_MAP_FILE.format(dtag=dtag))
    rprint(z_grid)

    # Get the raw xmap
    rprint('# Getting raw xmap...')
    raw_xmap_grid = read_mtz(dataset_dir / constants.PANDDA_MTZ_FILE.format(dtag))
    rprint(raw_xmap_grid)


    # Get the scoreing model
    rprint('# Getting build scoring model...')
    _, score_build, _, _ = get_scoring_models(GetScoringModelsArgs(False, True))

    # Get the event centroid
    rprint('# Getting event centroid...')
    centroid = event_info['Centroid']
    rprint(centroid)

    # Build
    rprint('# Generating builds...')
    builds, scores = {}, {}
    for ligand_key, ligand_conformers in conformers.items():
        for conformer_key, conformer in ligand_conformers.items():
            optimized_structure, score, centroid, arr = score_conformer(
                centroid,
                conformer,
                event_map_grid,
                score_build,
                z_grid,
                raw_xmap_grid,
            )
            builds[(ligand_key, conformer_key)] = optimized_structure
            scores[(ligand_key, conformer_key)] = score
    rprint(builds)
    rprint(scores)

    # Select the best structure
    rprint('# Getting best build...')
    best_build_id = max(scores, key=lambda _x: scores[_x])
    rprint(best_build_id)

    # Save the best structure
    rprint('# Saving structure...')
    save_structure(
        Structure(None, builds[best_build_id]),
        out_dir / f"{best_build_id[0]}_{best_build_id[1]}.pdb",
    )