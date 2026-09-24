import time
import yaml

import pandas as pd

try:
    from sklearnex import patch_sklearn

    patch_sklearn()
except ImportError:
    print('No sklearn-express available!')

from pandda_gemmi.interfaces import *
from pandda_gemmi import constants
from pandda_gemmi.site_model import HeirarchicalSiteModel, HeirarchicalSiteModelAlignedSequences, ResiduePainting, Site, get_sites
from pandda_gemmi.autobuild.merge import merge_autobuilds, MergeHighestBuildScore, MergeHighestEventScore
from pandda_gemmi.ranking import rank_events, RankHighEventScoreBySite
from pandda_gemmi.tables import output_tables
from pandda_gemmi import serialize
from pandda_gemmi.event_model.event import Event
from pandda_gemmi.serialize import read_residue_assignments, read_msa, output_residue_assignments, output_msa

from pandda_gemmi.metrics import get_hit_in_site_probabilities


def postrun(
        args,
        fs,
        console,
        datasets,
        pandda_events,
        autobuilds,
        datasets_to_process,
        event_score_quantiles,
        time_pandda_begin
):
    # TODO: Log properly time taken

    # Get existing site and event data (if it exists)
    inspect_table_file = fs.output.analyses_dir / constants.PANDDA_INSPECT_EVENTS_PATH
    inspect_sites_file = fs.output.analyses_dir / constants.PANDDA_INSPECT_SITES_PATH
    msa_file = fs.output.analyses_dir / constants.PANDDA_MSA_PATH
    sequence_assignment_file = fs.output.analyses_dir / constants.PANDDA_SEQ_ASSIGN_PATH

    if (not inspect_table_file.exists()) or (msa_file.exists()):
        print(f'New PanDDA or existing msa - using site model: Residue Painting')
        site_model = ResiduePainting(
                t=0.3, 
                debug=args.debug,
                distance=10.0
                )
    else:
        print(f'Old PanDDA with no msa - using site model: Hierarchical')
        site_model = HeirarchicalSiteModelAlignedSequences(t=args.max_site_distance_cutoff, debug=args.debug)
    


    if args.site_override_file:
        with open(args.site_override_file, 'r') as f:
            site_overrides_yaml = yaml.safe_load(f)
        existing_sites = {}
        for site_idx, site_info in site_overrides_yaml.items():
            existing_sites[site_idx] = Site(
                [],
                np.zeros(3),
                dtag=site_info['dtag'],
                residues=[(chain, res) for (chain, res) in site_info['residues']]
                )
            site_overrides = {_site_idx: _site for _site_idx, _site in existing_sites.items()}
    else:
        site_overrides = None

    # Get the existing events and sites
    if inspect_table_file.exists():
        inspect_events_table = pd.read_csv(inspect_table_file)
        inspect_sites_table = pd.read_csv(inspect_sites_file)
        print(f'Found existing sites')
        print(inspect_sites_table)
        existing_events = {
            (_row['dtag'], _row['event_idx']): Event(
                np.array([_row['x'], _row['y'], _row['z']]),
                None,
                0,
                np.array([_row['x'], _row['y'], _row['z']]),
                score=_row['z_peak'],
                site_idx = _row['site_idx']
            )
            for _idx, _row
            in inspect_events_table.iterrows()
        }
        if not existing_sites:
            existing_sites = {}
        existing_sites.update(
            {
                _row['site_idx']: Site(
                    [event_id for event_id, event in existing_events.items() if event.site_idx == _row['site_idx']],
                    _row['centroid'],
                    _row['Name'],
                    _row['Comment']
                ) 
                for _idx, _row 
                in inspect_sites_table.iterrows()
            }
        )

    else:
        print(f'Found no existing PanDDA Results at {inspect_table_file}')
        existing_events = None
        existing_sites = None

    # Get existing site residues and sequence alignments
    if msa_file.exists() & sequence_assignment_file.exists():
        residue_assignments = read_residue_assignments(sequence_assignment_file)
        msa = read_msa(msa_file)

        if existing_sites is not None:
            site_id_to_residues = {v: [_k for _k in residue_assignments if residue_assignments[_k] == v] for k, v in residue_assignments.items()}

            for site_id, site in sites.items():
                site.dtag = site_id_to_residues[site_id][0][0]
                site.residues = [(x[1], x[2]) for x in site_id_to_residues[site_id]]
    else:
        msa = None

    # Autobuild the best scoring event for each dataset
    console.start_autobuilding()

    autobuild_yaml_path = fs.output.path / "autobuilds.yaml"
    if autobuild_yaml_path.exists():
        autobuilds = serialize.unserialize_autobuilds(autobuild_yaml_path)
    else:
        # Merge the autobuilds into PanDDA output models
        if args.use_ligand_data & args.autobuild:
            merged_build_scores = merge_autobuilds(
                datasets,
                pandda_events,
                autobuilds,
                fs,
                MergeHighestEventScore()
            )

        #
        console.processed_autobuilds(autobuilds)

    print(f'Event scores: {[event.score for event_id, event in pandda_events.items()]}')

    # Get the sites
    sites, residue_assignments, msa = get_sites(
        datasets,
        pandda_events,
        datasets_to_process[
            min(
                datasets_to_process,
                key=lambda _dtag: datasets_to_process[_dtag].reflections.resolution()
            )
        ],
        site_model,
            existing_events,
            existing_sites,
            site_overrides,
            msa
    )
    output_residue_assignments(residue_assignments, sequence_assignment_file)
    output_msa(msa, msa_file)

    # TODO: Log properly sites

    # Rank the events for display in PanDDA inspect
    ranking, sorted_sites = rank_events(
        pandda_events,
        sites,
        autobuilds,
        RankHighEventScoreBySite(),
        existing_events
    )

    # Probabilities
    # Calculate the cumulative probability that a hit remains in the site using the event score quantile table
    hit_in_site_probabilities = get_hit_in_site_probabilities(pandda_events, ranking, sorted_sites, event_score_quantiles)

    # Output the event and site tables
    output_tables(datasets, pandda_events, ranking, sorted_sites, hit_in_site_probabilities, fs, existing_events, existing_sites, args.debug)
    time_pandda_finish = time.time()
    # TODO: Log properly pandda run time
    print(f"PanDDA ran in: {round(time_pandda_finish - time_pandda_begin, 2)} seconds!")
