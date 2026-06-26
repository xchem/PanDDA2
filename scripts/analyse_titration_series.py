import argparse
import yaml
from typing import TypedDict
from pathlib import Path
import time

import pandas as pd
from rich import print as rprint
import gemmi

from pandda_gemmi.args import PanDDATitrationSeriesArgs
from pandda_gemmi.interfaces import *
from pandda_gemmi.dataset import StructureArray
from pandda_gemmi.processor import ProcessLocalRay
from pandda_gemmi.pandda_logging import PanDDAConsole
from pandda_gemmi.dataset import XRayDataset

from pandda_gemmi.dmaps import (
    SparseDMapStream,
    TruncateReflections,
    SmoothReflections,
    RealSpaceSmoothReflections
)
from pandda_gemmi.fs.pandda_input import DatasetDir

from pandda_gemmi.alignment import Alignment, DFrame
from pandda_gemmi.processor import Partial
from pandda_gemmi.comparators import (
    get_comparators,
    FilterRFree,
    FilterSpaceGroup,
    FilterResolution,
    FilterCompatibleStructures,
)


class Input(TypedDict):
    dataset_dir: Path
    series: dict[str, dict[str, float]]
    output_path: Path

def parse_args(args):
    with open(args.input_yaml, 'r') as f:
        data = yaml.safe_load(f)

    return data 

def get_series(input_yaml):
    csv = pd.read_csv(input_yaml['series'])
    titration_series = {}
    for _idx, _row in csv.iterrows():
        series, crystal, concentration = _row['Series'], _row['Crystal'], _row['Concentration']
        if series not in titration_series:
            titration_series[series] = {}
        titration_series[series][crystal] = float(concentration)
    return titration_series

def get_datasets(args, input_yaml, series):
    all_dtags = [_dtag for _series_name in series for _dtag in series[_series_name]]

    dataset_dirs = {
        dataset_dir.dtag: dataset_dir
        for dataset_dir
        in [
            DatasetDir(
                path,
                None,
                args.pdb_regex,
                args.mtz_regex,
                args.ligand_dir_regex,
                args.ligand_cif_regex,
                args.ligand_smiles_regex,
                args.ligand_pdb_regex,
            )
            for path in Path(input_yaml['dataset_dir']).glob("*")
            if path.name in all_dtags
        ]
        if (dataset_dir.input_pdb_file and dataset_dir.input_mtz_file)
    }

    datasets: Dict[str, DatasetInterface] = {
        dataset_dir.dtag: XRayDataset.from_paths(
            dataset_dir.input_pdb_file,
            dataset_dir.input_mtz_file,
            dataset_dir.input_ligands,
            name=dataset_dir.dtag
        )
        for dataset_dir
        in dataset_dirs.values()
    }

    datasets_with_ligs = {}
    for dtag, dataset in datasets.items():
            lig_res = None
            for chain in dataset.structure.structure[0]:
                for res in chain:
                    if res.name == 'LIG':
                        lig_res = res
            if lig_res is not None:
                datasets_with_ligs[dtag] = dataset 


    datasets_to_process = []
    for series_name in series:
        series_datasets = [x for x in series[series_name] if x in datasets_with_ligs]
        rprint(f'{series_name} : {series_datasets}')
        if len(series_datasets) == 0:
            continue

        datasets_to_process.append(
            max(series_datasets, key=lambda _dtag: series[series_name][_dtag])
        )

    return datasets_with_ligs, datasets_to_process

# Get the sample grid
def get_sample_grid(dataset, radius=1.5):
    # Get the ligand
    lig_res = None
    for chain in dataset.structure.structure[0]:
        for res in chain:
            if res.name == 'LIG':
                lig_res = res

    # Get the heavy atom poss
    poss = []
    for atom in lig_res:
        if atom.element.name != 'H':
            pos = atom.pos
            poss.append([pos.x, pos.y, pos.z])
    pos_array = np.array(poss)

    # Define limits
    sample_min = np.min(pos_array, axis=0) - radius
    sample_max = np.max(pos_array, axis=0) + radius
    # rprint([sample_max, sample_min])
    # Scatter points
    rng = np.random.default_rng()
    initial_samples = rng.uniform(sample_min, sample_max, (10000, 3))
    # rprint('initial sample')
    # rprint(initial_samples.shape)
    # rprint(initial_samples[0])

    # Filter points
    deltas = initial_samples[:, np.newaxis, :, ] - pos_array[np.newaxis, :, :, ]
    # rprint('deltas')
    # rprint(deltas.shape)
    # rprint(deltas[0])
    distances = np.linalg.norm(deltas, axis=-1)
    # rprint('distances')
    # rprint(distances.shape)
    # rprint(distances[0])
    closest_distances = np.min(distances, axis=-1)
    # rprint('closest distances')
    # rprint(closest_distances.shape)
    # rprint(closest_distances[0])

    samples = initial_samples[closest_distances < radius]
    # rprint('samples')
    # rprint(samples.shape)
    # rprint(samples)
    return samples
    ...

def get_samples(dmaps_dict, reference_frame, sample_grid):
    samples = {}
    for dtag, sparse_xmap in dmaps_dict.items():
        mean, std = np.mean(sparse_xmap.data), np.std(sparse_xmap.data)
        sparse_xmap.data = (sparse_xmap.data - mean) / std
        xmap = reference_frame.unmask(sparse_xmap)
        # xmap_array = np.array(xmap)

        samples[dtag] = [xmap.interpolate_value(gemmi.Position(*sample_point)) for sample_point in sample_grid]

    return samples

def save_samples(samples, concentration_series, output_path):
    with open(output_path, 'w') as f:
        yaml.dump(samples, f,)
        
def plot_samples(samples, concentration_series, output_dir):
    ...

def main(args):
    # Parse the input yaml
    input_yaml = parse_args(args)
    rprint(input_yaml)

    # Get the titration series
    series = get_series(input_yaml)
    rprint('Titration series')
    rprint(series)

    # Record time at which PanDDA processing begins
    time_pandda_begin = time.time()

    # Create the console to print output throughout the programs run
    console = PanDDAConsole()

    # Print the PanDDA initialization message and the command line arguments
    console.start_pandda()

    # Get the processor to handle the dispatch of functions to multiple cores and the cache of parallel
    # processed objects
    # TODO: uses ray not mulyiprocessing_spawn
    console.start_initialise_multiprocessor()
    processor: ProcessorInterface = ProcessLocalRay(args.local_cpus)
    console.print_initialized_local_processor(args)

    # Get the datasets
    datasets, datasets_to_process = get_datasets(args, input_yaml, series)
    rprint('Datasets to process')
    rprint(datasets_to_process)
    rprint('Datasets')
    rprint(datasets)


    # Create processor references to datasets and structure arrays
    dataset_refs = {_dtag: processor.put(datasets[_dtag]) for _dtag in datasets}
    structure_array_refs = {
        _dtag: processor.put(StructureArray.from_structure(datasets[_dtag].structure)) 
        for _dtag 
        in datasets
    }

    # Process each dataset by identifying potential comparator datasets, constructing proposed statistical models,
    # calculating alignments of comparator datasets, locally aligning electron density, filtering statistical models
    # to the plausible set, evaluating those models for events, selecting a model to take forward based on those events
    # and outputing event maps, z maps and mean maps for that model
    pandda_events = {}
    autobuilds = {}

    for dtag in datasets_to_process:
        reference_series = {_dtag: _series_name  for _series_name in series for _dtag in series[_series_name]}[dtag]
        comparator_datasets: Dict[str, DatasetInterface] = {
    _dtag: datasets[_dtag] for _dtag in series[reference_series] if _dtag in series
}
        rprint(dtag)
        rprint(reference_series)
        rprint(comparator_datasets)

    time_begin_process_datasets = time.time()
    console.start_process_shells()
    for j, dtag in enumerate(datasets_to_process):
        # Record the time that dataset processing begins
        time_begin_process_dataset = time.time()

        # Get the dataset
        dataset = datasets[dtag]

        # Get the resolution of the dataset
        dataset_res = dataset.reflections.resolution()

        # Get the comparator datasets: these are filtered for reasonable data quality, space group compatability,
        # compatability of structural models and similar resolution
        reference_series = {_dtag: _series_name  for _series_name in series for _dtag in series[_series_name]}[dtag]
        comparator_datasets: Dict[str, DatasetInterface] = {
            _dtag: datasets[_dtag] for _dtag in series[reference_series] if _dtag in datasets
        }
        rprint('Comparator datasets')
        rprint(comparator_datasets)

        # Ensure the dataset itself is included in comparators
        if dtag not in comparator_datasets:
            comparator_datasets[dtag] = dataset

        # Get the resolution to process the dataset at
        processing_res = max(
            [_dataset.reflections.resolution() for _dataset in comparator_datasets.values()]
        )

        # Print basic information about the processing to be done of the dataset
        console.begin_dataset_processing(
            dtag,
            dataset,
            dataset_res,
            comparator_datasets,
            processing_res,
            j,
            datasets_to_process,
            time_begin_process_datasets
        )

        # Get the alignments, and save them to the object store
        time_begin_get_alignments = time.time()
        alignments: Dict[str, AlignmentInterface] = processor.process_dict(
            {_dtag: Partial(Alignment.from_structure_arrays).paramaterise(
                _dtag,
                structure_array_refs[_dtag],
                structure_array_refs[dtag],
            ) for _dtag in comparator_datasets}
        )
        alignment_refs = {_dtag: processor.put(alignments[_dtag]) for _dtag in comparator_datasets}
        time_finish_get_alignments = time.time()

        # Get the reference frame and save it to the object store
        time_begin_get_frame = time.time()
        reference_frame: DFrame = DFrame(dataset, processor, debug=args.debug)
        reference_frame_ref = processor.put(reference_frame)
        time_finish_get_frame = time.time()
        # TODO: Log properly
        # print(f"\t\tGot reference frame in: {round(time_finish_get_frame - time_begin_get_frame, 2)}")

        # Get the transforms to apply to the dataset before locally aligning and save them to the object store
        transforms = [
            TruncateReflections(
                comparator_datasets,
                processing_res,
            ),
            SmoothReflections(dataset)
        ]
        transforms_ref = processor.put(transforms)

        post_transforms = [
            RealSpaceSmoothReflections(dataset, fs=None, debug=args.debug)
        ]
        post_transforms_ref = processor.put(post_transforms)

        # Load the locally aligned density maps and construct an array of them
        time_begin_get_dmaps = time.time()
        # print('Datasets to transform')
        # print(sorted([_dtag for _dtag in comparator_datasets]))
        # print(sorted([_dtag for _dtag in dataset_refs]))
        dmaps_dict = processor.process_dict(
            {
                _dtag: Partial(SparseDMapStream.parallel_load).paramaterise(
                    dataset_refs[_dtag],
                    alignment_refs[_dtag],
                    transforms_ref,
                    post_transforms_ref,
                    reference_frame_ref,
                    args.debug
                )
                for _dtag
                in comparator_datasets
            }
        )
        # rprint(dmaps_dict)
        # for _dtag, _dmap in dmaps_dict.items():
        #     rprint(_dtag)
        #     rprint([np.mean(_dmap.data), np.std(_dmap.data)])

        # Get the sample grid
        sample_grid = get_sample_grid(datasets[dtag])

        # Perform sampling
        samples = get_samples(dmaps_dict, reference_frame, sample_grid)
        # rprint('samples')
        # rprint(samples)
        # rprint({_dtag: np.median(_samples) for _dtag, _samples in samples.items()})

        # Save
        save_samples(samples, series[reference_series], Path(input_yaml['output_dir']) / f'{reference_series}_samples.yaml', )
        
        # Plot
        plot_samples(samples, series[reference_series], input_yaml['output_dir'])


    rprint('Done')

        # 

if __name__ == "__main__":
    args = PanDDATitrationSeriesArgs.from_command_line()

    main(args)