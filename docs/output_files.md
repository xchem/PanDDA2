# PanDDA 2 Output File Description

## {output_dir}/input.yaml

A yaml file of the input arguments to the PanDDA

## {output_dir}/processed_datasets

A directory containing the datasets processed by the PanDDA

### {output_dir}/processed_datasets/{dtag}

A directory containing the results of processing the dataset {dtag}

#### {output_dir}/processed_datasets/{dtag}/processed_dataset.yaml

A yaml contianining various information about the processing of dataset {dtag}

#### {output_dir}/processed_datasets/{dtag}/events.yaml

A yaml containing information about the PanDDA Events in datasets {dtag}

#### {output_dir}/processed_datasets/{dtag}/{dtag}-pandda-input.pdb

The input ground state pdb for dataset {dtag}

#### {output_dir}/processed_datasets/{dtag}/{dtag}-pandda-input.mtz

The input reflection data mtz for dataset {dtag}

#### {output_dir}/processed_datasets/{dtag}/{dtag}-ground-state-average-map.native.ccp4

The ground state that generated the best scoring event for dataset {dtag} 

#### {output_dir}/processed_datasets/{dtag}/{dtag}-z_map.native.ccp4

The z-map that generated the best scoring event for dataset {dtag}

#### {output_dir}/processed_datasets/{dtag}/xmap.ccp4

The post-processed 2Fo-Fc map used by the best scoring event for dataset {dtag}

#### {output_dir}/processed_datasets/{dtag}/{dtag}-event_{event_num}_1-BDC_{bdc}_map.native.ccp4

A background corrected PanDDA event map for dataset {dtag}

#### {output_dir}/processed_datasets/{dtag}/ligand_files

The directory containing the ligand information for dataset {dtag}

#### {output_dir}/processed_datasets/{dtag}/modelled_structures

A directory containing user built models for dataset {dtag}

##### {output_dir}/processed_datasets/{dtag}/modelled_structures/{dtag}-pandda-model.pdb

The current user model for dataset {dtag}

## {output_dir}/analyses

The directory containing pan-dataset summary tables used in `pandda.inspect`

### {output_dir}/analyses/pandda_analyse_events.csv

The output table of PanDDA events used in `pandda.inspect`

### {output_dir}/analyses/pandda_analyse_sites.csv

The output table of PanDDA sites used in `pandda.inspect`

### {output_dir}/analyses/pandda_inspect_events.csv

The working table of PanDDA events used in `pandda.inspect`

### {output_dir}/analyses/pandda_inspect_sites.csv

The working table of PanDDA sites used in `pandda.inspect`