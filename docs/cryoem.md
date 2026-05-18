# PanDDA 2 For Cryo-Em

## Installing Dependencies

As this is an experimental feature additional dependencies are required:

```
python -m pip install mrcfile
```

## Preparing Data

Data should be arranged in the typical PanDDA 2 input format:

```
dataset-1
    - compound
        - <ligand_id>.cif
    - dataset-1.pdb
    - dataset-1.mrc
...
dataset-n
    ...
```

Once the dataset is prepared the `mrc` files can be converted to `mtzs` with the script:

```
python pandda2/scripts/prepare_pandda_from_cryoem.py <path to input data directory>
```

Which will update the input directory to be PanDDA 2 compatible:


```
dataset-1
    - compound
        - <ligand_id>.cif
    - dataset-1.pdb
    - dataset-1.mrc
    - dataset-1.mtz  # NEW!
...
dataset-n
    
```

## Running the PanDDA

PanDDA 2 needs some custom parameters to function reasonably on Cryo-EM data, and should be run like so:

```
python scripts/pandda.py --data_dirs=<prepared dataset dir> --out_dir=<desired output dir> --min_characterisation_datasets=3 --high_res_buffer=0.5
```

These two parameters control how characterization sets are selected. PanDDA will not analyze datasets for which credible ground states cannot be constructed (as defined by these parameters), and at the small sample sizes in cryoem they need to be adjusted to allow the program to run.
 - `--min_characterisation_datasets=3`: Require at least three datasets to characterize the ground state for a dataset. Normally this is much higher, with the default being `25`, however with small samlple sizes this constraint prevents any datasets being analyzed. 
 - `--high_res_buffer=0.5`: The maximum extent to which an analyzed dataset will have its analysis resolution dropped in order to get more comparators. The default is `0.1`, which can prevent the highest resolution datasets in a small sample from getting enough comparators to be analyzed.