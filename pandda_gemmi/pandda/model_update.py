import os
import inspect
from pathlib import Path

import gdown
import yaml
from zenodo_get import download


from pandda_gemmi import constants
from edanalyzer.models.event_scoring import LitEventScoring

def update_models():
    # Get model version from constants
    print(Path(os.path.dirname(inspect.getfile(LitEventScoring))) / "model_event_config.yaml")

    # Get model version from files
    configs = {
        "event": {
            "path": Path(os.path.dirname(inspect.getfile(LitEventScoring))) / "model_event_config.yaml", 
            "version": constants.EVENT_MODEL_VERSION,
            "files": [
                "model_event_config.yaml", 
                "event_score_quantiles.csv", 
                "model_event.ckpt",
                ]
                  },
        "event_no_ligand": {
            "path": Path(os.path.dirname(inspect.getfile(LitEventScoring))) / "model_event_no_ligand_config.yaml", 
            "version": constants.EVENT_MODEL_NO_LIGAND_VERSION,
            "files": [
                "event_score_no_ligand_config.yaml", 
                "event_score_no_ligand_quantiles.csv", 
                "model_event_no_ligand.ckpt",
                ]
            },
        "build": {
            "path": Path(os.path.dirname(inspect.getfile(LitEventScoring))) / "model_build_config.yaml", 
            "version": constants.BUILD_MODEL_VERSION,
            "files": [
                "model_build_config.yaml", 
                "model_build.ckpt",
            ]
            },
    }

    # If there are any differences (including missing files), download from zenodo
    for model_type, model in configs.items():
        path, version, files = model['path'], model['version'], model['files']
        do_update = False
        # Check if missing
        if not path.exists():
            do_update = True    
            print(f'No config file for model: {model_type} at {path}! Downloading new model!')
        else:
            with open(path, 'r') as f:
                model_config = yaml.safe_load(f)
            if model_config['version'] != version:
                print(f'Model version mismatch: {model_type} - code - {version} vs file - {model_config["version"]}! Downloading new model!')
                do_update = True

        if do_update:
            out_dir = Path(os.path.dirname(inspect.getfile(LitEventScoring)))
            for file in files:
                if (out_dir / file).exists():
                    os.remove(file)
            download(
                record_or_doi=constants.MODEL_ZENODO_DOI,
                output_dir=out_dir,
                file_glob=files,
            )
            print(f'Downloaded new files for from zendodo for {files}')

        

    ...