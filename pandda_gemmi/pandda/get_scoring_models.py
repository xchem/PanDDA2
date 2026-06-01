import os
import inspect

# try:
#     from sklearnex import patch_sklearn
#
#     patch_sklearn()
# except ImportError:
#     print('No sklearn-express available!')

import gdown
import yaml

import pandas as pd

from pandda_gemmi.interfaces import *
from pandda_gemmi import serialize
from pandda_gemmi.cnn import load_model_from_checkpoint, EventScorer, LitEventScoring, BuildScorer, LitBuildScoring


def get_scoring_models(args, ):
    # Get the method for scoring events
    if args.use_ligand_data:
        event_model_path = Path(os.path.dirname(inspect.getfile(LitEventScoring))) / "model_event.ckpt"
        event_config_path = Path(os.path.dirname(inspect.getfile(LitEventScoring))) / "model_event_config.yaml"
        event_score_quantiles_path = Path(
            os.path.dirname(inspect.getfile(LitEventScoring))) / "event_score_quantiles.csv"
        if not (event_model_path.exists() & event_config_path.exists()):
            print(f'No event model at {event_model_path}. Downloading event model...')
            with open(event_model_path, 'wb') as f:
                # gdown.download('https://drive.google.com/file/d/1b58MUIJdIYyYHr-UhASVCvIWtIgrLYtV/view?usp=sharing',
                #                f)
                gdown.download(id='1b58MUIJdIYyYHr-UhASVCvIWtIgrLYtV',
                               output=f)
            with open(event_config_path, 'wb') as f:
                gdown.download(id='1qyPqPylOguzXmt6XSFaXCKrnvb8gZ8E2',
                               output=f)
            with open(event_score_quantiles_path, 'wb') as f:
                gdown.download(id='15RnkrGtEmFvtBvIlwfaUE1QfQrD2npnu', output=f)

        with open(event_config_path, 'r') as f:
            event_model_config = yaml.safe_load(f)
        score_event_model = load_model_from_checkpoint(
            event_model_path,
            LitEventScoring(event_model_config),
        ).float().eval()
        score_event = EventScorer(score_event_model, event_model_config, debug=args.debug)

        # Get the method for scoring
        build_model_path = Path(os.path.dirname(inspect.getfile(LitBuildScoring))) / "model_build.ckpt"
        build_config_path = Path(os.path.dirname(inspect.getfile(LitBuildScoring))) / "model_build_config.yaml"

        if not (build_model_path.exists() & build_config_path.exists()):
            print(f'No build model at {build_model_path}.Downloading build model...')
            with open(build_model_path, 'wb') as f:
                # gdown.download('https://drive.google.com/file/d/17ow_rxuEvi0LitMP_jTWGMSDt-FfJCkR/view?usp=sharing',
                #                f
                #                )
                gdown.download(id='17ow_rxuEvi0LitMP_jTWGMSDt-FfJCkR',
                               output=f)
            with open(build_config_path, 'wb') as f:
                gdown.download(id='1HEXHZ6kfh92lQoWBHalGbUJ-iCsOIkFo',
                               output=f)
        with open(build_config_path, 'r') as f:
            build_model_config = yaml.safe_load(f)
        score_build_model = load_model_from_checkpoint(
            build_model_path,
            LitBuildScoring(build_model_config),
        ).float().eval()
        score_build = BuildScorer(score_build_model, build_model_config)
        event_score_quantiles = pd.read_csv(event_score_quantiles_path)
    else:  # use_ligand_data=False
        event_model_path = Path(os.path.dirname(inspect.getfile(LitEventScoring))) / "model_event_no_ligand.ckpt"
        event_config_path = Path(
            os.path.dirname(inspect.getfile(LitEventScoring))) / "model_event_no_ligand_config.yaml"
        event_score_quantiles_path = Path(
            os.path.dirname(inspect.getfile(LitEventScoring))) / "event_score_no_ligand_quantiles.csv"
        if not (event_model_path.exists() & event_config_path.exists()):
            print(f'No event model at {event_model_path}. Downloading event model...')
            with open(event_model_path, 'wb') as f:
                gdown.download(id='1ccUM3g6RKluxwz8hofqmXEH2iymMvjyy',
                               output=f)
            with open(event_config_path, 'wb') as f:
                gdown.download(id='1c_QyEjFD5DtYlbU-Gh1o79gkbtdSrkDl',
                               output=f)
            with open(event_score_quantiles_path, 'wb') as f:
                gdown.download(id='1kHtBtLgGBuSBO8Mrf9pn7kjokL6fRMP6', output=f)

        with open(event_config_path, 'r') as f:
            event_model_config = yaml.safe_load(f)
        score_event_model = load_model_from_checkpoint(
            event_model_path,
            LitEventScoring(event_model_config),
        ).float().eval()
        score_event = EventScorer(score_event_model, event_model_config, debug=args.debug)
        event_score_quantiles = pd.read_csv(event_score_quantiles_path)
        score_build = None
    if args.debug:
        print(f'Using ligand?: {score_event.model.ligand} / {score_event.model.ligand is True}')
        print(f'Score model path: {event_model_path}')

    return score_event, score_build, event_model_config, event_score_quantiles
