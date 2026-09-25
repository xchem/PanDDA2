import time

import numpy as np
from skimage.segmentation import expand_labels
import gemmi
import torch
from torch import nn
from torch.nn import functional as F
# import lightning as lt
import pytorch_lightning as lt

from edanalyzer.datasets.event_scoring import get_ligand_array_from_lig_frame

from .interfaces import *
from .base import transform_from_arrays, SampleFrame, grid_from_template, get_ligand_mask, get_structure_array, copy_map, _get_ed_mask_float
from .constants import SAMPLE_SIZE, SAMPLE_SPACING
from .resnets import resnet10, _resnet, BasicBlock


def get_sample_frame_from_event(event: EventI, sample_size, sample_spacing) -> SampleFrameI:
    return SampleFrame(
        transform_from_arrays(
            event.centroid - (0.5*sample_spacing*np.array(sample_size)),
            np.eye(3) * sample_spacing,
        ),
        SAMPLE_SIZE,
    )


class Event(EventI):
    def __init__(self, centroid: np.array):
        self.centroid = centroid


def mask_xmap_radial(xmap: GridI, x: float, y: float, z: float, radius: float = 5.5) -> GridI:
    mask = grid_from_template(xmap)
    pos = gemmi.Position(x, y, z)
    mask.set_points_around(
        pos,
        radius=radius,
        value=1.0,
    )

    mask_arr = np.array(mask, copy=False)
    xmap_arr = np.array(xmap, copy=False)

    xmap_arr[:, :, :] = xmap_arr[:, :, :] * mask_arr[:, :, :]
    return xmap


# class LitEventScoring(lt.LightningModule):
#     def __init__(self, ):
#         super().__init__()
#         self.z_encoder = resnet10(num_classes=2, num_input=2, headless=True).float()
#         self.mol_encoder = resnet10(num_classes=2, num_input=1, headless=True).float()
#         self.fc = nn.Sequential(
#             nn.Linear(512, 2),
#         )
#
#     def forward(self, z, m, ):
#         mol_encoding = self.mol_encoder(m)
#         z_encoding = self.z_encoder(z)
#
#         full_encoding = z_encoding * F.hardtanh(mol_encoding, min_val=-1.0, max_val=1.0)
#
#         score = F.softmax(self.fc(full_encoding))
#
#
#         return float(score[0][1])


class LitEventScoring(lt.LightningModule):
    def __init__(self, config):
        super().__init__()

        self.z_encoder = _resnet(
            'resnet10',
            BasicBlock,
            [config['blocks_1'], config['blocks_2'], config['blocks_3'], config['blocks_4'], ],
            False, False,
            num_classes=3, num_input=2, headless=True, drop_rate=config['drop_rate'], config=config).float()
        self.mol_encoder = _resnet(
            'resnet10',
            BasicBlock,
            [config['blocks_1'], config['blocks_2'], config['blocks_3'], config['blocks_4'], ],
            False, False,
            num_classes=2, num_input=1, headless=True, drop_rate=config['drop_rate'], config=config).float()

        self.fc = nn.Sequential(
            nn.Linear(config['planes_5'] , 3),

        )
        self.train_annotations = []
        self.test_annotations = []

        # self.output = output_dir
        self.lr = config['lr']
        self.wd = config['wd']
        self.batch_size = config['batch_size']
        self.ligand = config['ligand']

    def forward(self, z, m, ):
        z_encoding = self.z_encoder(z)
        if self.ligand is True:
            mol_encoding = self.mol_encoder(m)
        else:
            mol_encoding = torch.ones(z_encoding.size())
        full_encoding = z_encoding * mol_encoding
        score = F.softmax(self.fc(full_encoding))

        return score


class EventScorer:

    def __init__(self, model, config, table=None, debug=False):
        self.model = model.eval().float()
        self.config = config
        self.table = table
        self.debug=debug

    def __call__(self, event: EventI, ligand_conformation: StructureI, zmap: GridI, xmap: GridI) -> float:
        # Get the sample frame
        sample_frame = get_sample_frame_from_event(
            event,
            self.config['sample_size'],
            self.config['sample_spacing']
        )

        # Cut the xmap
        # mask = _get_ed_mask_float()

        # Get the xmap sample
        if self.config['rescale_x']:
            xmap_sample = sample_frame(xmap, scale=True)
        else:
            xmap_sample = sample_frame(xmap, scale=False)

        # Get the zmap sample
        zmap_sample = sample_frame(zmap, scale=False)

        if self.config['ligand']:
            # Get the ligand mask sample
            # ligand_mask = get_ligand_mask(ligand_conformation, zmap)
            # ligand_mask_sample = sample_frame(ligand_mask, scale=False)
            ligand_frame = {
                'n': self.config['sample_size'],
                'd': self.config['sample_spacing'],
                'orientation': np.eye(3)  
            }
            ligand_mask_sample = get_ligand_array_from_lig_frame(
                ligand_conformation[0][0][0],
                ligand_frame
                )

            #
            _density_mask = (zmap_sample > self.config['z_cutoff']).astype(int)
            density_mask = expand_labels(_density_mask, distance=self.config['z_mask_radius'] / 0.5)
            density_mask[density_mask != 1] = 0
        else:
            ligand_mask_sample = np.zeros((6, self.config['sample_size'], self.config['sample_size'], self.config['sample_size']), np.float32)
            density_mask = _get_ed_mask_float(self.config['xmap_radius'])

        if self.debug:
            centroid = event.centroid
            z_mean = round(float(np.mean(zmap_sample)), 2)
            z_std = round(float(np.std(zmap_sample)), 2)
            x_mean = round(float(np.mean(zmap_sample)), 2)
            x_std = round(float(np.std(zmap_sample)), 2)
            m_num = np.sum(density_mask)
            print(f'\tEvent {round(centroid[0], 2)} {round(centroid[1], 2)} {round(centroid[2], 2)}: z: {round(z_mean, 2)} {round(z_std, 2)}: x: {round(x_mean, 2)} {round(x_std, 2)}; mask {m_num}')
        # Run the model
        map_array = np.stack(
                    [
                        zmap_sample * density_mask,
                        xmap_sample * density_mask,
                    ],
            dtype=np.float32
                )[np.newaxis,:]
        mol_array = np.stack(
                [
                    ligand_mask_sample
                ],
            dtype=np.float32
            )
        # mol_array = ligand_mask_sample

        score = self.model(
            None,
            torch.from_numpy(map_array),
            torch.from_numpy(mol_array),
            None

        ).detach().numpy()[0][1]

        print(f'Raw score: {score}')

        if self.table is not None:
            higher_scores = self.table[self.table['x']>score]['p']
            if len(higher_scores) == 0:
                rescaled_score = 1.0
            else:
                rescaled_score = higher_scores.iloc[0]
        else:
            rescaled_score = score

        return rescaled_score, map_array, mol_array
