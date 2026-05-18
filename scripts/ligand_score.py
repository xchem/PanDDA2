import argparse

import gemmi 
import numpy as np

from pandda_gemmi.pandda.get_scoring_models import get_scoring_models


def read_xmap(path): 
    xmap_ccp4 = gemmi.read_ccp4_map(str(path), )
    xmap_ccp4.setup(0.0)
    xmap = xmap_ccp4.grid

    return xmap


def get_ligand_structure(structure, ligand_id):
    new_structure = structure.clone()
    chains_to_delete = []
    for model in new_structure:
        for chain in model:
            chains_to_delete.append((str(model.num), chain.name))

    for model_name, chain_name in chains_to_delete:
        del new_structure[int(model_name)-1][chain_name]

    for model in structure:
        new_model = gemmi.Model(1)
        for chain in model:
            if chain.name == ligand_id[0]:
                new_chain = gemmi.Chain(chain.name)
                for res in chain:
                    if str(res.seqid.num) == ligand_id[1]:
                        new_chain.add_residue(res)
                new_model.add_chain(new_chain)
        new_structure.add_model(new_model)
    return new_structure


def read_mtz(path, f, phi, template):
    mtz = gemmi.read_mtz_file(path)
    grid = mtz.transform_f_phi_to_map(
        f,  
        phi, 
        exact_size=[template.nu, template.nv, template.nw]
        )

    grid_np = np.array(grid, copy=False)
    z_np = np.array(template, copy=False)

    grid_np[z_np == 0.0] = 0.0

    return grid


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--structure_path',)
    parser.add_argument('--ligand_id')
    parser.add_argument('--zmap_path')
    parser.add_argument('--mtz_path',)
    parser.add_argument('--out_path',)    
    parser.add_argument('--f', default='FWT')
    parser.add_argument('--phi', default='PHWT')
    parser.add_argument('--use_ligand_data', default=True)
    parser.add_argument('--debug', default=False)
    args = parser.parse_args()

    # Get the event and build scores
    print(f'Getting scoring model...')
    _, score_build, _, _ = get_scoring_models(args)

    # Get the structure, pulling out everything but the ligand
    print(f'Getting ligand structure...')
    ligand_id = args.ligand_id.split('/')
    print(ligand_id)
    structure = gemmi.read_structure(args.structure_path)
    ligand_structure = get_ligand_structure(structure, ligand_id)
    print(ligand_structure)
    print(ligand_structure[0])
    print([chain.name for chain in ligand_structure[0]])
    print([res.seqid for chain in ligand_structure[0] for res in chain])

    # Get the z grid
    print(f'Geting zmap...')
    z_grid = read_xmap(args.zmap_path)

    # Get the xmap grid
    print(f'Getting xmap...')
    raw_xmap_grid = read_mtz(
        args.mtz_path, 
        args.f, 
        args.phi, 
        z_grid,
        )

    # Score
    print(f'Scoring...')
    score, arr = score_build(
        ligand_structure,
        z_grid,
        raw_xmap_grid,
    )
    print(score)

    # Write score
    print(f'Writing score...')
    with open(args.out_path, 'w') as f:
        f.write(str(score[0][0]))