"""
Score a PanDDA result vs a set of predicted ligand poses

This requires first matching predictions to compounds to experiment results 
i.e. determining which predicted compound pose is associated with which dataset

We assert a smiles match for this - the smiles PanDDA picked up MUST be canonically 
identical to the one associated with a pose

Chemistry/Pose data is supplies in a directroy format - each dir is a chemical
and its predicted pose

Sequence alignments are determined for each chain and then the ligand is transformed 
with that chain and density tested for RSCC. The best result over all alignments/transforms
is returned. 
"""
import argparse

def score_events_vs_pose_predictions(
    pandda_dir,
    prediction_dir,
    total=True  # Assert that every prediction has at least one experiment
):
    # Get PanDDA FS model
    
    # Get PanDDA Events

    # Get prediction FS model

    # For each dataset
        # Check if any smiles match any predictions
        
        # If so, align ligands based on each chain sequence alignment transform

        # Get best RSCC over transforms

    # Format data for output

    # 
    ...

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--')

    args = parser.parse_args()

    score_events_vs_pose_predictions(
        pandda_dir,
        prediction_dir,
    )