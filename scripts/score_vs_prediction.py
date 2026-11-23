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