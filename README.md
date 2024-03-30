# MBPathNCP

The MBPathNCP is designed for the prediction of metabolic pathways of chemicals and enzymes. The kernels for chemicals/enzymes and pathways were constructed using the interactions of chemicals and proteins, and the validated associations between chemicals/enzymes and pathways. The network consistency projection was applied to the kernels and association adjacency matrix to yield the association score for each pair of chemicals/enzymes and pathways.
![model](./model.jpg)

## Datasets
The 'Data' folder contains the raw data used in MBPathNCP. Their specific sources are detailed in the paper. The following is a brief introduction on each file:
- Pathways/Path_CID_matrix.csv：The association between 2329 chemicals and 176 metabolic pathways.
- Pathways/Path_CPC_matrix.csv：The association between 2329 chemicals,1098 enzymes and 176 metabolic pathways.
- Pathways/Path_ENSP_matrix.csv：The association between 1098 enzymes and 176 metabolic pathways.
- chemical_enzyme_interactions/DataProcess/CCI.csv：The confidence scores for interactions between 2329 chemicals.
- chemical_enzyme_interactions/DataProcess/CPI.csv：The confidence scores for interactions between 2329 chemicals and 1098 enzymes.
- chemical_enzyme_interactions/DataProcess/PPI.csv：The confidence scores for interactions between 1098 enzymes.
- chemical_enzyme_interactions/1098PPI_matrix.csv: The association between 1098 enzymes.
- chemical_enzyme_interactions/2329CCI_matrix.csv: The association between 2329 chemicals.
- chemical_enzyme_interactions/3427CPC_net.csv: The association between 2329 chemicals and 1098 enzymes.

## GIP
GIP_gate.py：The source code for the GIP (Gaussian Interaction Profile) kernel similarity algorithm.
## NCP
consistency_projection.py：The source code for the NCP (Network Consistency Projection) algorithm.
## WKNKN
wknkn.py：The source code for the WKNKN (Weighted K Nearest Known Neighbors) algorithm. 

## Requirements
- python = 3.9.19
- pandas = 2.2.1
- numpy = 1.26.4
- scikit-learn = 1.4.1
## Quick start

Run main.py to run MBPathNCP
