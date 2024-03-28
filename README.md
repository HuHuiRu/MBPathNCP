# MBPathNCP

![model](./model.jpg)

## Datasets
The 'Data' folder contains the raw used in MBPathNCP.Their specific sources are detailed in the paper. The following is a brief introduction on each file:
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
consistency_projection.py：The source code for the NCP (Network Consistency Projection) algorithm.
## NCP
GIP_gate.py：The source code for the GIP (Gaussian Interaction Profile) kernel similarity algorithm.
## WKNKN
wknkn.py：The source code for the WKNKN (Weighted K Nearest Known Neighbors) algorithm. 

## Requirements
- python = 3.9
- pandas = 1.5.1
- numpy = 1.24.2

## Quick start

Run main.py to run MBPathNCP
