# Gap-filling workflow for metabolic networks with BioMeneco

This project is an integrated workflow combining the capabilities of tools like BioISO, Meneco and other relevant developed methods to address gaps in Genome-scale Metabolic models (GEMs). These gaps often arise from database limitations and incorrect genome annotations. While existing tools that identify network deficiencies or suggest reactions to fill these gaps, their efficiency decreases with increasing network complexity. This workflow aims to streamline this process, offering a more efficient and automated approach

## Installation

Create a conda environment with the following command:

```bash
conda env create -n gapfilling
```

Activate the environment with the following command:

```bash
conda activate gapfilling
```

Install requirements with the following command:

```bash
pip install -r requirements.txt
```

Install cplex with the following command, considering that the cplex folder's name is "cplex":

```bash
cd cplex && python setup.py install
```





