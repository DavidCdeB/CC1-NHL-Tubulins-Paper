# Relevant scripts and all workflow instructions necessary to reproduce the results presented in _Two members of a conserved LEA_2-domain protein family direct root growth under salt stress_


This repository hosts the relevant scripts and all workflow instructions necessary to reproduce the results presented here:

Title 
- Two members of a conserved LEA_2-domain protein family direct root growth under salt stress 
- LEA_2 proteins steer root growth in salt stress. 

Authors:
Orianne Montulet1, David Carrasco-Busturia, Jordy Perez Gonzalez, Cheryl Quek Chai Hui5, Liu Wang, Anja Fröhlich1, Ghazanfar Abbas Khan4, René Schneider5, Patrick Norman, Staffan Persson,2,3* Arun Sampathkumar1*

## Table of Contents

- [Citation](#citation)
- [CC1-Tubulin Dimer](#cc1-tubulin-dimer)
- [NHL28-Tubulin Dimer](#nhl28-tubulin-dimer)
- [NHL30-Tubulin Dimer](#nhl30-tubulin-dimer)
- [All regions map](#all-regions-map)
- [CC1-Tubulin Lattice](#cc1-tubulin-lattice)

## Citation

If you use this code or our methods in your research, please cite our paper:

```bibtex
@article{Montulet_et_al,
}
```

The relevant molecular dynamics trajectories necessary to reprduce the findings in this paper can be found here:

[![DOI](https://zenodo.org/badge/DOI/xyz/zenodo.xyz.svg)](https://)

<!-- https://doi.org/10.5281/zenodo.17338213 -->






## CC1-Tubulin Dimer

The FASTA sequences of $\alpha$-tubulin and $\beta$-tubulin were extracted from the sT2R-vinblastine structure (PDB ID: 4EB6).
[α- and β-tubulin sequences](CC1-Tubulin-Dimer/tubulin-A-B.fasta)

The FASTA sequence for the CC1 protein was extracted from the AF-Q9C636-F1 model, with only the
first 120 residues included.
[CC1 sequence](CC1-Tubulin-Dimer/CC1-res_1-to-res_120.fasta)

The three structures were run in AlphaFold Server (AlphaFold 3). The resulting structure is:

- [Predicted structure (.cif)](CC1-Tubulin-Dimer/fold_a_b_cc1_2025_10_27_16_55_model_0.cif)
<!-- - [VIA file (.via)](CC1-Tubulin-Dimer/fold_a_b_cc1_2025_10_27_16_55_model_0.via) -->

The .cif file was then converted to .pdb

```
> obabel -icif  fold_a_b_cc1_2025_10_27_16_55_model_0.cif -opdb -O fold_a_b_cc1_2025_10_27_16_55_model_0.pdb
```

and then proceed with the system preparation as described in the paper: CHARMM-GUI, MM energy minimization, followed by MM-MD equilibration and MM-MD production.

Available files in the Zenodo repository:

```text
CC1-Tubulin-Dimer/
└── Prod_NPT_10ns
    ├── imgui.ini
    ├── index.ndx
    ├── step5_production.cpt
    ├── step5_production.gro
    ├── step5_production.mdp
    ├── step5_production.tpr
    ├── step5_production.via
    ├── step5_production.xtc
    ├── step5_production.xtc.cache
    ├── topol.top
    └── toppar
        ├── CLA.itp
        ├── forcefield.itp
        ├── PROA.itp
        ├── PROB.itp
        ├── PROC.itp
        ├── SOD.itp
        └── TIP3.itp
```



## NHL28-Tubulin Dimer

The FASTA sequence for the NHL28 is here
[NHL28 sequence](NHL28-Tubulin-Dimer/NHL28.fasta)

NHL28,  $\alpha$-tubulin and $\beta$-tubulin were run in AlphaFold Server (AlphaFold 3). The resulting structure is:

- [Predicted structure (.cif)](NHL28-Tubulin-Dimer/fold_a_b_emb_model_0.cif)

Available files:

```text
NHL28-Tubulin-Dimer/
└── Prod_NPT_10ns
    ├── index.ndx
    ├── step5_production.cpt
    ├── step5_production.gro
    ├── step5_production_10ns.gro
    ├── step5_production_10ns.mdp
    ├── step5_production_10ns.tpr
    ├── step5_production_10ns.via
    ├── step5_production_10ns.xtc
    ├── step5_production_10ns_via_from_dardel.via
    ├── topol.top
    └── toppar
        ├── CLA.itp
        ├── forcefield.itp
        ├── PROA.itp
        ├── PROB.itp
        ├── PROC.itp
        ├── PROD.itp
        ├── SOD.itp
        └── TIP3.itp
```

## NHL30-Tubulin Dimer

The FASTA sequence for NHL30 is here:
[NHL30 sequence](NHL30-Tubulin-Dimer/NHL30.fasta)

NHL30, $\alpha$-tubulin and $\beta$-tubulin were run in AlphaFold Server (AlphaFold 3). The resulting structure is:

- [Predicted structure (.cif)](NHL30-Tubulin-Dimer/fold_a_b_lea_model_0.cif)
<!-- - [VIA file (.via)](NHL30-Tubulin-Dimer/fold_a_b_lea_model_0.via) -->

From the `NHL30-Tubulin-Dimer` directory, convert the `.cif` file to `.pdb`:

```bash
obabel -icif fold_a_b_lea_model_0.cif -opdb -O fold_a_b_lea_model_0.pdb
```

This creates the PDB file referenced by the model VIA file. Then proceed with the system preparation as described in the paper: CHARMM-GUI, MM energy minimization, followed by MM-MD equilibration and MM-MD production.

Available files in the Zenodo repository:

```text
NHL30-Tubulin-Dimer/
└── Prod_NPT_10ns
    ├── imgui.ini
    ├── index.ndx
    ├── step5_production.cpt
    ├── step5_production.gro
    ├── step5_production_10ns.gro
    ├── step5_production_10ns.mdp
    ├── step5_production_10ns.tpr
    ├── step5_production_10ns.via
    ├── step5_production_10ns.xtc
    ├── step5_production_10ns.xtc.cache
    ├── step5_production_10ns_via_from_dardel.via
    ├── topol.top
    └── toppar
        ├── CLA.itp
        ├── forcefield.itp
        ├── PROA.itp
        ├── PROB.itp
        ├── PROC.itp
        ├── SOD.itp
        └── TIP3.itp
```



## All regions map

These files help map the binding regions for CC1, NHL28 and NHL30 into the  αβ-tubulin dimer (Fig. 5D in the main manuscript):

```text
All_regions/
├── step3_input.gro
└── tubulin_binding_regions_all.via
```


## CC1-Tubulin Lattice

To place CC1 or NHL28 on the convex side of the tubulin lattice, follow the [placement procedure in steps-to-reproduce-EMB-on-Tubulin-Lattice directory](steps-to-reproduce-EMB-on-Tubulin-Lattice/README.md). All scripts used in this procedure are available in the [steps-to-reproduce-EMB-on-Tubulin-Lattice directory](steps-to-reproduce-EMB-on-Tubulin-Lattice/).


We performed a series of pulling runs 1-10 producing Fig. S7 shown in the S.I. This figure was produced using this script: [plot-3.py](CC1-Tubulin-Lattice/plot-3.py).

Across the pulling runs, we went from run 7 to run 8 by extracting the last frame of run 7 like this:

```
> gmx check -f pull_7.xtc
Item        #frames Timestep (ps)
Step           142    10
Time           142    10
Lambda           0
Coords         142    10
Velocities       0
Forces           0
Box            142    10


t_last = (142 - 1) * 10 ps = 1410 ps

Extract last frame:

> gmx trjconv -s pull_7.tpr -f pull_7.xtc -o lastframe.gro -dump 1410

> gmx grompp -f pulling_8.mdp  -c lastframe.gro  -t pull_7.cpt  -p topol.top -n index_pull_all_and_backbones.ndx -o pull_8.tpr

```


Available files in the Zenodo repository:

```text
CC1-Tubulin-Lattice/
└── Pull_Run_9
    ├── index_pull_all_and_backbones.ndx
    ├── lastframe.gro
    ├── pull.via
    ├── pull_8.cpt
    ├── pull_9.tpr
    ├── pull_9.xtc
    ├── pulling_8.mdp
    ├── topol.top
    └── toppar
        ├── CLA.itp
        ├── forcefield.itp
        ├── PROA.itp
        ├── PROB.itp
        ├── PROC.itp
        ├── PROD.itp
        ├── PROE.itp
        ├── PROF.itp
        ├── PROG.itp
        ├── PROH.itp
        ├── PROI.itp
        ├── PROJ.itp
        ├── PROK.itp
        ├── PROL.itp
        ├── PROM.itp
        ├── PRON.itp
        ├── PROO.itp
        ├── PROP.itp
        ├── PROQ.itp
        ├── PROR.itp
        ├── PROS.itp
        ├── SOD.itp
        └── TIP3.itp
```
