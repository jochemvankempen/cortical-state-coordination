# cortical-state-coordination

This repository contains the code to reproduce the analyses decsribed in the paper "Top-down coordination of local cortical state during selective attention", van Kempen et al., (Neuron, 2020) [(doi: 10.1016/j.neuron.2020.12.013)](https://authors.elsevier.com/sd/article/S0896-6273(20)30995-8). 

The hidden Markov model (HMM) was originally described in "Selective modulation of cortical state during spatial attention", Engel et al., (Science, 2016) [(doi: 10.1126/science.aag1420)](https://science.sciencemag.org/content/354/6316/1140.full)

The data can be found [here](https://doi.org/10.12751/g-node.b0mnn2).

## Setup
After cloning the repository, set the path to the data folder in `./analyses/computer_setup.m`.

```matlab
paths.base ='[your-local-drive]\Thiele_attention_gratc_V1_V4_laminar\';
paths.server = '[your-remote-drive]\Thiele_attention_gratc_V1_V4_laminar\';
```

## Folders/components

### analyses
A collection of analyses. Analyses are run from `./analyses/batch_analyses.m`.

### hpc
A collection of bash/shell scripts to run analyses on a high-performance-computing cluster (SLURM). 

### lists
Lists with recording information and matlab scripts to load this info. 
- Thiele.Session_[subjectName]: list of recordings
- Thiele.PenGrid_[subjectName]: list of which areas (and accompanying channel indices) was recorded from for given session.

## Credits
Jochem van Kempen  
[Tatiana A Engel](https://www.cshl.edu/research/faculty-staff/tatiana-engel/#research-profile)

## Licence
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
