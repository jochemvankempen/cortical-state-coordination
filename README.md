# cortical-state-coordination

This repository contains the code to reproduce the analyses decsribed in the paper "Top-down coordination of local cortical state during selective attention", van Kempen et al., (biorXiv) [(doi: 10.1101/2020.03.26.009365)](https://www.biorxiv.org/content/10.1101/2020.03.26.009365v1).

The data can be found [here](https://gin.g-node.org/jochemvankempen/Thiele_attention_gratc_V1_V4_laminar).

## Setup
After cloning the repository, set the path to the data folder in `./analyses/computerSetup.m`.

```matlab
paths.base = '[your local data location]/Thiele_attention_gratc_V1_V4_laminar/'
paths.server = '[your remote data location]/Thiele_attention_gratc_V1_V4_laminar/'
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