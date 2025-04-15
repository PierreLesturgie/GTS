# GTS - Pairwise Gene Tree Statistics Tools
---
#### Author: Pierre Lesturgie (pierrelesturgie@outlook.fr)

### (1) Installation: 

#### 1. Create env conda (optional): 

```conda create --name GTS python=3.12```

```conda activate GTS```

#### 2. Install: 

```pip install . ```

###### Uninstall : pip uninstall gts

###### Note: if error occurs related to msprime: uninstall and reinstall using pip

### (2) Running: 

```gts <command> <arguments>``` (do ```gts --help``` for more information)
 
#### 1. STATS. computing summary and topology statistics for each tree in the .trees tree sequence

```gts stats -t <treefile> -Nanc <ancestral_effective_size> -r <number_of_sampling_runs> -R <number_of_runs> -D <demographic_scenario> -G <genomic_structure>```

##### Arguments: 
  -t, --treefile: Path to file with the tree sequence (omit '.trees')
  
  -Nanc, --ancestral_effective_size: Ancestral effective size (used for recapitation)
  
  -r, --runs: Number of sampling runs
  
  -R, --rho: Recombination rate (used for recapitation)
  
  -D, --demographic_scenario: Demographic scenario file
  
  -G, --genomic_map: Genomic map

##### Output: 
- **Topology** (out four kind of topologies; see appendix)
- **Coalescent times** (C4, C3, C2)
- **TMRCA**
- **Pi(s)**
- **1D-SFS(s)**
- **Total Branch Length**
- **DAFi** (only multi-demes)
- **dxy** (only multi-demes)
- **Hudson's FST** (multi-demes)
- **Distribution of TMRCAs** per pair of individuals (multi-demes)
- **2D-SFS** (multi-demes)

###### Note: expressed in units of branch lentgths (when applicable)

#### 2. SUMMARY. Summarize output stats from STATS. 

```gts summary --tag <treefile> -r <number_of_sampling_runs> -G <genomic_structure>```

##### Arguments: 

  -t, --treefile: Path to file with the tree sequence (omit '.trees')
  
  -r, --runs:  Number of runs ** Must be the same than in ```gts stats```
  
  -G, --genomic_map: Genomic map
  
##### Output: 
- Proportion of each **topology** per region
- Average **coalescent times** (C4, C3, C2) for each topology
- Average **TMRCA**
- Average **Pi(s)**
- Average **1D-SFS(s)**
- Average **Total Branch Length**
- Average **DAFi** (only multi-demes)
- Average **dxy** (only multi-demes)
- Average **Hudson's FST** (multi-demes)
- Average **Distribution of TMRCAs** per pair of individuals (multi-demes)
- Average **2D-SFS** (multi-demes)

###### Note 1 - expressed in units of branch lentgths (when applicable)
###### Note 2 - average statistics over all the genomic map, or subsetted by neutral and selection regions


### Appendix: type of topologies tested (N=2 diploids)

<img width="206" alt="image" src="https://github.com/user-attachments/assets/a2ba1ce3-f3c7-4315-9376-d88f857b5e5c" />

<img width="199" alt="image" src="https://github.com/user-attachments/assets/212cc0ba-fdec-4900-b058-d862c16d4af3" />

<img width="204" alt="image" src="https://github.com/user-attachments/assets/548dead5-9102-4dda-ab7d-d29be9ecd65b" />

<img width="204" alt="image" src="https://github.com/user-attachments/assets/f2d23fb9-19fc-41c9-9de1-7bb066e64f1a" />


###### Note - with all possible combinations of (haploid) individuals 
