# GTS - Pairwise Gene Tree Stats 
### Version: 0.1.1
---
#### Author: Pierre Lesturgie (pierrelesturgie@outlook.fr)

### (1) Installation: 

#### 1. Create env conda (optional): 

```conda create --name GTS python=3.12```

```conda activate GTS```

#### 2. Install: 

```pip install . ```

###### Uninstall : ```pip uninstall gts```

###### Note: if error occurs related to msprime: uninstall and reinstall using pip

### (2) Running: 

```gts <command> <arguments>``` (do ```gts --help``` for more information)
 
#### 1. STATS 
##### Computing summary and topology statistics for each tree in the .trees tree sequence.

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

#### 2. SUMMARY 
##### Summarize output stats from STATS. 

```gts summary --tag <treefile> -r <number_of_sampling_runs> -G <genomic_structure>```

##### Arguments: 

  -t, --treefile: Path to file with the tree sequence (omit '.trees')
  
  -r, --runs:  Number of runs: **must be the same** than in ```gts stats```
  
  -G, --genomic_map: Genomic map
  
##### Output: 
- Proportion of each **topology** per region
- Average **coalescent times** (C4, C3, C2) for each topology
- Distribution of **topologies** per quantile of **FST**
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

<img width="204" alt="image" src="https://github.com/user-attachments/assets/f9061e12-3041-4a18-a0e4-29eea7005291" />

<img width="204" alt="image" src="https://github.com/user-attachments/assets/3f05db78-6233-4151-a0b8-54eab50f8c65" />

#### Note - The total number of combinations under each topologies is computed with all possible combinations of N=4 (haploid) individuals irrespective of the deme of origin, EXCEPT for Topology 1 and 2 which have two subset topologies: 
With individuals 1 and 2 being from deme 1 and individuals 3 and 4 from deme 2: 
* Topology **Intra** groups indivuduals from _similar_ deme for the first event, i.e., ((1,2),(3,4))-like for Topology 1, ((3,(1,2)),4)-like for topology 2.
* Topology **Inter** groups individuals from _different_ demes for the first event, i.e., ((1,3),(2,4))-like for Topology 1, ((3,(1,4)),2)-like for topology 2.

