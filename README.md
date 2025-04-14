# GTS - Pairwise Gene Tree Statistics Tools
---
### Author: Pierre Lesturgie (pierrelesturgie@outlook.fr)

### (1) Installation: 

#### 1. Create env conda: 

```conda create --name GTS python=3.12```

```conda activate GTS```

#### 2. Install: 

```pip install . ```

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

##### Output ######(in branch lentgths): 
- Topology (out four kind of topologies; see below)
- Coalescent times: C4, C3, C2 
- TMRCA
- Pi(s)
- 1D-SFS(s)
- Total Branch Length
- DAFi (only multi-demes)
- dxy (only multi-demes) 


#### 2. SUMMARY. Summarize output stats from STATS. 

```gts summary --tag <treefile> -r <number_of_sampling_runs> -G <genomic_structure>```

##### Arguments: 

  -t, --treefile: Path to file with the tree sequence (omit '.trees')
  
  -r, --runs:  Number of runs ** Must be the same than in ```gts stats```
  
  -G, --genomic_map: Genomic map
  
