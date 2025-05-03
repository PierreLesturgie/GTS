#!usr/bin/env python3.13

"""
***treefun*** - Module for Gene Tree Statistics Tool
Author: Pierre Lesturgie
Version: 0.1.0
Last update: 2025-04-24

IMPORTANT: Functions based on Barbara's scripts: function (1), (2), (12) and (18). 
"""

from math import factorial
from random import shuffle
from pandas import DataFrame, concat
import itertools
import numpy as np

# ************************************************************
# ******************** DEFINING FUNCTIONS ********************
# ************************************************************

### THERE ARE 17 FUNCTIONS

# <<<<<< Function 1: subsample trees for given individuals >>>>>>
def simplify_tree(ts_chroms, keep_indivs):
    sts_chroms = []
    for tree in ts_chroms:
        keep_nodes = []
        for i in keep_indivs:
            keep_nodes.extend(tree.individual(i).nodes)
            # https://tskit.dev/pyslim/docs/latest/metadata.html#sec-metadata
        sts = tree.simplify(keep_nodes, keep_input_roots=False)
        #print(f"The tree sequence has {sts.num_trees} trees\n")
        sts_chroms.append(sts)
    check_root(sts_chroms)
    return sts

# <<<<<< Function 2: Verify that there is only one root >>>>>>
def check_root(sts_chroms):
    for tree in sts_chroms:
        max_roots = max(t.num_roots for t in tree.trees())
        if max_roots != 1:
            print(f"nroots: {max_roots}")


# <<<<<< Function 3: converts a combination to a topology >>>>>> 
def translation(combination):
    topology = ""
    for i in range(len(combination)): 
        if len(combination[i]) == 1:
            if i == 0:
                topology = topology + f"({combination[i][0]},"
            else:
                topology = "(" + topology + f",{combination[i][0]})"
        elif len(combination[i]) == 2:
            if i == 0: 
                topology = topology + f"({combination[i][0]},{combination[i][1]})"
            else:
                topology = "(" + topology + f",({combination[i][0]},{combination[i][1]}))"
        elif len(combination[i]) == 3:
            if i == 0:
                topology = topology + f"({combination[i][0]},{combination[i][1]},{combination[i][2]})"
            else:
                topology = "(" + topology + f"({combination[i][0]},{combination[i][1]},{combination[i][2]}))"
        elif len(combination[i]) == 4:
            topology = f"({combination[i][0]},{combination[i][1]},{combination[i][2]},{combination[i][3]})"
    return topology

# <<<<<< Function 4: converts a list of external branch lengths to a topology >>>>>> 
def coal_times_to_combination(list_coal_times,TBL):
    list_indivs, combination = [1, 2, 3, 4], []
    
    sorted_coal_times_with_indivs = sorted(zip(list_coal_times, list_indivs))
    
    # This is an exception when C4=C3 but not starlike
    if len(set(list_coal_times)) == 1 and sum(list_coal_times) != TBL:
        combination = [[1,2],[3,4]]
    
    else: 
        
        while sorted_coal_times_with_indivs:
            coal_time, individual = sorted_coal_times_with_indivs.pop(0)
                
            indiv_coal_event = [individual]
                
            # Process duplicates by checking if the next element is the same as the current one.
            while sorted_coal_times_with_indivs and sorted_coal_times_with_indivs[0][0] == coal_time:
                    
                coal_time, individual = sorted_coal_times_with_indivs.pop(0)
                indiv_coal_event.append(individual)
                
            combination.append(indiv_coal_event)

    return translation(combination)

# <<<<<< Function 5: returns topology from external Branch lengths >>>>>> 
def topology(list_coal_times,topologies,TBL):
    
    combination = coal_times_to_combination(list_coal_times,TBL)
    #topologies = possible_topologies()
    
    for i in topologies:
        if combination in topologies[i].tolist(): 
            res = i
    return res

# <<<<<< Function 6: initiate first combination of a topology >>>>>> 
def initiate(topology_char,sequence_c,sequence):
    topology = topology_char
    for i in range(len(sequence)): topology = topology.replace(f'{sequence_c[i]}', f'{sequence[i]}')
    return [str(topology)]
    
# <<<<<< Function 7: computes all possible combinations for a topology >>>>>> 
def combinations(topology_char):
    
    sequence, sequence_to_shuffle, sequence_c = [1,2,3,4], [1,2,3,4], ["a","b","c","d"]

    ### Initiating first combination
    sequence_char = initiate(topology_char=topology_char,sequence_c=sequence_c,sequence=sequence)
    sequence_char = DataFrame({'Combination':sequence_char})
    
    for s in range(factorial(4)-1):
        
        shuffle(sequence_to_shuffle)
        temp_char = initiate(topology_char=topology_char,sequence_c=sequence_c,sequence=sequence_to_shuffle)

        while temp_char[0] in sequence_char["Combination"].tolist():
            shuffle(sequence_to_shuffle)     
            temp_char = initiate(topology_char=topology_char,sequence_c=sequence_c,sequence=sequence_to_shuffle)
        
        new_row=DataFrame({'Combination':temp_char})
        sequence_char = concat([sequence_char,new_row], ignore_index=True)
        
    return sequence_char 

# <<<<<< Function 8a: computes all possible combinations for all 4 topologies >>>>>> 
# <<<<<< DEPRECACATED >>>>>> 
def possible_topologies():


    ### TOPO 1
    #                     *
    #                    * *
    #                   *   *
    #                  *     *
    #                 *       *
    #                * *     * *
    #               *   *   *   *
    #              a     b c     d

    combinaison_char_1 = "((a,b),(c,d))"
    T1 = combinations(combinaison_char_1)


    ### TOPO 2
    #                     *
    #                    * *
    #                   *   *
    #                  *     *
    #                 * *     *
    #                *   *     *
    #               *   * *     *
    #              *   *   *     *
    #             a   b     c     d
    #     NEED ALSO TO MAKE THE REVERSE ONE FOR COMPUATION PURPOSES

    #combinaison_char_2 = "((a,(b,c)),d)"
    #combinaison_char_2_reverse = "(a,((b,c),d))"
    #T2 = combinations(combinaison_char_2)
    #T2_reverse = combinations(combinaison_char_2_reverse)
    #T2 = concat([T2,T2_reverse], ignore_index=True)


    ### TOPO 3
    #                     *
    #                    * *
    #                   *   *
    #                  *   * *
    #                 *   *   *
    #                *   *   * *
    #               *   *   *   *
    #              a   b   c     d
    #     NEED ALSO TO MAKE THE REVERSE ONE FOR COMPUATION PURPOSES

    combinaison_char_2 = "(a,(b,(c,d)))"
    combinaison_char_2_reverse = "(((a,b),c),d)"
    T2 = combinations(combinaison_char_2)
    T2_reverse = combinations(combinaison_char_2_reverse)
    T2 = concat([T2,T2_reverse], ignore_index=True)


    ### TOPO 4
    #                     *
    #                    * *
    #                   *   *
    #                  *     *
    #                 * *     *
    #                * * *     *
    #               *  *  *     *
    #              a    b  c     d
    ### NEED ALSO TO MAKE THE REVERSE ONE FOR COMPUATION PURPOSES

    combinaison_char_3 = "((a,b,c),d)"
    combinaison_char_3_reverse = "(a,(b,c,d))"
    T3 = combinations(combinaison_char_3)
    T3_reverse = combinations(combinaison_char_3_reverse)
    T3 = concat([T3,T3_reverse], ignore_index=True)


    ### TOPO 5 -- full star
        
    #                   *
    #                 ** **
    #                * * * *
    #               *  * *  *
    #              *   * *   *
    #             a    b c    d

    combinaison_char_4 = "(a,b,c,d)"
    T4 = combinations(combinaison_char_4)


    result = DataFrame({"Topology_1":T1["Combination"],
                        "Topology_2":T2["Combination"],
                        "Topology_3":T3["Combination"],
                        "Topology_4":T4["Combination"]})

    return result

# <<<<<< Function 8b: computes all possible combinations for all 4 topologies with subtopo >>>>>> 
def possible_topologies_with_subgroups():
    

    ### TOPO 1
    #                     *
    #                    * *
    #                   *   *
    #                  *     *
    #                 *       *
    #                * *     * *
    #               *   *   *   *
    #              a     b c     d

    combinaison_char_1 = "((a,b),(c,d))"
    T1 = combinations(combinaison_char_1)


    ### TOPO 2
    #                     *
    #                    * *
    #                   *   *
    #                  *     *
    #                 * *     *
    #                *   *     *
    #               *   * *     *
    #              *   *   *     *
    #             a   b     c     d
    #     NEED ALSO TO MAKE THE REVERSE ONE FOR COMPUATION PURPOSES

    #combinaison_char_2 = "((a,(b,c)),d)"
    #combinaison_char_2_reverse = "(a,((b,c),d))"
    #T2 = combinations(combinaison_char_2)
    #T2_reverse = combinations(combinaison_char_2_reverse)
    #T2 = concat([T2,T2_reverse], ignore_index=True)


    ### TOPO 3
    #                     *
    #                    * *
    #                   *   *
    #                  *   * *
    #                 *   *   *
    #                *   *   * *
    #               *   *   *   *
    #              a   b   c     d
    #     NEED ALSO TO MAKE THE REVERSE ONE FOR COMPUATION PURPOSES

    combinaison_char_2 = "(a,(b,(c,d)))"
    combinaison_char_2_reverse = "(((a,b),c),d)"
    T2 = combinations(combinaison_char_2)
    T2_reverse = combinations(combinaison_char_2_reverse)
    T2 = concat([T2,T2_reverse], ignore_index=True)


    ### TOPO 4
    #                     *
    #                    * *
    #                   *   *
    #                  *     *
    #                 * *     *
    #                * * *     *
    #               *  *  *     *
    #              a    b  c     d
    ### NEED ALSO TO MAKE THE REVERSE ONE FOR COMPUATION PURPOSES

    combinaison_char_3 = "((a,b,c),d)"
    combinaison_char_3_reverse = "(a,(b,c,d))"
    T3 = combinations(combinaison_char_3)
    T3_reverse = combinations(combinaison_char_3_reverse)
    T3 = concat([T3,T3_reverse], ignore_index=True)


    ### TOPO 5 -- full star
        
    #                   *
    #                 ** **
    #                * * * *
    #               *  * *  *
    #              *   * *   *
    #             a    b c    d

    combinaison_char_4 = "(a,b,c,d)"
    T4 = combinations(combinaison_char_4)

    # --- SUBTOPOLOGY GROUPING ---

    # Topology 1 subgroups
    group_1A = {"((1,2),(3,4))", "((2,1),(3,4))", 
                "((1,2),(4,3))", "((2,1),(4,3))",
                "((3,4),(1,2))", "((3,4),(2,1))", 
                "((4,3),(1,2))", "((4,3),(2,1))"}
    SubT1_A, SubT1_B = filter_topologies(T1, group_1A)

    # Topology 2 subgroups
    group_2A = {
        "(1,(2,(3,4)))", "(1,(2,(4,3)))",
        "(2,(1,(4,3)))", "(2,(1,(3,4)))",
        "(3,(4,(1,2)))", "(3,(4,(2,1)))",
        "(4,(3,(2,1)))", "(4,(3,(1,2)))",
        "(((1,2),3),4)", "(((2,1),3),4)",
        "(((1,2),4),3)", "(((2,1),4),3)",
        "(((3,4),1),2)", "(((3,4),2),1)",
        "(((4,3),1),2)", "(((4,3),2),1)",
    }
    SubT2_A, SubT2_B = filter_topologies(T2, group_2A)

    # Final dataframe with all subgroups
    result = DataFrame({
        "Topology_1_INTRA": SubT1_A["Combination"],
        "Topology_1_INTER": SubT1_B["Combination"],
        "Topology_2_INTRA": SubT2_A["Combination"],
        "Topology_2_INTER": SubT2_B["Combination"],
        "Topology_3": T3["Combination"],
        "Topology_4": T4["Combination"]
    })

    return result

# <<<<<< Function 9: Used to filter topologies 1 and 2 in two subtopologies >>>>>> 
def filter_topologies(df, group_set):
    ### THE ~ inverts the is in (i.e., returns false)
    return df[df["Combination"].isin(group_set)].reset_index(drop=True), df[~df["Combination"].isin(group_set)].reset_index(drop=True)

# <<<<<< Function 10: computes coalescent times given Internal, External branches and a topology >>>>>> 
def coalescent_times_n4(IB,EB,topo):
    if topo == 'Topology_1_INTRA' or topo == 'Topology_1_INTER':
        c4 = min(EB)
        EB.pop(np.argmin(EB)); EB.pop(np.argmin(EB))
        c3 = min(EB)
        c2 = c3 + min(IB)
    if topo == 'Topology_2_INTRA' or topo == 'Topology_2_INTER': 
        c4 = min(EB)
        EB.pop(np.argmin(EB)); EB.pop(np.argmin(EB))
        c3 = min(EB); EB.pop(np.argmin(EB))
        c2 = max(EB)
    if topo == 'Topology_3': 
        c4 = min(EB)
        EB.pop(np.argmin(EB)); EB.pop(np.argmin(EB)); EB.pop(np.argmin(EB))
        c3 = 0
        c2 = c4 + max(IB)
    if topo == 'Topology_4': 
        c4 = min(EB)
        c3 = c2 = 0
    return c4, c3, c2

# <<<<<< Function 11: computes initial objects used to compute statistics >>>>>> 
#### Now works for panmictic population
def initiate_stats(ts_chroms,pop_size,pop_sample,alive,rng):
    if len(pop_size) == 1:
        keep_indivs = rng.choice(alive, pop_sample[0], replace=False)
    else: 
        INDV=[]
        A=0
        for a in range(len(pop_size)):
            alive_temp=alive[(A):(pop_size[a]+A)]
            
            A += pop_size[a]
            keep_indivs = rng.choice(alive_temp, pop_sample[a], replace=False)
            INDV.append(keep_indivs)
        
        keep_indivs = np.append(INDV[0], INDV[1])
    
    sts_temp = simplify_tree(ts_chroms,keep_indivs)
    
    #print(popsamplediplo,keep_indivs)
    #print(len(popsamplediplo),len(keep_indivs))

    return sts_temp

# <<<<<< Function 12: Computes the average TMRCA over all the trees >>>>>>
def get_TMRCA(sts):
    av_tMRCA = 0
    for tree in sts.trees():
        av_tMRCA += tree.time(tree.root) * tree.span/sts.sequence_length
    return av_tMRCA

# <<<<<< Function 13: Computes pairwise TMRCAs between individuals >>>>>>
def get_pairwise_tmrcas(sts_temp):
    RES, COLN = [],[]
    samples = sts_temp.samples()
    haplotype_pairs = {(pair[0], pair[1]): [] for pair in itertools.combinations(samples, 2)}
    for tree in sts_temp.trees():
        for pair in itertools.combinations(samples, 2):
            haplotype_pairs[(pair[0], pair[1])].append(tree.tmrca(pair[0], pair[1]))
            average_tmrcas = np.zeros((sts_temp.num_samples, sts_temp.num_samples))
    average_tmrcas = np.zeros((sts_temp.num_samples, sts_temp.num_samples))
    for pair, tmrcas in haplotype_pairs.items():
        average_tmrcas[pair[0], pair[1]] = np.exp(np.mean(np.log(tmrcas)))
    
    # Code below is simply to save all TMRCAs in list format
    for i in range(average_tmrcas.shape[0]):
        for j in range(i,average_tmrcas.shape[1]):
            if i is not j:
                RES.append(average_tmrcas[i,j].tolist())
                COLN.append(f"TMRCA_individual{i}_{j}")  
    return RES, COLN

# <<<<<< Function 14: Computes 2D sfs, returns as list >>>>>>
def derived_2D_sfs(sts_temp,popsamplediplo,pop0=0,pop1=1):
    group1,group2,SFS_ALL,SFS_COLNAMES = [],[],[],[]
    #getting the haploid individuals to compute sfs
    for r in range(0,popsamplediplo[pop0]): group1.append(r)
    for r in range(popsamplediplo[pop0],popsamplediplo[pop1] + popsamplediplo[pop0]): group2.append(r)
    #print(f"sample for group 1 is {group1}")
    #print(f"sample for group 2 is {group2}")
    ### DEFINE GROUPS BEFORE!! 
    # 2D-SFS
    SFS = sts_temp.allele_frequency_spectrum([group1,group2], 
                                                mode="branch",span_normalise=True,polarised=True)
    #print(SFS)
    for d in range(popsamplediplo[pop0]+1):
        for f in range(popsamplediplo[pop1]+1):
            #print(d,f)
            SFS_ALL.append(SFS[d,f].tolist())
            SFS_COLNAMES.append(f"2DSFS_pop{pop0}_{d}_pop{pop1}_{f}")
            #print(f"calculating 2DSFS_pop{a}_{d}_pop{b+1}_{f}")
    #print(SFS_COLNAMES)
    return SFS_ALL, SFS_COLNAMES

# <<<<<< Function 15: Computes 1D sfs, returns as list >>>>>>
def derived_1D_sfs(sts_temp):
    SFS = sts_temp.allele_frequency_spectrum(mode="branch",polarised=True)
    SFS_ALL,SFS_COLNAMES = [],[]
    
    for s in range(len(SFS)):
        class_SFS = SFS.tolist()[s]
        SFS_ALL.append(class_SFS)
        SFS_COLNAMES.append(f"SFS_{s}")
        
    return SFS_ALL, SFS_COLNAMES

# <<<<<< Function 16: computes all pairwise summary statistics >>>>>>
def pairwise_summary_statistics(sts_temp, popsamplediplo, pop0=0, pop1=1):
    
    ### get tmrca from distribution of TMRCAS
    TOT_BR_LEN, COLNAMES, RESULT = [],[],[]
                
    TMRCA_btw = get_TMRCA(sts_temp)
    RESULT.append(TMRCA_btw)
    COLNAMES.append(f"TMRCA_{pop0}_{pop1}")
    
    A,B = sts_temp.samples(population=0) , sts_temp.samples(population=1)
    
    div=sts_temp.diversity(sample_sets=[A,B],mode="branch")
    RESULT.append(div[0]);RESULT.append(div[1])
    COLNAMES.append(f"pi_{pop0}");COLNAMES.append(f"pi_{pop1}")
    
    RESULT.append(sts_temp.num_mutations)
    COLNAMES.append(f"mutations_{pop0}_{pop1}")
                
    # computing total branch length --> used to compute dafi
    TOT_BR_LEN = sum(sts_temp.allele_frequency_spectrum(mode="branch"))
    RESULT.append(TOT_BR_LEN)
    COLNAMES.append(f"total_branch_length_{pop0}_{pop1}")
                        
    ### Computing DAFI
    RESULT.append(TMRCA_btw/TOT_BR_LEN)
    COLNAMES.append(f"dafi_{pop0}_{pop1}")
                        
    ### Computing DXY and FST
    dxy = sts_temp.divergence(sample_sets=[A,B],mode="branch")
    RESULT.append(dxy)
    COLNAMES.append(f"dxy_{pop0}_{pop1}")
    
    ### FST IS from TSKIT, not really Hudson's?
    #RESULT.append(sts_temp.Fst(sample_sets=[A,B],mode="branch"))
    #COLNAMES.append(f"fst_{pop0}_{pop1}")
    
    #HUDSON'S FST
    FST = (dxy - (div[0] + div[1])/2)/dxy
    RESULT.append(FST)
    COLNAMES.append(f"Hudson_fst_{pop0}_{pop1}")

    ### GET PAIRWISE TMRCAs
    pTMRCAS = get_pairwise_tmrcas(sts_temp)
    
    for u in range(len(pTMRCAS[0])):
        RESULT.append(pTMRCAS[0][u])
        COLNAMES.append(pTMRCAS[1][u])

    ### GET 2D-SFS
    SFS = derived_2D_sfs(sts_temp,popsamplediplo,pop0=pop0,pop1=pop1)
                    
    for u in range(len(SFS[0])):
        RESULT.append(SFS[0][u])
        COLNAMES.append(SFS[1][u])
        
    return RESULT, COLNAMES

# <<<<<< Function 17: computes summary statistics for a single deme/population >>>>>>
def one_pop_summary_statistics(sts_temp):
    COLNAMES, RESULT = [],[]
    
    RESULT.append(get_TMRCA(sts_temp))
    COLNAMES.append(f"TMRCA")
    
    RESULT.append(sts_temp.diversity(mode="branch"))
    COLNAMES.append(f"pi")
    
    RESULT.append(sts_temp.num_mutations)
    COLNAMES.append(f"mutations")
                
    # computing total branch length --> used to compute dafi
    RESULT.append(sum(sts_temp.allele_frequency_spectrum(mode="branch")))
    COLNAMES.append(f"total_branch_length")

    SFS = derived_1D_sfs(sts_temp)
    for u in range(len(SFS[0])):
        RESULT.append(SFS[0][u])
        COLNAMES.append(SFS[1][u])
    
    return RESULT, COLNAMES

# <<<<<< Function 18: Write branch length info, and returns coordinates of each tree >>>>>>
def write_branch_lengths_extract_coord(sts_temp,treefile):
    coordinates=[]
    iblFile = open(f'{treefile}_ibl', 'w', encoding="utf-8")
    eblFile = open(f'{treefile}_ebl', 'w', encoding="utf-8")
    for t in sts_temp.trees():
        iblFile.write("{} {} {} {} ".format(t.num_roots, t.interval[0], t.interval[1], t.total_branch_length))
        eblFile.write("{} {} {} {} ".format(t.num_roots, t.interval[0], t.interval[1], t.total_branch_length))
        coordinates.append([t.interval[0],t.interval[1]])
        a = 0
        for u in t.nodes():
            if t.is_leaf(u):
               # print(f"{t.branch_length(u)}")
                eblFile.write("{} ". format(t.branch_length(u)))
            if not t.is_leaf(u): 
                iblFile.write("{} ".format(t.branch_length(u)))
                #print(t.branch_length(u))
                a += 1
        if a != 3:
            iblFile.write("{} ".format(0))
        iblFile.write("\n")
        eblFile.write("\n")
    iblFile.close()
    eblFile.close()
    return coordinates 



### THESE FUNCTIONS ARE FOR THE SUMMARIZING

# <<<<<< Function 19: Summarize summary statistics output by GTS >>>>>>
def summary_sumstats(data, span):
    data['diff'] = data['END'] - data['START']
    span_weights = data['diff'] / span
    data = data.drop(columns=['loc_type', 'START', 'END', 'Topology', 'C4', 'C3', 'C2', 'Unnamed: 6', 'Unnamed: 29', 'diff'], errors='ignore')
    
    result = (data.mul(span_weights, axis=0)).sum().tolist()
    return result, data.columns.tolist()

# <<<<<< Function 20: Summarize topology statistics output by GTS >>>>>>
def summary_topology_coalescence_times(data):
    result, names = [], []
    topologies = ["Topology_1_INTRA", "Topology_1_INTER","Topology_2_INTRA", "Topology_2_INTER", "Topology_3", "Topology_4"]
    
    base = data[["START", "END", "Topology", "C4", "C3", "C2"]].copy()
    base["diff"] = base["END"] - base["START"]
    base["SPAN"] = base["diff"] / base["diff"].sum()
    
    for topo in topologies:
        filtered = base[base["Topology"] == topo]
        prop = filtered["SPAN"].sum()
        result.extend([prop])
        names.append(f"{topo}")
        
        weighted = filtered.copy()
        for col in ["C4", "C3", "C2"]:
            weighted[col] *= weighted["SPAN"] / prop ## ADDED THIS ONE TO CORRECT BY TOPOLOGY SPAN
            result.append(weighted[col].sum())
            names.append(f"{topo}_{col}")
    
    return result, names

def reshape_results(topo, sumstat, qtfst=None):
    values = topo[0] + sumstat[0] 
    columns = topo[1] + sumstat[1] 
    if qtfst is not None:
        values = values + qtfst[0] 
        columns = columns + qtfst[1] 
    
    return DataFrame([values], columns=columns)

# <<<<<< Function 21: compute spans values >>>>>>
def compute_spans(genome):
    chrom_pos = [0] + genome.iloc[0].tolist()
    loc_types = ['0'] + genome.iloc[1].tolist()
    chrom_pos = [int(x) + 1 for x in chrom_pos]
    spans = [chrom_pos[i + 1] - chrom_pos[i] for i in range(len(chrom_pos) - 1)]
    spans[0] += 1
    return chrom_pos, loc_types, spans

# <<<<<< Function 22: compute total span for each region >>>>>>
def compute_total_span(spans, loc_types, target):
    return sum(spans[i - 1] for i, x in enumerate(loc_types) if x == target)

# <<<<<< Function 23: classify regions (i.e., neutral, selection, all) >>>>>>
def classify_regions(data, chrom_pos, loc_types):
    loc = []
    for j in range(data.shape[0]):
        pos = data.iloc[j, 1]
        for i in range(len(chrom_pos) - 1):
            if chrom_pos[i] <= pos < chrom_pos[i + 1]:
                loc.append(loc_types[i + 1])
                break
    return loc

# <<<<<< Function 24: summarize dataset >>>>>>
### HERE FIND A SOLUTION FOR QFST: not removing 0!! 
def summarize(df, label):
    topology_colnames = ["Topology_1_INTRA", "Topology_1_INTER","Topology_2_INTRA",
                            "Topology_2_INTER","Topology_3","Topology_4"]
    
    # filtering for FST. Normally, if no FST, returns nothing. 
    # This is to not replace the 0s in FST by NA
    df_fst = df.filter(regex='^FST')
    df = df[df.columns.drop(list(df.filter(regex='^FST')))]
    
    df_temp = df
    df = df.replace(0, np.nan)
    
    for i in topology_colnames:
        df[f"{i}"] =  df_temp[f"{i}"] 
    
    df_merged = concat([df,df_fst], ignore_index=True)
    means = df_merged.mean().to_frame().T
    means["LOC_TYPE"] = label
    return means

# <<<<<< Function 25: returns the final dataset >>>>>>
def get_stats(data,span):
    df = DataFrame()
    s_sum = summary_sumstats(data, span)
    s_topo = summary_topology_coalescence_times(data)
    if 'Hudson_fst_0_1' in data.columns:
        s_qt = qt_fst(data)
        df = concat([df, reshape_results(s_topo, s_sum, qtfst=s_qt)], ignore_index=True)
    else: 
        df = concat([df, reshape_results(s_topo, s_sum, qtfst=None)], ignore_index=True)
    return df


# Function to compute FST proportions and column names per topology
def topo_span_fst(df, tag):
    topologies = [
    "Topology_1_INTRA", "Topology_1_INTER", "Topology_2_INTRA",
    "Topology_2_INTER", "Topology_3", "Topology_4"
    ]
    total = df["diff"].sum()
    colnames = [f"FST_{tag}_{topo}" for topo in topologies]
    values = [df[df["Topology"] == topo]["diff"].sum() / total for topo in topologies]
    return values, colnames

# Main function to extract FST summaries
def qt_fst(data):
    # Calculate quantile thresholds
    thresholds = data["Hudson_fst_0_1"].quantile([0.05, 0.25, 0.5, 0.75, 0.95]).tolist()
    quantile_tags = ["0.05", "0.25", "0.5", "0.75", "0.95"]
    
    # Define bins for quantiles
    bins = [
        data[data["Hudson_fst_0_1"] <= thresholds[0]],
        data[(data["Hudson_fst_0_1"] > thresholds[0]) & (data["Hudson_fst_0_1"] <= thresholds[1])],
        data[(data["Hudson_fst_0_1"] > thresholds[1]) & (data["Hudson_fst_0_1"] <= thresholds[2])],
        data[(data["Hudson_fst_0_1"] > thresholds[2]) & (data["Hudson_fst_0_1"] <= thresholds[3])],
        data[data["Hudson_fst_0_1"] > thresholds[4]]
    ]
    
    # Compute results and column names
    result, colnames = [], []
    for tag, bin_df in zip(quantile_tags, bins):
        values, names = topo_span_fst(bin_df, tag)
        result.extend(values)
        colnames.extend(names)
    
    return result, colnames