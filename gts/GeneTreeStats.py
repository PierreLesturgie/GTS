#!usr/bin/env python3.13

"""
***GeneTreeStats*** - Script for computing statistics
Author: Pierre Lesturgie
Version: 0.1.1
Last update: 2025-05-07
"""


from pandas import read_csv
from numpy import unique, random, sum
from tqdm import tqdm
from subprocess import run
import pyslim, tskit
import gts.treefun as TF  # updated import
import argparse


def run_gene_tree_stats(treefile, Nanc, runs,rho,demographic_scenario,genomic_map):
    #print(f"[GeneTreeStats] Running with treefile: {treefile}, Nanc: {Nanc}")

    reco_rate = rho

    # Load population scenario
    populations = read_csv(demographic_scenario, sep=" ")
    populations = populations.loc[:, ~populations.columns.str.contains('^Unnamed')]

    # Genome map (deprecated)
    genome = read_csv(genomic_map, sep=" ")
    chrom_positions = [0, int(genome.iloc[0].tolist()[-1])]

    print("Chromosome positions:", chrom_positions)

    ts = tskit.load(f"{treefile}.trees")
    print(f"The tree sequence has {ts.num_trees} trees on a genome of length {ts.sequence_length}, "
          f"{ts.num_individuals} individuals, {ts.num_samples} samples, and {ts.num_mutations} mutations.")

    # Recapitate
    rts = pyslim.recapitate(ts, recombination_rate=reco_rate, ancestral_Ne=Nanc)

    # Extract chromosome intervals
    ts_chroms = []
    for j in range(len(chrom_positions) - 1):
        start, end = chrom_positions[j:j + 2]
        chrom_ts = rts.keep_intervals([[start, end]], simplify=False).trim()
        ts_chroms.append(chrom_ts)
        print(f"Chrom {j+1}: {chrom_ts.num_samples} samples, {chrom_ts.num_individuals} individuals, "
              f"{chrom_ts.num_sites} SNPs, length {chrom_ts.sequence_length}")

    chrom_ts = ts_chroms[-1]
    alive = pyslim.individuals_alive_at(chrom_ts, 0)
    print(f"{len(alive)} individuals alive in the final generation.")

    individual_times = chrom_ts.individuals_time
    for t in unique(individual_times):
        print(f"{sum(individual_times == t)} individuals from time {t}.")

    rng = random.default_rng()

    # Population setup
    pop_sample = populations.iloc[0].tolist()
    pop_size = populations.iloc[1].tolist()
    number_of_demes = len(pop_size)

    if number_of_demes > 2:
        print("Warning: code not tested for more than 2 demes.")

    popsamplediplo = [x * 2 for x in pop_sample]
    topologies = TF.possible_topologies_with_subgroups()

    # Sampling loop
    for i in tqdm(range(runs), desc="Sampling"):
        sts_temp = TF.initiate_stats(ts_chroms, pop_size, pop_sample, alive, rng)
        coordinates = TF.write_branch_lengths_extract_coord(sts_temp, treefile)

        run(f"paste {treefile}_ibl {treefile}_ebl > {treefile}_temp", shell=True)
        run(f"cut -d' ' -f2,3,4,6,7,12,13,14,15,16,17,18,19 {treefile}_temp > {treefile}_BL", shell=True)

        # Topology analysis
        with open(f"{treefile}_BL") as bl, open(f"{treefile}_topology.stats", 'w') as topo_out:
            topo_out.write("START END Topology C4 C3 C2\n")
            for line in bl:
                L = [float(x) for x in line.split()]
                TBL, IBL = L[2], L[3:5]
                a, b =  [5,7,9,11], [6,8,10,12]
                r1, r2 = [L[i] for i in a], [L[i] for i in b]
                s=sorted(zip(r1,r2))
                EBL = [s[0][1],s[1][1],s[2][1],s[3][1]]
                topo = TF.topology(EBL, topologies, TBL)
                CT = TF.coalescent_times_n4(IBL, EBL, topo)
                topo_out.write(f"{L[0]} {L[1]} {topo} {CT[0]} {CT[1]} {CT[2]}\n")

        run(f"rm {treefile}_ibl {treefile}_ebl {treefile}_temp {treefile}_BL", shell=True)

        # Summary statistics
        with open(f"{treefile}_summary.stats", 'w') as summary:
            for w in range(len(coordinates)):
                interval = [coordinates[w]]
                ts_window = sts_temp.keep_intervals(interval, simplify=False).trim()
                if len(pop_sample) == 1:
                    stat = TF.one_pop_summary_statistics(ts_window)
                else:
                    stat = TF.pairwise_summary_statistics(ts_window, popsamplediplo)
                summary.write(" ".join(map(str, stat[0])) + "\n")

        with open(f"{treefile}_header.stats", 'w') as header:
            header.write(" ".join(stat[1]) + "\n")

        run(f"cat {treefile}_header.stats {treefile}_summary.stats > {treefile}_sumstat.temp", shell=True)
        run(f"paste -d ' ' {treefile}_topology.stats {treefile}_sumstat.temp > {treefile}_stats_{i}.detailed", shell=True)

        run(f"rm {treefile}_topology.stats {treefile}_summary.stats {treefile}_sumstat.temp {treefile}_header.stats", shell=True)
    pass

