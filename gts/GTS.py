#!usr/bin/env python3.13

"""
GTS - Gene Tree Statistics Tool
Author: Pierre Lesturgie
Version: 0.1.0
Last update: 2025-04-14
"""

import argparse
from .GeneTreeStats import run_gene_tree_stats
from .SummarizeGTS import run_summary
from ._version import __version__, __last_update__

def main():
    print(f"GTS - Gene Tree Statistics Tool v{__version__} (Last update {__last_update__})")
    
    parser = argparse.ArgumentParser(description="GeneTree Toolkit")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # GeneTreeStats args
    stats_parser = subparsers.add_parser("stats", description='GeneTreeStats (<< stats >>): computation of pairwise gene tree-based statistics',help="Run gene tree stats")
    
    stats_parser.add_argument("-t", "--treefile", required=True,
                        help="Path to file with the tree sequence (omit '.trees')")
    stats_parser.add_argument("-Nanc", "--ancestral_effective_size", required=True, type=int,
                        help="Ancestral effective size")
    stats_parser.add_argument("-r", "--runs", required=True,
                        help="Number of sampling runs", default=100, type=int)
    stats_parser.add_argument("-R", "--rho", required=True,
                        help="Recombination rate", default=2.5e-8, type=float)
    stats_parser.add_argument("-D", "--demographic_scenario", required=True,
                        help="Demographic scenario file", default="scenario.txt")
    stats_parser.add_argument("-G", "--genomic_map", required=True,
                        help="Genomic map (DEPRECATED)", default="genome.txt")

    # SummarizeGTS args
    summary_parser = subparsers.add_parser("summary", description='SummarizeGTS (<< summary >>): summarize GeneTreeStats output',help="Summarize GTS results")
    summary_parser.add_argument("-t", "--tag", required=True, help="Tag of the replicate")
    summary_parser.add_argument("-r", "--runs", required=True, help="Number of runs", default=100)
    summary_parser.add_argument("-G", "--genomic_map", required=True, help="Genomic map", default="genome.txt")

    args = parser.parse_args()

    if args.command == "stats":
        run_gene_tree_stats(args.treefile, args.ancestral_effective_size, args.runs, args.rho, args.demographic_scenario, args.genomic_map)
    elif args.command == "summary":
        run_summary(args.tag, args.runs, args.genomic_map)

if __name__ == "__main__":
    main()