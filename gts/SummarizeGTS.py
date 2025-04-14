#!usr/bin/env python3.13
# SummarizeGTS.py

from pandas import DataFrame, concat, read_csv
import warnings
from tqdm import tqdm
import argparse
import gts.treefun as TF  # updated import

warnings.filterwarnings("ignore")

def run_summary(tag, runs, genomic_map):
    print(f"[SummarizeGTS] Tag: {tag}, Runs: {runs}, Genome map: {genomic_map}")
    
    nruns = int(runs)
    genome = read_csv(genomic_map, sep=" ")

    chrom_pos, loc_types, spans = TF.compute_spans(genome)

    # Check which types are present
    has_sel = "S" in loc_types
    has_neu = "N" in loc_types

    for r in tqdm(range(nruns), desc="Processing runs"):
        
        data = read_csv(f"{tag}_stats_{r}.detailed", sep=" ")
        data.insert(0, "loc_type", TF.classify_regions(data, chrom_pos, loc_types))

        span_sel = TF.compute_total_span(spans, loc_types, "S") if has_sel else 0
        span_neu = TF.compute_total_span(spans, loc_types, "N") if has_neu else 0
        span_all = span_sel + span_neu if has_sel or has_neu else sum(spans)

        sel_data = data[data['loc_type'] == "S"] if has_sel else DataFrame()
        neu_data = data[data['loc_type'] == "N"] if has_neu else DataFrame()

        final = TF.summarize(TF.get_stats(data,span_all), "ALL")

        if has_sel and not sel_data.empty:
            final = concat([final, TF.summarize(TF.get_stats(sel_data,span_sel), "SELECTION")], ignore_index=True)

        if has_neu and not neu_data.empty:
            final = concat([final, TF.summarize(TF.get_stats(neu_data,span_neu), "NEUTRAL")], ignore_index=True)


    final.to_csv(f"{tag}.sumstat", index=False)
    final.to_csv(f"{tag}.rawsumstat", header=False, index=False)

    print(f"Summary statistics saved to: {tag}.sumstat and {tag}.rawsumstat")



def main():
    parser = argparse.ArgumentParser(description='SummarizeGTS: summarize GeneTreeStats output')
    parser.add_argument("-t", "--tag", required=True, help="Tag of the replicate")
    parser.add_argument("-r", "--runs", required=True, help="Number of runs", default=100)
    parser.add_argument("-G", "--genomic_map", required=True, help="Genomic map", default="genome.txt")
    args = parser.parse_args()

    run_summary(args.tag, args.runs, args.genomic_map)

if __name__ == "__main__":
    main()