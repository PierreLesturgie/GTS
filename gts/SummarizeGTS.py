#!usr/bin/env python3.13

"""
***SummarizeGTS*** - Script to summarize computed statistics by GeneTreeStats
Author: Pierre Lesturgie
Version: 0.1.0
Last update: 2025-04-22
"""


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

    all_final,sel_final,neu_final = DataFrame(),DataFrame(),DataFrame()

    for r in tqdm(range(nruns), desc="Processing runs"):
        
        data = read_csv(f"{tag}_stats_{r}.detailed", sep=" ")
        data.insert(0, "loc_type", TF.classify_regions(data, chrom_pos, loc_types))

        span_sel = TF.compute_total_span(spans, loc_types, "S") if has_sel else 0
        span_neu = TF.compute_total_span(spans, loc_types, "N") if has_neu else 0
        span_all = span_sel + span_neu if has_sel or has_neu else sum(spans)

        sel_data = data[data['loc_type'] == "S"] if has_sel else DataFrame()
        neu_data = data[data['loc_type'] == "N"] if has_neu else DataFrame()

        all_final = concat([all_final,TF.get_stats(data,span_all)])

        if has_sel and not sel_data.empty:
            sel_final = concat([sel_final,TF.get_stats(sel_data,span_sel)])

        if has_neu and not neu_data.empty:
            neu_final = concat([neu_final,TF.get_stats(neu_data,span_neu)])

    final = concat([TF.summarize(all_final,"ALL"),TF.summarize(sel_final,"SELECTION"),TF.summarize(neu_final,"NEUTRAL")])
    
    final.to_csv(f"{tag}.sumstat", index=False)
    final.to_csv(f"{tag}.rawsumstat", header=False, index=False)

    print(f"Summary statistics saved to: {tag}.sumstat and {tag}.rawsumstat")