import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm

from MBTools.context_utils import create_zone_column, create_tag_column, merge_overlapping_genes, concatenate_overlapping_genes, concatenate_overlapping_gene_names
from MBTools.context_utils import merge_overlapping_genes, concatenate_overlapping_genes, concatenate_overlapping_gene_names, translator

def context(infile, window):

    samples = None

    # Read the source file
    df = pd.read_csv(infile, low_memory=False, header=0)

    # Create a tag column
    df["TAG"] = create_tag_column(df)

    # Add zones
    df["ZONE_ID"] = create_zone_column(df, window)
    
    for index, row in df.iterrows():
        df.at[index, 'SHORT_ANNOTATION'] = translator(row['DIAMOND_ANNOTATION'], row['PREDICTOR']).replace(",", "")

    total_length = []
    consensus_tags = []
    tags = []
    genes = []
    gene_count = []
    start = []
    end = []
    annotation_count = []
    annotation_id = []
    samples = []
    zone_id_list = df["ZONE_ID"].unique().tolist()
    contigs = []

    for zone in tqdm(zone_id_list):

        zone_df = df[df["ZONE_ID"]==zone]
        zone_df = zone_df.sort_values(by=["Q_START"]).reset_index(drop=True)
        total_length.append(zone_df.iloc[-1]["Q_END"] - zone_df.loc[0]["Q_START"])
        consensus_gene_tags = '@'.join([x for x in merge_overlapping_genes(zone_df)])
        consensus_tags.append(consensus_gene_tags)
        gene_tags = '@'.join([x for x in concatenate_overlapping_genes(zone_df)])
        tags.append(gene_tags)

        gene_names = '@'.join([x for x in concatenate_overlapping_gene_names(zone_df)])
        genes.append(gene_names)
        annotation_count.append(len(zone_df))
        gene_count.append(len(consensus_gene_tags.split('@')))
        start.append(zone_df.loc[0]["Q_START"])
        end.append(zone_df.iloc[-1]["Q_END"])

        contigs.append(zone_df["CONTIG"].unique().tolist()[0])
        samples.append(zone_df["SAMPLE"].unique().tolist()[0])

    zone_df = {
        "ZoneID": zone_id_list,
        "TotalLength": total_length,
        "GeneCount": gene_count,
        "AnnotationCount": annotation_count,
        "Q_START": start,
        "Q_END": end,
        "ConsensusTags": consensus_tags,
        "Tags": tags,
        "Genes": genes,
        "Sample": samples,
        "Contig": contigs
    }

    zone_df = pd.DataFrame(zone_df)
    infile = infile.split(".")[0]
    zone_df.to_csv(f"{infile}_zones.csv", index=False)
    df.to_csv(f"{infile}_tags.csv", index=False)
