import dask.dataframe as dd
from typing import Dict, Any

def filter_overlapping_hits(hits_df: dd.DataFrame, schema_cols: Dict[str, Any]) -> dd.DataFrame:
    """Filters overlapping genomic hits within sample/contig groups, keeping the first non-overlapping entries.
    
    Groups by 'Sample' and 'Contig', sorts by start, and greedily selects non-overlapping entries.

    Args:
        hits_df: Input dask DataFrame with genomic hits
        schema_cols: Column schema with types

    Returns:
        dask DataFrame with filtered non-overlapping hits
    """
    def discard_overlaps(group_df: dd.DataFrame) -> dd.DataFrame:
        group_df = group_df.sort_values("Start").drop_duplicates(["Sample", "Contig", "Start"], keep="first")
        kept = []
        current_end = 0
        for idx, row in group_df.iterrows():
            if row["Start"] >= current_end:
                kept.append(idx)
                current_end = row["Stop"]
        return group_df.loc[kept]

    return hits_df.groupby(["Sample", "Contig"]).apply(discard_overlaps, meta=schema_cols)


def compute_zones_df(zones_df: dd.DataFrame, window: int) -> dd.DataFrame:
    """Computes ZoneIDs for genomic hits based on positional windows within sample/contig groups.
    
    Zones increment when the previous entry spans at least `window` length.

    Args:
        zones_df: Input dask DataFrame with genomic hits
        window: Minimum length to trigger a new zone

    Returns:
        dask DataFrame with an added 'ZoneID' column
    """
    zones_df = zones_df.reset_index(drop=True)
    zone = 1
    zone_list = [1]
    for i in range(1, len(zones_df)):
        if (zones_df.loc[i, "Sample"] != zones_df.loc[i-1, "Sample"]) or \
           (zones_df.loc[i, "Contig"] != zones_df.loc[i-1, "Contig"]):
            zone = 1
        else:
            if zones_df.loc[i-1, "Stop"] >= zones_df.loc[i-1, "Start"] + window:
                zone += 1
        zone_list.append(zone)
    zones_df["ZoneID"] = zone_list
    return zones_df


def context(hits_df: dd.DataFrame, window: int, schema_cols: Dict[str, Any]) -> dd.DataFrame:
    """Applies zone computation across dask DataFrame partitions while preserving order.

    Args:
        hits_df: Input dask DataFrame with genomic hits
        window: Minimum length to trigger a new zone
        schema_cols: Column schema (updated with 'ZoneID')

    Returns:
        dask DataFrame with 'ZoneID'
    """
    hits_df = hits_df.reset_index(drop=True)
    updated_cols = dict(schema_cols)
    updated_cols["ZoneID"] = "int64"
    return hits_df.map_partitions(lambda part: compute_zones_df(part, window), meta=updated_cols)

