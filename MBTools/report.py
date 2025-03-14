import dask.dataframe as dd
from typing import Dict, Any

def compute_zones_df(zones_df: dd.DataFrame) -> dd.DataFrame:
    """Computes ZoneIDs for genomic hits based on positional windows within sample/contig groups.
    
    Zones increment when the previous entry spans at least `window` length.

    Args:
        zones_df: Input dask DataFrame with genomic hits
        window: Minimum length to trigger a new zone

    Returns:
        dask DataFrame with added 'DrawStart' and 'DrawStop' columns
    """

    zones_df = zones_df.compute() if isinstance(zones_df, dd.DataFrame) else zones_df.copy()
    zones_df = zones_df.reset_index(drop=True)
    
    start = []
    stop = []
    
    if not zones_df.empty:
        # Handle the first row
        start.append(zones_df.loc[0, "Start"])
        stop.append(zones_df.loc[0, "Stop"])
        
        # Iterate from the second row to the end
        for i in range(1, len(zones_df)):

            # Previous entry's stop + 1
            new_start = stop[-1] + 10
            current_length = zones_df.loc[i, "Stop"] - zones_df.loc[i, "Start"]
            new_stop = new_start + current_length

            start.append(new_start)
            stop.append(new_stop)
        
    zones_df["DrawStart"] = start
    zones_df["DrawStop"] = stop
    return zones_df

def compute_draw_coordinates(hits_df: dd.DataFrame, schema_cols: Dict[str, Any]) -> dd.DataFrame:
    """Applies zone computation across dask DataFrame partitions while preserving order.

    Args:
        hits_df: Input dask DataFrame with genomic hits
        schema_cols: Column schema (updated with 'DrawStart' and 'DrawStop')

    Returns:
        dask DataFrame with 'DrawStart' and 'DrawStop'
    """
    updated_cols = dict(schema_cols)
    updated_cols.update({"ZoneID": "int64", "DrawStart": "int64", "DrawStop": "int64"})
    
    # Process each partition ensuring the schema is maintained
    return hits_df.map_partitions(
        compute_zones_df,
        meta=updated_cols
    )