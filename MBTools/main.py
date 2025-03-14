import argparse
import os
import dask.dataframe as dd
import shutil

from jinja2 import Template
from context import filter_overlapping_hits, context
from report import compute_draw_coordinates

columns = {
    "Sample": "str",
    "Contig": "str",
    "Start": "int64",
    "Stop": "int64",
    "Strand": "str",
    "Gene": "str",
    "Gaps": "int64",
    "Identity": "float64",
    "Coverage": "float64",
    "Category": "str",
    "Product": "str"
}

def main():

    ###################################################################
    # Setup
    ###################################################################

    # Add command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', help="Comma separated list of input files")
    parser.add_argument('-o','--output', help="Output folder")
    parser.add_argument('-f','--filter', type=float, help="The cutoff for %Identity (default: 90.00)")
    parser.add_argument('-r', '--report', default=False, help='if to create report', action="store_true")
    parser.add_argument('-v', '--visualiser', type=bool, default=False, help='if to create visualiser')
    parser.add_argument('-w', '--window', type=int, default=1000, help='window size')
    parser.add_argument('--overwrite', default=False, help='if to overwrite output folder', action="store_true")


    # Parse arguments
    args = parser.parse_args()
    input = args.input
    output = args.output
    window = args.window
    overwrite = args.overwrite
    id_filter = args.filter

    # Validate input and output paths
    if not input:
        raise ValueError("Input files are required")
    if not output:
        raise ValueError("Output file is required")
    
    if not os.path.exists(input):
        raise FileNotFoundError(f"Input file not found: {input}")
    
    if os.path.exists(output) and not overwrite:
        raise FileExistsError(f"Output folder already exists: {output}")
    
    if os.path.exists(output) and overwrite:
        shutil.rmtree(output)
    
    ###################################################################
    # Read input file csv
    ###################################################################
    
    os.mkdir(output)

    # Get the columns from the input file by reading the first row
    input_cols = dd.read_csv(input, dtype=str).head(1).columns

    for col in ["Sample", "Contig", "Start", "Stop", "Strand", "Gene"]:
        assert col in input_cols, f"Column {col} not found in input file"
    
    cols_to_read = set(columns.keys()).intersection(input_cols)

    # Read the input csv and ignore columns that are not in the columns dict
    df = dd.read_csv(input, dtype=columns, usecols=cols_to_read)
    df = df.dropna() # TO BE REMOVED 
    assert not df.isnull().values.any(), "Input file contains NaN values"

    df_cols = {col: df[col].dtype for col in df.columns}

    if id_filter:
        # Chech for Identity column in df
        assert "Identity" in df.columns, "Filter option selected but Identity column was not found in input file"

    ###################################################################
    # Analyze context
    ###################################################################

    # Filter identity hits
    df = df[df["Identity"] >= id_filter] if id_filter else df

    # Compute context
    df = filter_overlapping_hits(df, df_cols)
    df = context(df, window, df_cols)

    if not args.report:
        df.to_csv(f"{output}/results.csv", index=False, single_file=True)
        
    else:

        df = compute_draw_coordinates(df, df_cols)
        df.to_csv(f"{output}/results.csv", index=False, single_file=True)

        # Compute zones of length at least 3
        # Give each zone a unique ID
        df["ZoneUID"] = df["Sample"].astype(str) + "_" + df["Contig"].astype(str) + "_" + df["ZoneID"].astype(str)

        # Find which zones have more than 3 genes
        zones = df.groupby("ZoneUID").size().compute()
        zones = zones[zones > 3].index

        # Keep only the zones with more than 3 genes
        df = df[df["ZoneUID"].isin(zones)]

        # Save the results
        df.to_csv(f"{output}/results_top.csv", index=False, single_file=True)


        with open(f"{output}/results_top.csv") as f:
            csv_str = f.read()
        
        with open("index.html") as f:
            template_str = f.read()

        template = Template(template_str)
        report = template.render(df=csv_str)

        with open(f"{output}/report.html", "w") as f:
            f.write(report)


if __name__ == "__main__":
    main()