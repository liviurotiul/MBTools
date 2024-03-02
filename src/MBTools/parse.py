import pandas as pd
import argparse

def parse(infiles, outfile, cutoff):

    # Check if the infiles exist
    infiles = infiles.split(",")

    for file in infiles:
        try:
            open(file)
        except FileNotFoundError:
            print(f"The file {file} does not exist")
            exit()

    df_final = pd.DataFrame()

    for file in infiles:
            
        # 007_contig00046	T6081_AIA61620.1_KJ534568:c4122-3826_[Acinetobacter_baumannii]	100.0	98	0	0	351	644	1	98	1.2e-54	206.8

        df = pd.read_csv(file, sep='\t')

        # Name the columns
        df.columns = [
            "SAMPLE_CONTIG",
            "DIAMOND_ANNOTATION",
            "%IDENTITY",
            "ALIGNMENT_LENGTH",
            "MISMATCHES",
            "GAPS",
            "Q_START",
            "Q_END",
            "S_START",
            "S_END",
            "E_VALUE",
            "BIT_SCORE"
            ]

        df = df[df['%IDENTITY'] >= cutoff]

        df = df[['SAMPLE_CONTIG', 'DIAMOND_ANNOTATION', '%IDENTITY', 'Q_START', 'Q_END', 'MISMATCHES', 'ALIGNMENT_LENGTH', 'GAPS']]

        # Split SAMPLE_CONTIG into two columns
        df[['SAMPLE', 'CONTIG']] = df['SAMPLE_CONTIG'].str.split("_contig", expand=True)

        # If there is something after contig, remove it
        df['CONTIG'] = df['CONTIG'].str.split("_contig").str[0]

        # Drop the SAMPLE_CONTIG column
        df = df.drop(columns=["SAMPLE_CONTIG"])

        # Convert contig data to int from contig000001 to 1
        df['CONTIG'] = df['CONTIG'].astype(int)

        df["PREDICTOR"] = file.replace(".tsv", "")

        # Add the df to the final df
        df_final = pd.concat([df_final, df])

    # Set the data types
    dbtype_mapping = {
        'DIAMOND_ANNOTATION': 'str',
        '%IDENTITY': 'float64',
        'Q_START': 'int64',
        'Q_END': 'int64',
        'MISMATCHES': 'int64',
        'ALIGNMENT_LENGTH': 'int64',
        'GAPS': 'int64',
        'SAMPLE': 'str',
        'CONTIG': 'int64',
        'PREDICTOR': 'str'
    }

    df = df.astype(dbtype_mapping)

    df.to_csv(outfile, index=False)