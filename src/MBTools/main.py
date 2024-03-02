import pandas as pd
import argparse
from MBTools.parse import parse
from MBTools.context import context

def main() -> None:

    # Add command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--parse', default=None, help="The parse functionality is used")
    parser.add_argument('-i','--infiles', help="Comma separated list of input files")
    parser.add_argument('-o','--outfile')
    parser.add_argument('-f','--filter', default=90.00, type=float, help="The cutoff for %IDENTITY (default: 90.00)")
                        
    parser.add_argument('-c','--context', default=None, help="The context functionality is used")
    parser.add_argument('-v', '--visualiser', type=bool, default=False, help='if to create visualiser')
    parser.add_argument('-w', '--window', type=int, default=1000, help='window size')

    args = parser.parse_args()
    infiles = args.infiles
    outfile = args.outfile
    cutoff = args.filter
    parser = args.parse
    context_analysis = args.context
    window = args.window

    if parser:

        parse(infiles, outfile, cutoff)
    
    if context_analysis:

        context(infiles, window)


if __name__ == "__main__":
    main()