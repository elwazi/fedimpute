#!/usr/bin/env python

import argparse
import pandas as pd
import os

def generate_html(samplefiles, output):
    # Start the HTML file
    with open(output, "w") as f:
        f.write("""
        <!DOCTYPE html>
        <html>
        <head>
        <title>Imputation report</title>
        <style>
        table {
          border-collapse: collapse;
          width: 50%;
        }

        th, td {
          border: 1px solid black;
          padding: 8px;
          text-align: left;
        }
        
        thead > tr > th.rotate:not(:first-child) {
          /* Tilt the headers */
          writing-mode: vertical-lr;
          /* Set text orientation */
          transform: rotate(180deg);
        }
        </style>
        </head>
        <body>
        <h1>Imputation report</h1>
        <h3>Sample Membership</h2>
        """)

    table_html = generate_html_crosstab(samplefiles)
    # add "rotate" class for the crosstab headers
    with open(output, "a") as f:
        f.write(table_html.replace('<th>', '<th class="rotate">'))

    # End HTML file
    with open(output, "a") as f:
        f.write("</body></html>")

    print("HTML file generated:", output)

def generate_html_crosstab(samplefiles):
    dfs = []
    for samplefile in samplefiles:
        df = pd.read_csv(samplefile, names=["SampleID"])
        filename = samplefile[:samplefile.index(".vcf")]
        df["Filename"] = os.path.basename(filename)
        dfs.append(df[["Filename","SampleID"]])
    dfx = pd.concat(dfs)
    dfh = pd.crosstab(dfx['SampleID'], dfx['Filename'])
    # dfh = pd.crosstab(dfx['Filename'], dfx['SampleID'])
    return dfh.to_html()    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate HTML file from a list of text files.")
    parser.add_argument("--samples", nargs="+", type=str, help="Space-separated list of sample file names.")
    parser.add_argument("--output", type=str, help="Name of the output HTML file.")
    args = parser.parse_args()

    if not args.samples or not args.output:
        print("Please provide both --samples and --output arguments.")
    else:
        generate_html(args.samples, args.output)
