#!/usr/bin/env python

import argparse
import pandas as pd

# do any processing of the samplesheet here
def samplesheet(samplesheet):
    ss = pd.read_csv(samplesheet, sep=',')
    ss.to_csv('samplesheet_utf8.csv', index=False, encoding='utf-8-sig')
    
    stats = ss.describe()
    stats.to_csv('samplesheet_stats.csv', index=True, encoding='utf-8-sig')
    
    print(ss.head())

def main():
    # initialize parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--samplesheet', 
                        metavar='samplesheet', 
                        type=str, 
                        help='sample metadata passed in through samples CSV file')

    args = parser.parse_args()

    samplesheet(args.samplesheet)

if __name__ == "__main__":
    main()