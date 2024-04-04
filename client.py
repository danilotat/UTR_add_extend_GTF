#!/usr/bin/python3

from gtftools.workers import *
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Process some integers.")
    parser.add_argument("--input", "-i", help="Input GTF file", required=True)
    parser.add_argument("--output", "-o", help="Output GTF file", required=True)
    parser.add_argument("--length", "-l", help="Extension length", required=True)
    parser.add_argument(
        "--min_dist", "-m", help="Minimum distance from the next gene", required=True
    )
    parser.add_argument("--logs", help="Output for logs", required=False)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    chromosome_wise(args.input, args.output, args.length, int(args.min_dist), args.logs)
