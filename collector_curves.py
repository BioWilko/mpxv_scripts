import pysam
import random
import argparse
from pathlib import Path
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--infile")
parser.add_argument("-t", "--threshold", type=int, default=20)
args = parser.parse_args()

collector_curve = []

with pysam.AlignmentFile(args.infile, "rb") as bam_fh:
    coverage = np.zeros(bam_fh.lengths[0], dtype=int)

    read_counter = 0

    segments = [x for x in bam_fh]
    random.shuffle(segments)

    for segment in segments:
        if segment.is_unmapped or segment.is_supplementary:
            continue

        read_counter += 1

        coverage[segment.reference_start : segment.reference_end] += 1

        collector_curve.append(
            (str(read_counter), str(np.sum(coverage >= args.threshold)))
        )

print("reads\tcovered_bases")
for i in collector_curve:
    print(f"{i[0]}\t{i[1]}")
