"""Process ATAC-seq files for visualization in igv.js

PREREQUISITE:
- Install htslib: `brew install htslib` on macOS (TODO: Create Docker file)

"""
import os
import shutil
import subprocess

def index_fragments(tsv_path):
    """Bgzip and TBI-index a raw ATAC-seq fragments file
    """
    print("Index fragments file from DSP Pipelines team")

    print("Change suffix from .tsv to .bed, as expected by igv.js")
    tmp_tsv_path = f"{tsv_path}.tmp"
    shutil.copy2(tsv_path, tmp_tsv_path) # Preserve original file
    bed_path = tmp_tsv_path.replace(".tsv.tmp", ".bed")
    os.rename(tmp_tsv_path, bed_path)

    print("Sort BED file by genomic position")
    # The "-k1.4V" argument ensures `sort` uses column "1" and breaks ties
    # by ordering character "4" using version-sort, i.e. "V", which does a
    # natural sort of (version) numbers within text.  This ensures e.g. "12"
    # appears after "2".
    sorted_bed_path = bed_path.replace(".bed", ".possorted.bed")
    sort_command = (f"sort -k1.4V {bed_path}").split(" ")
    sorted_bed_file = open(sorted_bed_path, 'w')
    subprocess.call(sort_command, stdout=sorted_bed_file)

    print("Compress sorted BED file with blocked gzip")
     # bgzip enables requesting small indexed chunks of a gzip'd file
    bgzip_command = (f"bgzip -kf {sorted_bed_path}").split(" ")
    subprocess.call(bgzip_command)

    print("Make a TBI index of the bgzipped BED file")
    # tabix creates an index for the tab-delimited files, used for getting small chunks
    tabix_command = (f"tabix -s 1 -b 2 -e 3 {sorted_bed_path}.gz").split(" ")
    subprocess.call(tabix_command)

index_fragments("Library-8-20230710_atac.fragments.tsv")

