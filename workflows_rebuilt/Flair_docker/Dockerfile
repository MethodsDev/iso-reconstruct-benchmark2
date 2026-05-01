# Use miniconda3 as base image
FROM continuumio/miniconda3

RUN apt-get update && apt-get install -y --no-install-recommends build-essential zlib1g-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install mamba
RUN conda install -c conda-forge mamba -y && conda clean -a -y

# Create a new environment and install required system tools
RUN mamba create -n flair_env -c bioconda -c conda-forge python=3.12 samtools bedtools minimap2 pip && \
    conda clean -a -y

# Activate the environment and set it as the default
SHELL ["conda", "run", "-n", "flair_env", "/bin/bash", "-c"]

# Set the working directory
WORKDIR /usr/src/app

# Set the PATH to include the conda environment
ENV PATH=/opt/conda/envs/flair_env/bin:$PATH

# Install FLAIR from source release tag to avoid broken conda bam2Bed12 entrypoint
RUN pip install --no-cache-dir git+https://github.com/BrooksLabUCSC/flair.git@v3.0.0

# Hotfix FLAIR v3.0.0 bug in ssPrep: juncsToBed12 can return None tuples for
# malformed corrected junctions; downstream code calls min(sizes) and crashes.
# Skip those reads instead of aborting the whole run.
RUN python - <<'PY'
from pathlib import Path
ssprep = Path('/opt/conda/envs/flair_env/lib/python3.12/site-packages/flair/ssPrep.py')
txt = ssprep.read_text()
needle = "        blocks, sizes, starts = juncsToBed12(bedObj.name, bedObj.chrom, bedObj.start, bedObj.end, newJuncs)\n\n        if correctStrand:\n"
patch = (
    "        blocks, sizes, starts = juncsToBed12(bedObj.name, bedObj.chrom, bedObj.start, bedObj.end, newJuncs)\n"
    "        if sizes is None:\n"
    "            # malformed corrected junction chain; skip read instead of crashing\n"
    "            continue\n\n"
    "        if correctStrand:\n"
)
if needle not in txt:
    raise SystemExit('Could not find ssPrep patch location')
ssprep.write_text(txt.replace(needle, patch, 1))
print('Patched', ssprep)
PY

# Provide working bam2Bed12 (FLAIR v3 package ships a dangling entrypoint)
RUN cat > /opt/conda/envs/flair_env/bin/bam2Bed12 <<'PY' && \
    chmod +x /opt/conda/envs/flair_env/bin/bam2Bed12
#!/opt/conda/envs/flair_env/bin/python
import argparse
import re
import sys
import pysam


def juncs_to_bed12(start, end, coords):
    sizes = []
    starts = []
    if len(coords) > 0:
        for num, junc in enumerate(coords):
            ss1, ss2 = junc
            if num == 0:
                st = 0
                size = ss1 - start
            else:
                st = coords[num - 1][1] - start
                size = ss1 - (st + start)
            starts.append(st)
            sizes.append(size)
        st = coords[-1][1] - start
        size = end - (st + start)
        starts.append(st)
        sizes.append(size)
        return len(starts), sizes, starts
    return 1, [end - start], [0]


def infer_mm2_junc_strand(read):
    orientation = read.flag
    try:
        junc_dir = read.get_tag("ts")
    except Exception:
        junc_dir = None

    if not junc_dir:
        left, right = read.cigar[0], read.cigar[-1]
        s1, s2 = read.query_sequence[:50], read.query_sequence[-50:]
        if ("T" * 10 in s1 and left[0] == 4 and left[1] >= 10) and (
            "A" * 10 in s2 and right[0] == 4 and right[1] >= 10
        ):
            junc_dir = "ambig"
        elif "T" * 10 in s1 and left[0] == 4 and left[1] >= 10:
            junc_dir = "-" if orientation == 16 else "+"
        elif "A" * 10 in s2 and right[0] == 4 and right[1] >= 10:
            junc_dir = "+" if orientation == 16 else "-"
        else:
            junc_dir = "ambig"
    else:
        if orientation == 0 and junc_dir == "+":
            junc_dir = "+"
        elif orientation == 0 and junc_dir == "-":
            junc_dir = "-"
        elif orientation == 16 and junc_dir == "+":
            junc_dir = "-"
        elif orientation == 16 and junc_dir == "-":
            junc_dir = "+"
    return junc_dir


def bam2bed12(alignment_file, keep_supplementary=False):
    strand_ok = {0: "+", 16: "-"}
    supplementary_ok = {2048: "+", 2064: "-"}
    positive_txn = "27,158,119"
    negative_txn = "217,95,2"
    unknown_txn = "99,99,99"

    bam = pysam.AlignmentFile(alignment_file, "rb")
    for read in bam.fetch():
        if read.flag not in strand_ok:
            if not (keep_supplementary and read.flag in supplementary_ok):
                continue

        qname = read.query_name
        chrom = read.reference_name
        ref_end = read.pos
        junctions = []
        for op, length in read.cigar:
            if op not in (0, 2, 3, 7, 8):
                continue
            if op == 3:
                junctions.append((ref_end, ref_end + length))
            ref_end += length

        start_pos = int(read.pos)
        end_pos = ref_end
        junc_dir = infer_mm2_junc_strand(read)
        flags = str(read.flag)
        score = read.mapq

        rgb = unknown_txn
        if junc_dir == "+":
            rgb = positive_txn
        elif junc_dir == "-":
            rgb = negative_txn
        else:
            junc_dir = "+" if flags == "0" else "-"

        blocks, sizes, starts = juncs_to_bed12(start_pos, end_pos, junctions)
        # BED12 blockSizes/blockStarts should be comma-separated lists without
        # an extra trailing comma for parser compatibility.
        sizes_str = ",".join(str(x) for x in sizes)
        starts_str = ",".join(str(x) for x in starts)
        outitems = [
            chrom,
            start_pos,
            end_pos,
            qname,
            score,
            junc_dir,
            start_pos,
            end_pos,
            rgb,
            blocks,
            sizes_str,
            starts_str,
        ]
        print("\t".join(str(x) for x in outitems))


def main():
    parser = argparse.ArgumentParser(description="Convert BAM to BED12 for FLAIR.")
    parser.add_argument("-i", "--input_bam", required=True, help="Input BAM file")
    parser.add_argument(
        "--keep_supplementary",
        action="store_true",
        default=False,
        help="Keep supplementary alignments",
    )
    args = parser.parse_args()
    bam2bed12(args.input_bam, keep_supplementary=args.keep_supplementary)
    return 0

if __name__ == "__main__":
    sys.argv[0] = re.sub(r"(-script\\.pyw?|\\.exe)?$", "", sys.argv[0])
    sys.exit(main())
PY

# Verify installations
RUN command -v flair && \
    command -v bam2Bed12 && \
    command -v minimap2 && \
    command -v samtools && \
    command -v bedtools


WORKDIR /usr/local/src
RUN wget https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz && \
    tar xvf gffread-0.12.7.Linux_x86_64.tar.gz && \
    mv gffread-0.12.7.Linux_x86_64/gffread /usr/local/bin


COPY flair-runner.py /usr/local/bin
