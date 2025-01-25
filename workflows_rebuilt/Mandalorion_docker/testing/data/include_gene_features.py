#!/usr/bin/env python3

from collections import defaultdict
import csv


# Debugging: Adjusting the parsing logic for attributes
def add_gene_and_transcript_rows_debug(input_path, output_path):
    gene_data = defaultdict(lambda: {"start": float("inf"), "end": 0})
    transcript_data = defaultdict(
        lambda: {"start": float("inf"), "end": 0, "gene_id": ""}
    )

    gene_id_to_transcript_ids = defaultdict(set)
    transcript_id_to_exon_rows = defaultdict(list)

    # Parse the GTF file and collect data for genes and transcripts
    with open(input_path, "r") as infile:
        reader = csv.reader(infile, delimiter="\t")
        exon_rows = []

        for row in reader:
            if len(row) < 9 or row[2] != "exon":
                continue  # Skip invalid rows or non-exon features

            chrom, source, feature, start, end, score, strand, frame, attributes = row

            # Safely parse the attributes field
            attributes_dict = {}
            for item in attributes.split(";"):
                if not item.strip():
                    continue
                parts = item.strip().split(" ", 1)
                if len(parts) == 2:
                    key, value = parts
                    attributes_dict[key.strip()] = value.strip('"')

            gene_id = attributes_dict.get("gene_id")
            transcript_id = attributes_dict.get("transcript_id")
            start, end = int(start), int(end)

            # Update gene and transcript bounds
            if gene_id:
                gene_data[gene_id]["start"] = min(gene_data[gene_id]["start"], start)
                gene_data[gene_id]["end"] = max(gene_data[gene_id]["end"], end)
                gene_data[gene_id]["chrom"] = chrom
                gene_data[gene_id]["strand"] = strand

            if transcript_id:
                transcript_data[transcript_id]["start"] = min(
                    transcript_data[transcript_id]["start"], start
                )
                transcript_data[transcript_id]["end"] = max(
                    transcript_data[transcript_id]["end"], end
                )
                transcript_data[transcript_id]["gene_id"] = gene_id
                transcript_data[transcript_id]["chrom"] = chrom
                transcript_data[transcript_id]["strand"] = strand

                gene_id_to_transcript_ids[gene_id].add(transcript_id)

            # exon_rows.append(row)
            transcript_id_to_exon_rows[transcript_id].append(row)

    # Write the updated GTF with new gene and transcript rows
    with open(output_path, "w", newline="") as outfile:
        writer = csv.writer(
            outfile, delimiter="\t", quoting=csv.QUOTE_NONE, escapechar="\\"
        )

        # Add gene rows
        for gene_id, bounds in gene_data.items():
            writer.writerow(
                [
                    bounds["chrom"],
                    "AddedFeature",
                    "gene",
                    bounds["start"],
                    bounds["end"],
                    ".",
                    bounds["strand"],
                    ".",
                    f'gene_id "{gene_id}"; gene_name "{gene_id}";',
                ]
            )

            for transcript_id in gene_id_to_transcript_ids[gene_id]:
                bounds = transcript_data[transcript_id]

                writer.writerow(
                    [
                        bounds["chrom"],
                        "AddedFeature",
                        "transcript",
                        bounds["start"],
                        bounds["end"],
                        ".",
                        bounds["strand"],
                        ".",
                        f'gene_id "{bounds["gene_id"]}"; transcript_id "{transcript_id}"; gene_name "{gene_id}";',
                    ]
                )

                exon_rows = transcript_id_to_exon_rows[transcript_id]

                # Add original exon rows
                for row in exon_rows:
                    row[-1] += f' gene_name "{gene_id}";'
                    writer.writerow(row)


# Retry processing the input file with the debug function
input_file_path = "SIRVs1-7.annot.reduced.gtf"
output_file_path = "ladeda.gtf"
add_gene_and_transcript_rows_debug(input_file_path, output_file_path)
