version 1.0

task FilterTranscripts {
    input {
        String fasta_path
        String gtf_path
        String expr_file_path
        String output_gtf_path
        Float threshold
        Int memoryGB
        Int diskSizeGB
    }

    command <<<
        python -c "
import os
from Bio import SeqIO
from Bio.Seq import Seq
import glob

def load_genomic_sequences(fasta_path):
    genomic_sequences = {}
    fasta_path = os.path.expanduser(fasta_path)
    try:
        for record in SeqIO.parse(fasta_path, 'fasta'):
            genomic_sequences[record.id] = record.seq
    except FileNotFoundError:
        print(f'Error: The file {fasta_path} was not found.')
    except Exception as e:
        print(f'An error occurred while loading genomic sequences: {e}')
    return genomic_sequences

def extract_post_transcript_sequence(chromosome, end_position, strand, genomic_sequences, length=20):
    if strand not in {'+', '-'}:
        raise ValueError('Strand must be \'+\' or \'-\'')
    if chromosome not in genomic_sequences:
        raise ValueError(f'Chromosome {chromosome} not found in genomic sequences.')
    sequence = genomic_sequences[chromosome]
    sequence_length = len(sequence)
    if strand == '+':
        start = end_position
        end = min(start + length, sequence_length)
    else:
        start = max(0, end_position - length)
        end = min(end_position, sequence_length)
    extracted_sequence = sequence[start:end]
    if strand == '-':
        extracted_sequence = extracted_sequence.reverse_complement()
    return str(extracted_sequence)

def check_sequence_for_pattern(sequence, strand):
    target = 'A' if strand == '+' else 'T'
    return sequence.count(target) >= 12 or target*6 in sequence

def calculate_total_transcripts_in_gtf(gtf_path):
    transcript_ids = set()
    try:
        with open(gtf_path, 'r') as gtf_file:
            for line in gtf_file:
                if line.startswith('#') or line.strip() == '':
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                attributes = fields[8]
                transcript_id = [attr for attr in attributes.split(';') if 'transcript_id' in attr]
                if not transcript_id:
                    continue
                transcript_id = transcript_id[0].split('\"')[1]
                transcript_ids.add(transcript_id)
    except FileNotFoundError:
        print(f'Error: The file {gtf_path} was not found.')
    except Exception as e:
        print(f'An error occurred while calculating total transcripts: {e}')
    total_transcripts = len(transcript_ids)
    print(f'Total transcripts in {gtf_path}: {total_transcripts}')
    return total_transcripts

def analyze_gtf_and_count_transcripts(gtf_path, genomic_sequences, output_gtf_path, tpm, threshold):
    summary = {
        'initial_count': 0,
        'ip_count': 0,
        'mono_exonic_count': 0,
        'mono_exonic_ip_count': 0,
        'non_ip_mono_exonic_count': 0,
        'low_tpm_mono_exonic_non_ip_count': 0,
        'final_count': 0,
    }
    initial_transcripts = set()
    ip_transcripts = set()
    mono_exonic_transcripts = set()
    mono_exonic_ip_transcripts = set()
    transcript_lines = {}
    exon_counts = {}
    try:
        with open(gtf_path, 'r') as gtf_file:
            for line in gtf_file:
                if line.startswith('#') or line.strip() == '':
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                feature_type, attributes = fields[2], fields[8]
                chromosome, end_position, strand = fields[0], int(fields[4]), fields[6]
                transcript_id = [attr for attr in attributes.split(';') if 'transcript_id' in attr][0].split('\"')[1]
                if transcript_id not in transcript_lines:
                    transcript_lines[transcript_id] = []
                transcript_lines[transcript_id].append(line)
                if feature_type == 'transcript':
                    initial_transcripts.add(transcript_id)
                    post_transcript_sequence = extract_post_transcript_sequence(chromosome, end_position, strand, genomic_sequences)
                    if check_sequence_for_pattern(post_transcript_sequence, strand):
                        ip_transcripts.add(transcript_id)
                if feature_type == 'exon':
                    exon_counts[transcript_id] = exon_counts.get(transcript_id, 0) + 1
        for transcript_id, count in exon_counts.items():
            if count == 1:
                mono_exonic_transcripts.add(transcript_id)
                if transcript_id in ip_transcripts:
                    mono_exonic_ip_transcripts.add(transcript_id)
        non_ip_mono_exonic_count = len(mono_exonic_transcripts - ip_transcripts)
        low_tpm_mono_exonic_non_ip = {t for t in mono_exonic_transcripts - ip_transcripts if tpm.get(t, 0) < threshold}
        final_transcripts = initial_transcripts - mono_exonic_ip_transcripts - low_tpm_mono_exonic_non_ip
        with open(output_gtf_path, 'w') as output_file:
            for transcript_id in final_transcripts:
                for line in transcript_lines.get(transcript_id, []):
                    output_file.write(line)
    except FileNotFoundError:
        print(f'Error: The file {gtf_path} was not found.')
    except Exception as e:
        print(f'An error occurred while analyzing GTF: {e}')
    final_count = len(final_transcripts)
    summary.update({
        'initial_count': len(initial_transcripts),
        'ip_count': len(ip_transcripts),
        'mono_exonic_count': len(mono_exonic_transcripts),
        'mono_exonic_ip_count': len(mono_exonic_ip_transcripts),
        'non_ip_mono_exonic_count': non_ip_mono_exonic_count,
        'low_tpm_mono_exonic_non_ip_count': len(low_tpm_mono_exonic_non_ip),
        'final_count': final_count,
    })
    print(f'Filtered GTF file saved to: {output_gtf_path}')
    print(f'Initial Count: {len(initial_transcripts)}')
    print(f'Filtered Count (IP): {len(ip_transcripts)}')
    print(f'Mono-exonic Count: {len(mono_exonic_transcripts)}')
    print(f'Mono-exonic IP Count: {len(mono_exonic_ip_transcripts)}')
    print(f'Mono-exonic Non-IP Count: {non_ip_mono_exonic_count}')
    print(f'Mono-exonic Non-IP with < {threshold} TPM: {len(low_tpm_mono_exonic_non_ip)}')
    print(f'Total Final Transcripts (after all filters): {final_count}')

def load_expression_values(expr_file_path):
    expr_values = {}
    with open(expr_file_path, 'r') as file:
        next(file)  # Skip header
        for line in file:
            parts = line.strip().split()
            transcript_id = parts[1]
            all_reads = float(parts[3])
            if transcript_id in expr_values:
                expr_values[transcript_id] += all_reads
            else:
                expr_values[transcript_id] = all_reads
    return expr_values

def calculate_tpm(expr_values):
    total_reads = sum(expr_values.values())
    tpm = {transcript_id: (reads / total_reads * 1e6) for transcript_id, reads in expr_values.items()}
    return tpm

def process_files(gtf_path, expr_file_path, fasta_path, output_gtf_path, threshold=1.0):
    genomic_sequences = load_genomic_sequences(fasta_path)
    expr_values = load_expression_values(expr_file_path)
    tpm = calculate_tpm(expr_values)
    analyze_gtf_and_count_transcripts(gtf_path, genomic_sequences, output_gtf_path, tpm, threshold)

process_files('~{fasta_path}', '~{gtf_path}', '~{expr_file_path}', '~{output_gtf_path}', ~{threshold})
        "
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 30
        memory: "~{memoryGB} GiB"
        disks: "local-disk ~{diskSizeGB} HDD"
    }

    output {
        File filtered_gtf = output_gtf_path
    }
}


workflow TranscriptFiltering {
    input {
        String fasta_path
        String gtf_path
        String expr_file_path
        String output_gtf_path
        Float threshold = 1.0
        Int memoryGB = 32
        Int diskSizeGB = 1024
        String docker = "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:latest"
    }

    call FilterTranscripts {
        input:
            fasta_path = fasta_path,
            gtf_path = gtf_path,
            expr_file_path = expr_file_path,
            output_gtf_path = output_gtf_path,
            threshold = threshold,
            docker = docker,
            memoryGB = memoryGB,
            diskSizeGB = diskSizeGB
    }

    output {
        File filtered_gtf = FilterTranscripts.filtered_gtf
    }
}
