import os
import argparse
from subprocess import call
import pandas as pd

__version__ = '0.0'

def parse_args():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-V', '--version', action='version', version="%(prog)s ("+__version__+")")
    parser.add_argument(
        '--genome', '-g', type=str, default='',
        help='genome - fasta'
    )
    parser.add_argument(
        '--annotation', '-a', type=str, default='',
        help='annotation - gtf'
    )
    parser.add_argument(
        '--reads', '-r', type=str,
        help='reads - fasta/fastq'
    )
    parser.add_argument(
        '--output_prefix', '-p', type=str, default='LRQuant_OUT',
        help='output prefix'
    )
    parser.add_argument(
        '--output_path', '-o', type=str, default=os.getcwd(),
        help='working directory'
    )
    parser.add_argument(
        '--threads', '-t', type=int, default=4,
        help='threads'
    )
    return parser.parse_args()
    return parser.parse_args()

def create_directory(dir_name):
    """
    Create a directory if it does not exist.
    """
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

def execute_commands(args):
    """
    Execute a series of commands for transcriptome analysis.
    """
    # Define paths to tools
    gffread_path = '/opt/conda/envs/mas_seq/bin/gffread'
    paftools_path = '/home/jupyter/projects/morf/analysis/reads/Quantification_bechmarking/k8-0.2.4/minimap2/misc/paftools.js'
    minimap2_path = '/home/jupyter/projects/morf/analysis/reads/Quantification_bechmarking/minimap2/minimap2'
    gff_compare_path = '/opt/conda/envs/mas_seq/bin/gffcompare'
    k8_path = '/home/jupyter/projects/iso_quantification_benchmark/mouse_simulation/k8-0.2.4/node-v18.17.0/k8/k8'

    # Make junction.bed file
    print('Making junction.bed file')
    call(f'cd intermediate && {k8_path} {paftools_path} gff2bed {args.annotation} > junction.bed && cd ..', shell=True)    
    
    # Call minimap2
    print('Calling minimap2 - mapped to genome sam')
    call(f'{minimap2_path} -uf -a -y -x splice:hq --junc-bed intermediate/junction.bed -t {args.threads} {args.genome} {args.reads} > intermediate/splicing.mapping.sam', shell=True)
    
    # Call convert_SAM_to_GTF
    print('Calling convert_SAM_to_GTF')
    call('docker pull us-central1-docker.pkg.dev/methods-dev-lab/lrtools-sqanti3/lrtools-sqanti3-plus', shell=True)
    call(f'docker run --rm -v /:/mnt/mydata/ --entrypoint=/bin/sh us-central1-docker.pkg.dev/methods-dev-lab/lrtools-sqanti3/lrtools-sqanti3-plus -c "cd /mnt/mydata{args.output_path} && convert_SAM_to_GTF_for_SQANTI3.py --sam_file intermediate/splicing.mapping.sam --output_prefix intermediate/splicing.mapping --reference_genome /mnt/mydata{args.genome} --allow_non_primary && cd /mnt/mydata{args.output_path}/intermediate"', shell=True)
         
    # Call gffcompare
    print('Calling gffcompare')
    call(f'cd intermediate && {gff_compare_path} splicing.mapping.gtf -o splicing.mapping.gff_compare -r {args.annotation}', shell=True)


def read_gffcompare_tmap(tmap_path):
    """
    Read the gffcompare tmap file into a pandas DataFrame.
    """
    df = pd.read_csv(tmap_path, sep='\t')
    return df

def generate_expression_matrix(df):
    """
    Generate an expression matrix from the DataFrame.
    """
    expression_data = df[df['class_code'].isin(['c','='])]
    expression_matrix = expression_data.groupby('ref_id').size().reset_index(name='expression')
    return expression_matrix

def save_expression_matrix(expression_matrix, output_path):
    """
    Save the expression matrix to a file.
    """
    expression_matrix.to_csv(output_path, sep='\t', index=False)

def main():
    """
    Main function to run the analysis.
    """
    args = parse_args()

    create_directory('intermediate')
    execute_commands(args)

    tmap_path = 'intermediate/splicing.mapping.gff_compare.splicing.mapping.gtf.tmap'
    expression_matrix = generate_expression_matrix(read_gffcompare_tmap(tmap_path))
    save_expression_matrix(expression_matrix, f'{args.output_prefix}_expression_matrix.tsv')

    print(f'Expression matrix saved to {args.output_prefix}_expression_matrix.tsv')

if __name__ == '__main__':
    main()
