'''
LRQuant version 0.2

Usage:
    LRQuant.py -r <reads> -o <output_prefix>
    SIRV genome and annotation are set default
    for example:
    python3 LRQuant.py -r 1s_3dg.refined_without_polya.bam.fq -o 1s_3dg
'''
__version__ = "v0.2"
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from subprocess import call
from matplotlib_venn import venn3
import argparse
import shutil


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-V', '--version', action='version', version="%(prog)s ("+__version__+")")
    parser.add_argument(
        '--genome', '-g', type=str, default='',
        help='genome - fasta'
    )
    parser.add_argument(
        '--transcriptome', '-t', type=str,
        help='transcriptome - fasta, optional, will be made from other refrerence files if not provided', required = False, default='not_exist'
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
        help='output path'
    )
    parser.add_argument(
        '--output_path', '-o', type=str, default=os.getcwd()
,
        help='working directory'
    )
    return parser.parse_args()

def create_directory(dir_name):
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

def execute_commands(transcriptome, reads, genome, annotation, output_path):

    gffread_path = 'gffread'
    paftools_path = 'paftools.js'
    minimap2_path = 'minimap2'
    gff_compare_path = 'gffcompare'
    k8_path = 'k8'

    print('making juction.bed file')
    call(f'cd intermediate && {k8_path} {paftools_path} gff2bed {annotation} > junction.bed && cd ..', shell = True)    
    
    print('Calling minimap2 - mapped to genome sam')
    call(f'{minimap2_path} -uf -a -y -x splice:hq --junc-bed intermediate/junction.bed -t 64 {genome} {reads} > intermediate/splicing.mapping.sam', shell=True)
    
    print('Calling squanti3')
    call(f'convert_SAM_to_GTF_for_SQANTI3.py --sam_file intermediate/splicing.mapping.sam --output_prefix intermediate/splicing.mapping --reference_genome {genome} --allow_non_primary && cd {output_path}/intermediate && sqanti3_qc.py --force_id_ignore splicing.mapping.gtf {annotation} {genome} -o squanti3_OUT --skipORF --report skip --isoform_hits', shell = True)
    
    print('Calling gffcompare')
    call(f'cd intermediate && {gff_compare_path} splicing.mapping.gtf -o splicing.mapping.gff_compare -r {annotation}', shell = True)


#### LRQuant v0.0
def read_gffcompare_tmap(tmap_path):
    df = pd.read_csv(tmap_path, sep='\t')
    return df

def generate_expression_matrix(df):
    # Example: Consider transcripts with class code '=' for direct matches
    expression_data = df[df['class_code'].isin(['c','='])]
    expression_data = df

    # Aggregate transcripts by qry_gene_id (or ref_gene_id if you prefer)
    # and count them for expression levels. Adjust as needed for your analysis.
    expression_matrix = expression_data.groupby('ref_id').size().reset_index(name='expression')
    
    return expression_matrix

def save_expression_matrix(expression_matrix, output_path):
    expression_matrix.to_csv(output_path, sep='\t', index=False)


#### LRQuant later versions
def process_mapping_stats(file_path):
    mapping_stat = pd.read_csv(file_path, sep='\t', header=None)
    mapping_stat_filtered = mapping_stat[mapping_stat[1] != '*'].copy()
    mapping_stat_filtered['alignment_score1'] = mapping_stat_filtered[3].apply(lambda x: float(x.split(':')[2]))
    mapping_stat_filtered['alignment_score2'] = mapping_stat_filtered[5].apply(lambda x: float(x.split(':')[2]))
    return mapping_stat_filtered

def identify_assignments(mapping_stat_filtered):
    ambiguous_assignment = mapping_stat_filtered[mapping_stat_filtered[0].isin(mapping_stat_filtered[0][mapping_stat_filtered[0].duplicated()])]
    unique_assignment = mapping_stat_filtered[~mapping_stat_filtered[0].isin(mapping_stat_filtered[0][mapping_stat_filtered[0].duplicated()])]
    return unique_assignment, ambiguous_assignment

def save_read_ids(read_ids, file_path):
    pd.DataFrame(read_ids).to_csv(file_path, sep='\t', index=False, header=False)

def plot_alignment_scores(unique_assignment, ambiguous_assignment, file_prefix):
    file_prefix = "Figure"
    for score_column, fig_number in zip(['alignment_score1', 'alignment_score2'], [1, 2]):
        plt.hist(unique_assignment[score_column], bins=100, log=True)
        plt.savefig(f'qc_figures/{file_prefix}_{fig_number}_unique.png')
        plt.close()
        
        plt.hist(ambiguous_assignment[score_column], bins=100, log=True)
        plt.savefig(f'qc_figures/{file_prefix}_{fig_number}_ambiguous.png')
        plt.close()
        
        plt.hist(unique_assignment[score_column], bins=100, log=True, alpha=0.5, label='Unique')
        plt.hist(ambiguous_assignment[score_column], bins=100, log=True, alpha=0.5, label='Ambiguous')
        plt.legend()
        plt.savefig(f'qc_figures/{file_prefix}_{fig_number}_combined.png')
        plt.close()


def process_squanti_output(file_path):
    squanti = pd.read_csv(file_path, sep='\t')
    ambiguous_assignment_sq = squanti[squanti['Isoform'].isin(squanti['Isoform'][squanti['Isoform'].duplicated()])].copy()
    unique_assignment_sq = squanti[~squanti['Isoform'].isin(squanti['Isoform'][squanti['Isoform'].duplicated()])].copy()
    return unique_assignment_sq, ambiguous_assignment_sq

def plot_venn_diagrams(ambigous_reads_transcriptome, unique_reads_transcriptome, ambigous_reads_squanti, unique_reads_squanti):
    plt.figure(figsize=(8, 8))
    
    # Plot for ambiguous reads from transcriptome and SQANTI
    plt.subplot(1, 2, 1)  # 1 row, 2 columns, first plot
    venn3([set(ambigous_reads_transcriptome), set(unique_reads_transcriptome), set(ambigous_reads_squanti)],
          set_labels=('Ambiguous Transcriptome', 'Unique Transcriptome', 'Ambiguous SQANTI'))
    plt.title("Overlap: Ambiguous & Unique Reads from Transcriptome with Ambiguous SQANTI")
    
    # Adjust font size for clarity
    plt.gca().set_facecolor('white')
    plt.gca().patch.set_alpha(0)
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)

    # Plot for unique reads from transcriptome and SQANTI
    plt.subplot(1, 2, 2)  # 1 row, 2 columns, second plot
    venn3([set(ambigous_reads_transcriptome), set(unique_reads_transcriptome), set(unique_reads_squanti)],
          set_labels=('Ambiguous Transcriptome', 'Unique Transcriptome', 'Unique SQANTI'))
    plt.title("Overlap: Ambiguous & Unique Reads from Transcriptome with Unique SQANTI")
    
    # Adjust font size for clarity
    plt.gca().set_facecolor('white')
    plt.gca().patch.set_alpha(0)
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)
    
    plt.savefig('qc_figures/venn_diagrams.png')
    plt.close()


def output_overlap_read_ids(ambigous_reads_transcriptome, unique_reads_transcriptome, ambigous_reads_squanti, unique_reads_squanti):
    # Convert numpy arrays or lists to sets for efficient operation
    ambigous_reads_transcriptome_set = set(ambigous_reads_transcriptome)
    unique_reads_transcriptome_set = set(unique_reads_transcriptome)
    ambigous_reads_squanti_set = set(ambigous_reads_squanti)
    unique_reads_squanti_set = set(unique_reads_squanti)
    
    # Calculate intersections
    atas = ambigous_reads_transcriptome_set & ambigous_reads_squanti_set  # Ambiguous in both transcriptome and SQANTI
    atus = ambigous_reads_transcriptome_set & unique_reads_squanti_set   # Ambiguous in transcriptome but unique in SQANTI
    utas = unique_reads_transcriptome_set & ambigous_reads_squanti_set   # Unique in transcriptome but ambiguous in SQANTI
    utus = unique_reads_transcriptome_set & unique_reads_squanti_set     # Unique in both transcriptome and SQANTI
    
    # Save to files
    pd.DataFrame(list(atas)).to_csv('reads_id_extraction/atas_reads.txt', index=False, header=False)
    pd.DataFrame(list(atus)).to_csv('reads_id_extraction/atus_reads.txt', index=False, header=False)
    pd.DataFrame(list(utas)).to_csv('reads_id_extraction/utas_reads.txt', index=False, header=False)
    pd.DataFrame(list(utus)).to_csv('reads_id_extraction/utus_reads.txt', index=False, header=False)
        
    return  atas, atus, utas, utus


def filter_isoforms(df, match_type=None, diff_to_tts=None):
    filtered = df
    if match_type is not None:
        filtered = filtered[filtered['Matching_type'] == match_type]
    if diff_to_tts is not None:
        filtered = filtered[filtered['Diff_to_TTS'] == diff_to_tts]
    return filtered

def plot_diff_to_tts_histogram(data, range_vals=None, bins=100, fig_name="figure.png", output_dir="qc_figures"):
    plt.hist(data, bins=bins, range=range_vals)
    plt.savefig(f"{output_dir}/{fig_name}")
    plt.close()


def generate_all_plots(matrix_a_ism_only):
    conditions = [
        {'match_type': 'primary', 'diff_to_tts': 0, 'range_vals': (-50, 50), 'fig_name': 'figure3_1.png'},
        {'match_type': 'primary', 'diff_to_tts': 0, 'range_vals': None, 'fig_name': 'figure3_2.png', 'secondary_filter': True},
        {'match_type': 'secondary', 'diff_to_tts': 0, 'range_vals': None, 'fig_name': 'figure3_3.png', 'secondary_filter': True},
        {'match_type': 'primary', 'diff_to_tts': None, 'range_vals': (-10, 20), 'fig_name': 'figure3_4.png'},
        {'match_type': 'secondary', 'diff_to_tts': None, 'range_vals': (-10, 20), 'fig_name': 'figure3_5.png'},
        {'match_type': 'primary', 'diff_to_tts': None, 'range_vals': (-10000, 10000), 'fig_name': 'figure3_6.png'},
        {'match_type': 'secondary', 'diff_to_tts': None, 'range_vals': (-10000, 10000), 'fig_name': 'figure3_7.png'},
    ]
    for cond in conditions:
        filtered_data = filter_isoforms(matrix_a_ism_only, cond.get('match_type'), cond.get('diff_to_tts'))
        if cond.get('secondary_filter', False):
            filtered_data = filter_isoforms(filtered_data, 'secondary')
        data_to_plot = filtered_data['Diff_to_TTS'] if 'diff_to_tts' in cond or cond.get('secondary_filter', False) else filtered_data
        plot_diff_to_tts_histogram(data_to_plot, cond['range_vals'], fig_name=cond['fig_name'])


def main():
    args = parse_args()
    transcriptome = args.transcriptome
    reads = args.reads
    genome = args.genome
    annotation = args.annotation
    output_prefix = args.output_prefix
    output_path = args.output_path
                                                
                         
    directories = ['reads_id_extraction', 'lrquant_out', 'intermediate', 'qc_figures', 'gffcompare_out']
    for dir_name in directories:
        create_directory(dir_name)

    execute_commands(transcriptome, reads, genome, annotation, output_path)
    
    mapping_stat_filtered = process_mapping_stats('intermediate/transcript.mapping_stat.tsv')
    unique_assignment_tr, ambiguous_assignment_tr = identify_assignments(mapping_stat_filtered)
    
    ambigous_reads_transcriptome = ambiguous_assignment_tr[0].unique()
    unique_reads_transcriptome = unique_assignment_tr[0].values
    
    unique_assignment_sq, ambiguous_assignment_sq = process_squanti_output('intermediate/squanti3_OUT_isoform_hits.txt')
    ambigous_reads_squanti = ambiguous_assignment_sq['Isoform'].unique()
    unique_reads_squanti = unique_assignment_sq['Isoform'].values
    
    file_prefix = "analysis" 

    plot_alignment_scores(unique_assignment_tr, ambiguous_assignment_tr, file_prefix)
    plot_venn_diagrams(ambigous_reads_transcriptome, unique_reads_transcriptome, ambigous_reads_squanti, unique_reads_squanti)
    atas, atus, utas, utus = output_overlap_read_ids(ambigous_reads_transcriptome, unique_reads_transcriptome, ambigous_reads_squanti, unique_reads_squanti)
    

    ###unique assignment matrixes
    matrix1 = unique_assignment_tr
    matrix2 = unique_assignment_sq[unique_assignment_sq['Isoform'].isin(atus)]
    
    
    print('********************\n'*10)
    print('ambiguous transcriptome mapping but unique squanti support read total count\nthese reads are uniquely assigned')
    print(len(matrix2))

    ### print reads (that were mapped) to be discarded
    print('reads that falls into ambiguous transcriptome but does not have any squanti category at all \nthese are input for gffcompare to clean up')
    
    squanti = pd.read_csv('intermediate/squanti3_OUT_isoform_hits.txt', sep = '\t')
    
    gff_input_count = np.logical_not(pd.DataFrame(ambiguous_assignment_tr[0].unique())[0].isin(squanti['Isoform'])).sum()
    print(gff_input_count)

    gff_reads = pd.DataFrame(ambiguous_assignment_tr[0].unique())[0][np.logical_not(pd.DataFrame(ambiguous_assignment_tr[0].unique())[0].isin(squanti['Isoform']))]
    (gff_reads + '\t').to_csv('reads_id_extraction/ambiguous_tr_no_squanti_reads.txt', sep = '\n', index=False, header = False)

    ## combine unique assginements before gffcompare clean up
    u_assignment2 = pd.DataFrame(matrix2).groupby('Hit').count()['Isoform']
    u_assignment2 = pd.DataFrame({"Isoform":u_assignment2.index, "u_assignment2": u_assignment2})
    u_assignment1 = pd.DataFrame(matrix1).groupby(1).count()[0]
    u_assignment1 = pd.DataFrame({"Isoform":u_assignment1.index, "u_assignment1": u_assignment1})
    u_assignment = u_assignment1.merge(u_assignment2, how = 'outer', on = 'Isoform').fillna(0)
    u_assignment['count'] = u_assignment['u_assignment1'] + u_assignment['u_assignment2']
    u_assignment = u_assignment[['Isoform','count']]

    ### gffcompare clean up
    gff = pd.read_csv('intermediate/splicing.mapping.gff_compare.splicing.mapping.gtf.tmap', sep = '\t')
    gff = gff[gff['qry_gene_id'].isin(gff_reads)]

    print('gffcompare uniquely assinged read count')
    print(gff['class_code'].isin(['c','=']).sum())
    print('discarded reads as gffcompare could not rescue')
    print(gff_input_count-gff['class_code'].isin(['c','=']).sum())

    rescued = gff[gff['class_code'].isin(['c','='])]
    gff_next = gff[np.logical_not(gff['class_code'].isin(['c','=']))]
    rescued = pd.DataFrame(rescued.groupby('ref_id').count()['class_code'])
    rescued['Isoform'] = rescued.index
    (gff_next['qry_gene_id'] + '\t').to_csv('reads_id_extraction/gff_not_rescued_reads.txt', sep = '\n', index=False, header = False)

    ### combine all unique assignment
    expression_combine = u_assignment[['Isoform','count']].merge(rescued, on='Isoform', how = 'outer').fillna(0)
    expression_combine['count'] = expression_combine['count'] + expression_combine['class_code']
    expression_combine['transcript_name'] =  expression_combine['Isoform']
    expression_combine = expression_combine[['transcript_name','count']]


    ### ambigous competibility matrix
    matrix_a = ambiguous_assignment_sq[ambiguous_assignment_sq['Isoform'].isin(atas)]
    print('ambigous assigned from both transcriptome mapping and squanti read total count')
    print(len(matrix_a['Isoform'].unique()))

    #### ambiguous table to work out
    match = pd.DataFrame(matrix_a.groupby('Isoform')['Match'].unique())
    match['Isoform'] = match.index
    match["Match_string"] = list(map(lambda x: np.array2string(x), list(match['Match'])))

    ism_only = match['Isoform'][match['Match_string']=="['ISM']"]
    fsm_only = match['Isoform'][match['Match_string']=="['FSM']"]
    ism_fsm = match['Isoform'][match['Match_string']=="['ISM' 'FSM']"]

    (ism_only + '\t').to_csv('reads_id_extraction/ism_only_reads.txt', sep = '\n', index=False, header = False)
    (fsm_only + '\t').to_csv('reads_id_extraction/sm_only_reads.txt', sep = '\n', index=False, header = False)
    (ism_fsm + '\t').to_csv('reads_id_extraction/ism_fsm_reads.txt', sep = '\n', index=False, header = False)


    matrix_a_ism_only = matrix_a[matrix_a['Isoform'].isin(ism_only)]
    matrix_a_fsm_only = matrix_a[matrix_a['Isoform'].isin(fsm_only)]
    matrix_a_ism_fsm = matrix_a[matrix_a['Isoform'].isin(ism_fsm)]


    cm = pd.DataFrame(matrix_a_fsm_only)[['Hit','Isoform']]
    a_assignment = 1/cm.groupby('Isoform').count()
    a_assignment['Isoform'] = a_assignment.index
    a_assignment.index.names = ['index']
    a_assignment = a_assignment.merge(cm, how = 'outer', on = 'Isoform').fillna(0)
    a_assignment.columns = ["Expression", "Reads", "Isoform"]
    a_assignment = a_assignment.groupby('Isoform').sum('Expression')
    a_assignment["Isoform"] = a_assignment.index
    a_assignment.index.names = ['index']
    fsm_only_a_assignment = a_assignment
    matrix_a_ism_only.to_csv('reads_id_extraction/ism_only_reads.csv')

    generate_all_plots(matrix_a_ism_only)



    ####unique
    matrix_a_ism_only = matrix_a_ism_only[['Isoform','Hit','Match','Diff_to_TTS','Matching_type']]
    matrix_a_ism_only['Matching_type']=='primary'
    matrix_a_ism_only['Matching_type']=='secondary'
    matrix_a_ism_only.loc[:, 'Diff_to_TTS_within_10'] = abs(matrix_a_ism_only['Diff_to_TTS']) <= 10
    matrix_a_ism_only.loc[:, 'Diff_to_TTS_beyond_30'] = abs(matrix_a_ism_only['Diff_to_TTS']) > 30

    matrix_a_ism_only_primary_check = matrix_a_ism_only[matrix_a_ism_only['Matching_type']=='primary'][matrix_a_ism_only['Diff_to_TTS_within_10']]['Isoform'].unique()
    matrix_a_ism_only_secondary_check = matrix_a_ism_only[matrix_a_ism_only['Matching_type']=='secondary'][matrix_a_ism_only['Diff_to_TTS_beyond_30']]['Isoform'].unique()
    matrix_a_ism_only_assigned = matrix_a_ism_only[matrix_a_ism_only['Matching_type']=='primary']
    matrix_a_ism_only_assigned = matrix_a_ism_only_assigned[matrix_a_ism_only_assigned['Isoform'].isin(np.intersect1d(matrix_a_ism_only_primary_check,matrix_a_ism_only_secondary_check))]
    ism_only_a_assignment1 = pd.DataFrame(matrix_a_ism_only_assigned.groupby('Hit').count()['Match'])
    ism_only_a_assignment1['Isoform']=ism_only_a_assignment1.index

    #####ambiguous
    ism_only_a_assignment2 = matrix_a_ism_only[np.logical_not(matrix_a_ism_only['Isoform'].isin(ism_only_a_assignment1['Isoform']))]
    cm = pd.DataFrame(ism_only_a_assignment2)[['Hit','Isoform']]
    a_assignment = 1/cm.groupby('Isoform').count()
    a_assignment['Isoform'] = a_assignment.index
    a_assignment.index.names = ['index']
    a_assignment = a_assignment.merge(cm, how = 'outer', on = 'Isoform').fillna(0)
    a_assignment.columns = ["Expression", "Reads", "Isoform"]
    a_assignment = a_assignment.groupby('Isoform').sum('Expression')
    a_assignment["Isoform"] = a_assignment.index
    a_assignment.index.names = ['index']
    ism_only_a_assignment2 = a_assignment


    matrix_a_ism_fsm_FSM = matrix_a_ism_fsm[matrix_a_ism_fsm['Match'] == 'FSM']
    matrix_a_ism_fsm_FSM = matrix_a_ism_fsm_FSM[np.array(abs(matrix_a_ism_fsm_FSM['Diff_to_TSS'])<=2) & np.array(abs(matrix_a_ism_fsm_FSM['Diff_to_TSS'])<=2)]
    matrix_a_ism_fsm_FSM = matrix_a_ism_fsm_FSM[matrix_a_ism_fsm_FSM['Matching_type'] == 'primary']
    ism_fsm_a_assignment1 = pd.DataFrame(matrix_a_ism_fsm_FSM.groupby('Hit').count()['Match'])
    ism_fsm_a_assignment1['Isoform']=ism_fsm_a_assignment1.index
    print('reads assigned by stricker FSM rule under FSM-ISM')
    len(matrix_a_ism_fsm_FSM['Isoform'].unique())
    print('********************\n'*10)

    #####ambiguous
    ism_fsm_a_assignment2 = matrix_a_ism_fsm[np.logical_not(matrix_a_ism_fsm['Isoform'].isin(ism_fsm_a_assignment1['Isoform']))]
    cm = pd.DataFrame(ism_fsm_a_assignment2)[['Hit','Isoform']]
    a_assignment = 1/cm.groupby('Isoform').count()
    a_assignment['Isoform'] = a_assignment.index
    a_assignment.index.names = ['index']
    a_assignment = a_assignment.merge(cm, how = 'outer', on = 'Isoform').fillna(0)
    a_assignment.columns = ["Expression", "Reads", "Isoform"]
    a_assignment = a_assignment.groupby('Isoform').sum('Expression')
    a_assignment["Isoform"] = a_assignment.index
    a_assignment.index.names = ['index']
    ism_fsm_a_assignment2 = a_assignment


    mat = matrix_a_ism_fsm_FSM['Diff_to_TSS']
    plt.hist(mat[mat.isin(range(-10, 10))], bins=20)

    plt.savefig('qc_figures/figure4_1.png')
    plt.close()

    plt.hist(mat[mat.isin(range(-100, 101))], bins=300)

    plt.savefig('qc_figures/figure4_2.png')
    plt.close()

    mat = matrix_a_ism_fsm_FSM['Diff_to_TTS']
    plt.hist(mat[mat.isin(range(-10, 10))], bins=20)

    plt.savefig('qc_figures/figure4_3.png')
    plt.close()

    plt.hist(mat[mat.isin(range(-100, 101))], bins=300)

    plt.savefig('qc_figures/figure4_4.png')
    plt.close()


    expression_combine1 = expression_combine.rename(columns={"transcript_name": "Isoform"}).merge(ism_only_a_assignment1, on='Isoform', how = 'outer').fillna(0).merge(ism_only_a_assignment2, on='Isoform', how = 'outer').fillna(0)
    expression_combine2 = expression_combine1.merge(ism_fsm_a_assignment1, on='Isoform', how = 'outer').fillna(0).merge(ism_fsm_a_assignment2, on='Isoform', how = 'outer').fillna(0)
    expression_combine3 = expression_combine2.merge(fsm_only_a_assignment, on='Isoform', how = 'outer').fillna(0)
    expression_combine3['count']=expression_combine3['count']+expression_combine3['Match_x']+expression_combine3['Expression_x']+expression_combine3['Match_y']+expression_combine3['Expression_y']+expression_combine3['Expression']
    expression_combine = expression_combine3.rename(columns={"Isoform": "transcript_name"})[["transcript_name","count"]]


    ####final output for benchmarking comparison
    expression_combine[['transcript_name','count']].to_csv(f"lrquant_out/{output_prefix}.lrquant.tsv", sep = '\t', index= False)
    
 
    ## Output LRquant-gffcompare
    # Assuming gffcompare has been run and its output is available
    tmap_path = 'intermediate/splicing.mapping.gff_compare.splicing.mapping.gtf.tmap'
    expression_matrix_output_path = f'gffcompare_out/{output_prefix}.gffcompare.tsv'

    # Read the .tmap file
    tmap_df = read_gffcompare_tmap(tmap_path)

    # Generate the expression matrix based on gffcompare results
    expression_matrix = generate_expression_matrix(tmap_df)

    # Save the expression matrix to a file
    save_expression_matrix(expression_matrix, expression_matrix_output_path)


    print('moving all output folders into output directory')
    call(f'sudo rm -r {output_prefix}', shell=True)
    call(f'mkdir {output_prefix}', shell=True)
    call(f'mv -t {output_prefix} intermediate lrquant_out gffcompare_out qc_figures reads_id_extraction', shell=True)

if __name__ == '__main__':
    main()
