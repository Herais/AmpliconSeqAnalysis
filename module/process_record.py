import re
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import motifs
import subprocess
import pickle



DPATH_TOOLS = "/content/tools"
DPATH_bowtie2 = f"{DPATH_TOOLS}/bowtie2-2.4.2-sra-linux-x86_64"
DPATH_NGmerge=f"{DPATH_TOOLS}/NGmerge"
DPATH_fastqjoin=f"{DPATH_TOOLS}/fastq-join"
DPATH_samtools=f"{DPATH_TOOLS}/samtools-1.17"
DPATH_clustalw=f"{DPATH_TOOLS}/clustalw-2.1"


# load local variables
from fio import fio

class Process_Record(object):
 
    def __init__(self, size:int=20):
        """
        Parameters
        ----------
        size : float
            radiu
        """
        super(Process_Record, self).__init__()

    @staticmethod
    def execute_shell_code(code):
        try:
            completed_process = subprocess.run(code, shell=True, check=True,
                                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                            universal_newlines=True)
            # Get the output and error messages
            output = completed_process.stdout
            error = completed_process.stderr

            # Process the output or error messages as needed
            if output:
                print("Output:")
                print(output)
            if error:
                print("Error:")
                print(error)

        except subprocess.CalledProcessError as e:
            print(f"An error occurred: {e}")

    @staticmethod
    def convert_quality_to_probability_of_error(quality:str):
        # gives a list of probabilty of error
        # 0.05 means a 5% chance that the base is called erroneously
        ls_qs = [10**((ord(q)-33)/-10) for q in list(quality)]
        return ls_qs
    
    @staticmethod
    def clean_ends(x, ls5, ls3, trunc3=False):
        ret = {}
        ret['seq'] = x
        ret['beg'] = 0
        ret['end'] = len(x)
        ret['bracket5'] = ''
        ret['bracket3'] = ''
        for f5 in ls5:
            search = re.search(f5, x)
            if search:
                ret['seq'] = x[search.span()[0]:]
                ret['beg'] = search.span()[0]
                ret['bracket5'] = f5
            for f3 in ls3:
                search = re.search(f3, x)
                if search:
                    ret['bracket3'] = f3
                if search.span()[0] > ret['beg']:
                    ret['seq'] = ret['seq'][:search.span()[1]+1]
                    ret['end'] = search.span()[1]+1
                if trunc3:
                    ret['seq'] = ret['seq'][:len(ret['seq'])-len(ret['seq'])%3] # clean the ends to end in triplicates
                    ret['end'] = len(ret['seq'])-len(ret['seq'])%3

        return ret
    
    @staticmethod
    def match_by_brackets(x, ls5, ls3):
        for f5 in ls5:
            if re.search(f5, x):
                for f3 in ls3:
                    if re.search(f3, x):
                        return True
        return False
    
    @staticmethod
    def create_SeqRecord_from_Nmerge(row):
        seq = Seq(row['seg_for_alignment']) # 'sequence'
        id = str(row['name'])
        name = str(row['name'])
        quality = row['quality']
        record = SeqRecord(seq,id, name, quality)
        return record
    
    @staticmethod
    def Nmerge_R1R2(
        fp_R1:str, 
        fp_R2:str, 
        fp_stitched:str, 
        fp_unstitched:str, 
        fp_adapter:str, 
        fp_R1R2_log:str,
        fp_stats:str,
        P_error_mean_max:float=0.01,
        savepickle:bool=False):
  
        ret = {}
        ret['fp_R1'] = fp_R1
        ret['fp_R2'] = fp_R2
        ret['fp_stitched'] = fp_stitched
        ret['fp_unstitched'] = fp_unstitched
        ret['fp_adapter'] = fp_adapter
        ret['fp_R1R2_log'] = fp_adapter
        ret['fp_stats'] = fp_stats
        ret['fp_pickle'] = fp_stitched.replace('.fastq', '_dict.pkl')

        str_command_NGmerge = "{}/NGmerge".format(DPATH_NGmerge)
        shell_code = str_command_NGmerge + ' '
        shell_code += '-1 {} '.format(fp_R1)
        shell_code += '-2 {} '.format(fp_R2)
        shell_code += '-o {} '.format(fp_stitched)
        shell_code += '-f {} '.format(fp_unstitched)
        shell_code += '-c {} '.format(fp_adapter)
        shell_code += '-l {} '.format(fp_R1R2_log)
        shell_code += '-g '
        shell_code += '-d '
        Process_Record.execute_shell_code(code=shell_code)

        ret['stats'] = ''
        # df_R1
        df_R1 = fio.process_fastq_to_df(fp_R1)
        num_all_reads = df_R1.shape[0]
        del df_R1

        # df_R1R2_stitched
        ret['stitched'] = fio.process_fastq_to_df(fp_stitched)

        ret['stats'] += 'stitched R1R2 / All Read Pairs = {}/{} = {}'.format(
                            ret['stitched'].shape[0],
                            num_all_reads,
                            ret['stitched'].shape[0] / num_all_reads
                        )
        ret['stats'] += "\n"

        """
        # df_R1R2_log
        df_R1R2_log = pd.read_csv(fp_R1R2_log, sep='\t')
        df_R1R2_stitch_properties = df_R1R2_log[~df_R1R2_log['OverlapLen'].isna()]
        df_R1R2_stitch_properties = df_R1R2_stitch_properties.rename(columns={'Read': 'name'})
        ret['log'] = df_R1R2_log

        # Merge Stitched Properties from df_R1R2_log to ret['stitched']
        df_R1R2_log = pd.read_csv(fp_R1R2_log, sep='\t')
        df_R1R2_stitch_properties = df_R1R2_log[~df_R1R2_log['OverlapLen'].isna()]
        df_R1R2_stitch_properties = df_R1R2_stitch_properties.rename(columns={'Read': 'name'})
        ret['stitched'] = ret['stitched'].merge(df_R1R2_stitch_properties, on='name', how='left')
        num_stitched = ret['stitched'].shape[0]
        """

        ls_filter2stat = []
        filter2stat = {}
        filter2stat['description'] = 'R1R2_stitched'
        filter2stat['n_input'] = num_all_reads
        filter2stat['n_input_unique'] = np.NaN
        filter2stat['n_output'] = ret['stitched'].shape[0]
        filter2stat['n_output_unique'] = np.NaN
        ls_filter2stat.append(filter2stat)



        # Filter by P_error_mean_max
        filter2stat = {}
        filter2stat['description'] = 'P_error_mean_max_of_stitch'
        filter2stat['n_input'] = ret['stitched'].shape[0]
        filter2stat['n_input_unique'] = np.NaN

        ret['stitched'] = ret['stitched'].loc[ret['stitched']['P_error_mean'] < P_error_mean_max, :]
        num_stitched_filtered = ret['stitched'].shape[0]

        ret['stats'] += 'stitched R1R2 filtered by P_error_mean < {} / All Read Pairs = {}/{} = {}'.format(
                            P_error_mean_max,
                            ret['stitched'].shape[0],
                            num_all_reads,
                            ret['stitched'].shape[0] / num_all_reads
                        )
        ret['stats'] += "\n"

        filter2stat['n_output'] = ret['stitched'].shape[0]
        filter2stat['n_output_unique'] = np.NaN
        ls_filter2stat.append(filter2stat)

        """
        # df_unstitched_R1
        fp_unstitched_1 = "{}_1.fastq".format(fp_unstitched.rstrip('.fastq'))
        fp_unstitched_2 = "{}_2.fastq".format(fp_unstitched.rstrip('.fastq'))
        ret['unstitched_R1'] = process_fastq_to_df(fp_unstitched_1)
        ret['unstitched_R2'] = process_fastq_to_df(fp_unstitched_2)

        ret['stats'] += 'unstitched R1 / All R1 = {}/{} = {}'.format(
                        ret['unstitched_R1'].shape[0],
                        num_all_reads,
                        ret['unstitched_R1'].shape[0] / num_all_reads
                    ) 
        ret['stats'] += "\n"
        """
    
        # unique sequence
        ret['seq'] = ret['stitched'].groupby('sequence')['name'].count().rename('count').reset_index().sort_values('count', ascending=False)
        ret['stats'] += '# unique seq {}'.format(ret['seq'].shape[0]) 
        ret['stats'] += "\n"
        ret['filters'] = ls_filter2stat

        with open(fp_stats, 'w') as f:
            f.write(ret['stats'])
            f.close()

        if savepickle:
            with open(ret['fp_pickle'], 'wb') as f:
                pickle.dump(ret, f)
            f.close()

        return ret