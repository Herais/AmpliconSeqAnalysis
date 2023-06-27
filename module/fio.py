import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import motifs
import pandas as pd
import numpy as np
import os
import gzip
import shutil
import re

# local imports
from process_record import Process_Record

class fio(object):
 
    def __init__(self, size:int=20):
        """
        Parameters
        ----------
        size : float
            radiu
        """
        super(io, self).__init__()

    @staticmethod
    def decompress_file(input_file, output_file):
        with gzip.open(input_file, 'rb') as gz_file:
            with open(output_file, 'wb') as output:
                shutil.copyfileobj(gz_file, output)

    @staticmethod
    def process_fastq_lines(lines=None):
        ks = ['name', 'sequence', 'optional', 'quality']
        return {k: v for k, v in zip(ks, lines)}
    
    @staticmethod
    def process_fastq_to_df(path_fastq:str):
        records = []
        n = 4
        with open(path_fastq, 'r') as fh:
            lines = []
            for line in fh:
                lines.append(line.rstrip())
                if len(lines) == n:
                    record = fio.process_fastq_lines(lines)
                    #sys.stderr.write("Record: %s\n" % (str(record)))
                    lines = []
                    records.append(record)
        fh.close()
        
        df = pd.DataFrame(records)
        df['len_nt'] = df['sequence'].apply(len)
        df['P_error'] = df['quality'].apply(Process_Record.convert_quality_to_probability_of_error)
        df['P_error_mean'] = df['P_error'].apply(np.mean)
        df['name'] = df['name'].apply(lambda x: x.lstrip('@'))
        

        cols_ordered = ['name', 'sequence', 'len_nt', 'P_error_mean', 'P_error', 'optional', 'quality']
        df = df[cols_ordered]

        return df.copy()
    
    @staticmethod
    def create_record_nt(row):
        seq = Seq(row['seq_nt'])
        id = str(row['s#'])
        name = "v{}_N{}_C{}".format(row['v#'], row['CpxA N-term linkage'], row['CpxA C-term linkage'])
        description = row['Binding design']
        record = SeqRecord(seq,id, name, description)
        return record

    @staticmethod
    def create_record_aa(row):
        seq = Seq(row['seq_aa'])
        id = str(row['s#'])
        name = "v{}_N{}_C{}".format(row['v#'], row['CpxA N-term linkage'], row['CpxA C-term linkage'])
        description = row['Binding design']
        record = SeqRecord(seq,id, name, description)
        return record

    @staticmethod
    def generate_paths_for_Nmerge(
        batch_id:str,
        dirpath_source:str='/content/ngs',
    )->dict:

        """
        Sample Usage
        ----
        batch_id = '30-879258766'
        name2fp = generate_paths_for_Nmerge(batch_id)
        """

        for dirpath, dirnames, filenames in os.walk(
            "{}/{}/00_fastq/".format(dirpath_source, batch_id)):
            break

        name2fp = {}
        for fn in filenames:
            fp = "{}/{}/00_fastq/{}".format(dirpath_source, batch_id, fn.rstrip('.gz'))
            if re.search('.gz', fn):
                fio.decompress_file(fp + '.gz', fp)

            if re.search('.fastq', fn):
                # name_lib
                name_ngs = fn.split('_')[0]
                if name_ngs not in name2fp:
                    name2fp[name_ngs] = {}

                # fp to fastq file
                if re.search('_R1', fn):
                    name2fp[name_ngs]['R1'] = {}
                    name2fp[name_ngs]['R1']['fastq'] = fp

                if re.search('_R2', fn):
                    name2fp[name_ngs]['R2'] = {}
                    name2fp[name_ngs]['R2']['fastq'] = fp

            # fp Nmerge
            name2fp[name_ngs]['Nmerge'] = {}

            ls = ['fastq', 'unstitched', 'adapter', 'log', 'stats', 'pkl']
            for item in ls:
                fp = "{}/{}/00_fastq/{}_{}_{}".format(dirpath_source, batch_id, name_ngs, 'Nmerge', item)
                if item == 'fastq': fp += '.fastq'
                elif item == 'log': fp += '.log'
                elif item == 'pkl': fp += '.pkl'
                name2fp[name_ngs]['Nmerge'][item] = fp

        return name2fp
