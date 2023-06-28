import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import motifs
import pandas as pd
import numpy as np


class Amplicon(object):
 
    def __init__(self):
        """
        Parameters
        ----------
        size : float
            radiu
        """
        super(Amplicon, self).__init__()

    @staticmethod
    def get_df_amplicon_with_count(
            dfseq, 
            primer_F5, 
            primer_F3
        ):
        """
        """

        dfseq = dfseq.copy()

        ret = {}
        ret['filters'] = []
        ret['stats'] = 'Stitched Sequence -> Amplicon bracketed by primers\n'
        ret['stats'] += 'Input: # sequences {}, # unique {}\n'.format(
            dfseq['count'].sum(), dfseq.shape[0],)

        filter2stat = {}
        filter2stat['description'] = 'amplicon'
        filter2stat['n_input'] = dfseq['count'].sum()
        filter2stat['n_input_unique'] = dfseq.shape[0]


        dfseq = dfseq[(dfseq['sequence'].str.contains(primer_F5)) &
                    (dfseq['sequence'].str.contains(primer_F3))
                    ].copy()
        dfseq['amplicon'] = dfseq['sequence'].apply(
        lambda x: clean_ends(x, [primer_F5], [primer_F3])['seq']).copy()

        dfseq = dfseq.groupby('amplicon')['count'].sum().rename('count').reset_index()
        dfseq = dfseq.sort_values('count', ascending=False)

        ret['stats'] += 'Output: # sequences {}, # unique {}\n'.format(
            dfseq['count'].sum(), dfseq.shape[0])

        filter2stat['n_output'] = dfseq['count'].sum()
        filter2stat['n_output_unique'] = dfseq.shape[0]
        ret['filters'].append(filter2stat)

        ret['df'] = dfseq.copy()
        print(ret['stats'])

        return ret
    
    @staticmethod
    def nt_to_aa(
            seq_nt:str, 
            orf_offset:int=0,
        )->str:
        """
        """

        left = orf_offset
        right = len(seq_nt) - orf_offset -(len(seq_nt)-orf_offset)%3 + 1
        seq_to_translate = seq_nt[left:right]
        aa = Seq(seq_to_translate).translate()

        return str(aa)
    
    @staticmethod
    def drop_seq_with_stop_codon(dfseq):

        dfseq = dfseq.copy()
        ret = {}
        ret['filters'] = []

        filter2stat = {}
        filter2stat['description'] = 'stop codon'

        ret['stats'] = 'Translated Sequences -> Drop Sequences with Stop Codon\n'
        ret['stats'] += 'Input: # sequences {}, # unique {}\n'.format(
        dfseq['count'].sum(), dfseq.shape[0],)

        filter2stat['n_input'] = dfseq['count'].sum()
        filter2stat['n_input_unique'] = dfseq.shape[0]

        dfseq = dfseq[~dfseq['amplicon_aa'].str.contains('\*')]
        ret['stats'] += 'Output: # sequences {}, # unique {}\n'.format(
            dfseq['count'].sum(), dfseq.shape[0])

        filter2stat['n_output'] = dfseq['count'].sum()
        filter2stat['n_output_unique'] = dfseq.shape[0]

        ret['filters'].append(filter2stat)
        ret['df'] = dfseq
        print(ret['stats'])

        return ret