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
    
    @staticmethod
    def bracket_by_NC(
            df, 
            ls_N_terminus, 
            ls_C_terminus
        ):
        """
        """
        
        df = df.copy()

        ret = {}
        ret['filters'] = []

        filter2stats = {}
        filter2stats['description'] = 'NC linkage'
        filter2stats['n_input'] = df['count'].sum()
        filter2stats['n_input_unique'] = df.shape[0]

        ret['stats'] = '-> Narrow down to sequences with NC linkages\n'
        ret['stats'] += 'Input: # sequences {}, # unique {}\n'.format(
            filter2stats['n_input'], filter2stats['n_input_unique'])

        df = df[df.amplicon_aa.apply(lambda x: match_by_brackets(x, ls_N_terminus, ls_C_terminus))]
        filter2stats['n_output'] = df['count'].sum()
        filter2stats['n_output_unique'] = df.shape[0]
        ret['stats'] += 'Output: # sequences {}, # unique {}\n'.format(
            filter2stats['n_output'], filter2stats['n_output_unique'])

        ret['filters'].append(filter2stats)

        df['beg'] = df['amplicon_aa'].apply(lambda x: clean_ends(x, ls_N_terminus, ls_C_terminus)['beg'])
        df['end'] = df['amplicon_aa'].apply(lambda x: clean_ends(x, ls_N_terminus, ls_C_terminus)['end']).copy()
        df['seq_aa_NtoC'] = df['amplicon_aa'].apply(lambda x: clean_ends(x, ls_N_terminus, ls_C_terminus)['seq'])
        df['Nterm'] = df['amplicon_aa'].apply(lambda x: clean_ends(x, ls_N_terminus, ls_C_terminus)['bracket5'])
        df['Cterm'] = df['amplicon_aa'].apply(lambda x: clean_ends(x, ls_N_terminus, ls_C_terminus)['bracket3'])
        df['CpxA N-term linkage'] = df['Nterm'].apply(lambda x: Npat2link[x])
        df['CpxA C-term linkage'] = df['Cterm'].apply(lambda x: Cpat2link[x])

        # with and without histidine
        df['H'] = 1
        df.loc[~df['amplicon_aa'].str.contains('H'), 'H'] = 0

        ret['df'] = df.copy()
        print(ret['stats'])

        return ret
    
    @staticmethod
    def assign_ref_id(
            row, 
            dfref
        ):
        """
        """

        cols = ['s#', 'CpxA N-term linkage', 'CpxA C-term linkage']
        dfref_wH = dfref[:49][cols]
        dfref_woH = dfref[49:-1][cols]
        dfref_wt = dfref[-1:][cols]
        Nt = row['CpxA N-term linkage']
        Ct = row['CpxA C-term linkage']

        if row['H']:
            cols = ['CpxA N-term linkage', 'CpxA C-term linkage']
            NC2sn_wH = dfref_wH.set_index(cols).to_dict()['s#']
            NC2sn_wt = dfref_wt.set_index(cols).to_dict()['s#']
            if (Nt, Ct) in NC2sn_wH:
            return NC2sn_wH[(Nt, Ct)]
            elif (Nt, Ct) in NC2sn_wt:
            return NC2sn_wt[(Nt, Ct)]
            else:
            return -1
        elif not row['H']:
            cols = ['CpxA N-term linkage', 'CpxA C-term linkage']
            NC2sn_woH = dfref_woH.set_index(cols).to_dict()['s#']
            if (Nt, Ct) in NC2sn_woH:
            return NC2sn_woH[(Nt, Ct)]
            else:
            return -2
        return -3
    
    @staticmethod
    def assign_nc(dfref):
        """
        """

        dfref = dfref.copy()

        dfref['N_seqpat'] = dfref['seq_aa'].apply(lambda x: x[26:34])
        dfref['C_seqpat'] = dfref['seq_aa'].apply(lambda x: x[-16-8:-16+1])

        Npat2link = dfref[['CpxA N-term linkage', 'N_seqpat']].drop_duplicates().set_index('N_seqpat').to_dict()['CpxA N-term linkage']
        Cpat2link = dfref[['CpxA C-term linkage', 'C_seqpat']].drop_duplicates().set_index('C_seqpat').to_dict()['CpxA C-term linkage']
        
        dfref['NC_linkage'] = dfref[['N_seqpat', 'C_seqpat']].apply(lambda x: (x[0], x[1]), axis=1)

        ls_N_terminus = dfref['N_seqpat'].unique()
        ls_C_terminus = dfref['C_seqpat'].unique()
        print('N\n', '\n'.join(list(ls_N_terminus)))
        print('C\n', '\n'.join(list(ls_C_terminus)))

        return dfref.copy()
