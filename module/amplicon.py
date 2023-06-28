import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import motifs
import pandas as pd
import numpy as np
import re


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
    # Functions for sequence truncation
    def clean_ends(
            x, 
            ls5, 
            ls3, 
            trunc3=False
        ):
        """
        """
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
                #return ret
        return ret

    @staticmethod
    def match_by_brackets(
            x, 
            ls5, 
            ls3
        ):
        """
        """
        for f5 in ls5:
            if re.search(f5, x):
                for f3 in ls3:
                    if re.search(f3, x):
                        return True
        return False

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
        lambda x: Amplicon.clean_ends(x, [primer_F5], [primer_F3])['seq']).copy()

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
            ls_C_terminus,
            dfref,
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

        df = df[df.amplicon_aa.apply(lambda x: Amplicon.match_by_brackets(x, ls_N_terminus, ls_C_terminus))]
        filter2stats['n_output'] = df['count'].sum()
        filter2stats['n_output_unique'] = df.shape[0]
        ret['stats'] += 'Output: # sequences {}, # unique {}\n'.format(
            filter2stats['n_output'], filter2stats['n_output_unique'])

        ret['filters'].append(filter2stats)

        df['beg'] = df['amplicon_aa'].apply(lambda x: Amplicon.clean_ends(x, ls_N_terminus, ls_C_terminus)['beg'])
        df['end'] = df['amplicon_aa'].apply(lambda x: Amplicon.clean_ends(x, ls_N_terminus, ls_C_terminus)['end']).copy()
        df['seq_aa_NtoC'] = df['amplicon_aa'].apply(lambda x: Amplicon.clean_ends(x, ls_N_terminus, ls_C_terminus)['seq'])
        df['Nterm'] = df['amplicon_aa'].apply(lambda x: Amplicon.clean_ends(x, ls_N_terminus, ls_C_terminus)['bracket5'])
        df['Cterm'] = df['amplicon_aa'].apply(lambda x: Amplicon.clean_ends(x, ls_N_terminus, ls_C_terminus)['bracket3'])
        
        dfref = Amplicon.assign_nc(dfref=dfref)
        Npat2link = dfref[['CpxA N-term linkage', 'N_seqpat']].drop_duplicates().set_index('N_seqpat').to_dict()['CpxA N-term linkage']
        Cpat2link = dfref[['CpxA C-term linkage', 'C_seqpat']].drop_duplicates().set_index('C_seqpat').to_dict()['CpxA C-term linkage']        
        
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

    @staticmethod
    def graph_abudance_bar(
            df_seq,
            primer_F5,
            primer_F3,
            ls_patN,
            ls_patC,
            dfref,
            orf_shift=1,
        ):
        """
        """

        df_L = df_seq[df_seq['sequence'].str.contains(primer_F5)]
        df_LR = df_L[df_L['sequence'].str.contains(primer_F3)]
        df_LR['amplicon'] = df_LR['sequence'].apply(lambda x: Amplicon.clean_ends(x, [primer_F5], [primer_F3])['seq'])
        df_LRU = df_LR.groupby('amplicon')['count'].sum().rename('count').reset_index()
        df_LRU['translated'] = df_LRU['amplicon'].apply(lambda x: str(Seq(x[orf_shift:len(x)-(len(x)-1)%3]).translate()))

        df_LRU_nostar = df_LRU[~df_LRU['translated'].str.contains('\*')]

        df_LRU_nostar_NC = df_LRU_nostar[df_LRU_nostar['translated'].apply(lambda x: Amplicon.match_by_brackets(x, ls_patN, ls_patC))]
        df_LRU_nostar_NC['beg'] = df_LRU_nostar_NC['translated'].apply(lambda x: Amplicon.clean_ends(x, ls_patN, ls_patC)['beg'])
        df_LRU_nostar_NC['end'] = df_LRU_nostar_NC['translated'].apply(lambda x: Amplicon.clean_ends(x, ls_patN, ls_patC)['end'])

        dfref = Amplicon.assign_nc(dfref=dfref)
        Npat2link = dfref[['CpxA N-term linkage', 'N_seqpat']].drop_duplicates().set_index('N_seqpat').to_dict()['CpxA N-term linkage']
        Cpat2link = dfref[['CpxA C-term linkage', 'C_seqpat']].drop_duplicates().set_index('C_seqpat').to_dict()['CpxA C-term linkage']
        
        df_LRU_nostar_NC['Nterm'] = df_LRU_nostar_NC['translated'].apply(lambda x: Amplicon.clean_ends(x, ls_patN, ls_patC)['bracket5'])
        df_LRU_nostar_NC['Cterm'] = df_LRU_nostar_NC['translated'].apply(lambda x: Amplicon.clean_ends(x, ls_patN, ls_patC)['bracket3'])
        df_LRU_nostar_NC['CpxA N-term linkage'] = df_LRU_nostar_NC['Nterm'].apply(lambda x: Npat2link[x])
        df_LRU_nostar_NC['CpxA C-term linkage'] = df_LRU_nostar_NC['Cterm'].apply(lambda x: Cpat2link[x])

        return df_LRU

    @staticmethod
    def process_merged_to_df(
            ret_merged, 
            dfref, 
            primer_F5, 
            primer_F3,
        ):
        """
        """
        track_filters = []
        track_filters.extend(ret_merged['filters'])

        ret = Amplicon.get_df_amplicon_with_count(dfseq=ret_merged['seq'],
                                primer_F5=primer_F5,
                                primer_F3=primers_F3,
                                )
        track_filters.extend(ret['filters'])

        df = ret['df']
        df['amplicon_aa'] = df['amplicon'].apply(lambda x: Amplicon.nt_to_aa(x, 1))

        # narrow to sequences w/o stop codon
        ret = Amplicon.drop_seq_with_stop_codon(df)
        track_filters.extend(ret['filters'])
        df = ret['df']

        # narrow to sequences w/ defined NC terminus
        ls_N_terminus = dfref['N_seqpat'].unique()
        ls_C_terminus = dfref['C_seqpat'].unique()
        ret = Amplicon.bracket_by_NC(df, ls_N_terminus, ls_C_terminus)
        track_filters.extend(ret['filters'])
        df = ret['df']

        # assign s# identity
        df['s#'] = df.apply(lambda x: Amplicon.assign_ref_id(x, dfref), axis=1)

        # adjust NGS amplicon size to size in ref
        df['len_amplicon_aa'] = df['amplicon_aa'].apply(len)
        df['len_seq_aa_NtoC'] = df['seq_aa_NtoC'].apply(len)
        df['len_amplicon_aa_adj'] = df['len_amplicon_aa'] +3-9 # NGS is 3 less at N, 9 more at

        return df.copy()