import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import motifs

import matplotlib as mp
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import display, HTML
import matplotlib.font_manager as fm

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
    def nt_to_aa_3orf(row, ls_pat):
        for pat in ls_pat:
            m = re.search(pat, row['amplicon_aa_0'])
            if m: 
                return  row['amplicon_aa_0']

            m = re.search(pat, row['amplicon_aa_1'])
            if m: 
                return  row['amplicon_aa_1']

            m = re.search(pat, row['amplicon_aa_2'])
            if m: 
                return row['amplicon_aa_2']
        return np.NAN
    
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
        #print('N\n', '\n'.join(list(ls_N_terminus)))
        #print('C\n', '\n'.join(list(ls_C_terminus)))

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
    def Pipeline_001(
            ret_merged, 
            dfref, 
            primer_F5,
            primer_F3,
            ls_sensor_pat,
        ):
        """
        """
        track_filters = []
        track_filters.extend(ret_merged['filters'])

        ret = Amplicon.get_df_amplicon_with_count(dfseq=ret_merged['seq'],
                                primer_F5=primer_F5,
                                primer_F3=primer_F3,
                                )
        track_filters.extend(ret['filters'])

        df = ret['df']
        df['amplicon_aa_0'] = df['amplicon'].apply(lambda x: Amplicon.nt_to_aa(x,0))
        df['amplicon_aa_1'] = df['amplicon'].apply(lambda x: Amplicon.nt_to_aa(x,1))
        df['amplicon_aa_2'] = df['amplicon'].apply(lambda x: Amplicon.nt_to_aa(x,2))
        S_aa = df.apply(lambda x: Amplicon.nt_to_aa_3orf(x, ls_sensor_pat), axis=1)
        df['amplicon_aa'] = S_aa

        


        filter = {}
        filter['description'] = 'amplicon_orf'
        filter['n_input'] = df['count'].sum()
        filter['n_input_unique'] = df.shape[0]
        toprint = 'amplicon_orf -> Drop Sequences w/o 5wlk sensor domain\n'
        toprint += 'Input: # sequences {}, # unique {}\n'.format(
                    filter['n_input'], filter['n_input_unique'])

        df = df[~df['amplicon_aa'].isna()]
        filter['n_output'] = df['count'].sum()
        filter['n_output_unique'] = df.shape[0]
        toprint += 'Output: # sequences {}, # unique {}\n'.format(
            filter['n_input'], filter['n_input_unique'])
        track_filters.extend([filter])
        print(toprint)


        # narrow to sequences w/o stop codon
        ret = Amplicon.drop_seq_with_stop_codon(df)
        track_filters.extend(ret['filters'])
        df = ret['df']

        # narrow to sequences w/ defined NC terminus
        ls_N_terminus = dfref['N_seqpat'].unique()
        ls_C_terminus = dfref['C_seqpat'].unique()
        ret = Amplicon.bracket_by_NC(df, ls_N_terminus, ls_C_terminus, dfref)
        track_filters.extend(ret['filters'])
        df = ret['df']

        # assign s# identity
        df['s#'] = df.apply(lambda x: Amplicon.assign_ref_id(x, dfref), axis=1)

        # adjust NGS amplicon size to size in ref
        df['len_amplicon_aa'] = df['amplicon_aa'].apply(len)
        df['len_seq_aa_NtoC'] = df['seq_aa_NtoC'].apply(len)
        df['len_amplicon_aa_adj'] = df['len_amplicon_aa'] +3-9 # NGS is 3 less at N, 9 more at

        return df.copy(), track_filters

    @staticmethod
    def Pipeline_002(
            ret_merged, 
            dfref, 
            primer_F5,
            primer_F3,
            ls_sensor_pat,
        ):
        """
        """
        track_filters = []
        track_filters.extend(ret_merged['filters'])

        ret = Amplicon.get_df_amplicon_with_count(dfseq=ret_merged['seq'],
                                primer_F5=primer_F5,
                                primer_F3=primer_F3,
                                )
        track_filters.extend(ret['filters'])

        df = ret['df']
        df['amplicon_aa_0'] = df['amplicon'].apply(lambda x: Amplicon.nt_to_aa(x,0))
        df['amplicon_aa_1'] = df['amplicon'].apply(lambda x: Amplicon.nt_to_aa(x,1))
        df['amplicon_aa_2'] = df['amplicon'].apply(lambda x: Amplicon.nt_to_aa(x,2))
        S_aa = df.apply(lambda x: Amplicon.nt_to_aa_3orf(x, ls_sensor_pat), axis=1)
        df['amplicon_aa'] = S_aa

        filter = {}
        filter['description'] = 'amplicon_orf'
        filter['n_input'] = df['count'].sum()
        filter['n_input_unique'] = df.shape[0]

        df = df[~df['amplicon_aa'].isna()]
        filter['n_output'] = df['count'].sum()
        filter['n_output_unique'] = df.shape[0]
        track_filters.append(filter)


        # narrow to sequences w/o stop codon
        ret = Amplicon.drop_seq_with_stop_codon(df)
        track_filters.extend(ret['filters'])
        df = ret['df']

        return df.copy(), track_filters
    
    @staticmethod
    def barplot_abundance(
            df,
            xlabel='strain id',
            ylabel='occurrences',
            title=None,
            color='skyblue',
            fontsize=12,
            alpha=0.7,
            rot=0,
            figsize=(16,2),
        ):
        """
        """

        df = df.copy()
        fg, ax = plt.subplots(figsize=figsize)
        S_plot = df.groupby(['s#'])['count'].sum().rename('count').astype('int')
        S_plot.plot(
            kind='bar',
            ax=ax,
            xlabel=xlabel,
            ylabel=ylabel,
            title=title,
            fontsize=fontsize,
            color=color,
            alpha=alpha,
            rot=rot,
        )
        ax.set_xticks(ticks=np.arange(S_plot.shape[0]),
                                        labels=S_plot.index.astype('int')
                                    )
        plt.show()

    @staticmethod
    def heatmap(data, row_labels, col_labels, ax=None,
                cbar_kw=None, cbarlabel="", **kwargs):
        """
        Create a heatmap from a numpy array and two lists of labels.

        Parameters
        ----------
        data
            A 2D numpy array of shape (M, N).
        row_labels
            A list or array of length M with the labels for the rows.
        col_labels
            A list or array of length N with the labels for the columns.
        ax
            A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
            not provided, use current axes or create a new one.  Optional.
        cbar_kw
            A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
        cbarlabel
            The label for the colorbar.  Optional.
        **kwargs
            All other arguments are forwarded to `imshow`.
        """

        if ax is None:
            ax = plt.gca()

        if cbar_kw is None:
            cbar_kw = {}

        # Plot the heatmap
        im = ax.imshow(data, **kwargs)

        # Create colorbar
        cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
        cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

        # Show all ticks and label them with the respective list entries.
        ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
        ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)

        # Let the horizontal axes labeling appear on top.
        ax.tick_params(top=True, bottom=False,
                    labeltop=True, labelbottom=False)

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
                rotation_mode="anchor")

        # Turn spines off and create white grid.
        ax.spines[:].set_visible(False)

        ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
        ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
        ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
        ax.tick_params(which="minor", bottom=False, left=False)

        return im, cbar

    @staticmethod
    def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                        textcolors=("black", "white"),
                        threshold=None, **textkw):
        """
        A function to annotate a heatmap.

        Parameters
        ----------
        im
            The AxesImage to be labeled.
        data
            Data used to annotate.  If None, the image's data is used.  Optional.
        valfmt
            The format of the annotations inside the heatmap.  This should either
            use the string format method, e.g. "$ {x:.2f}", or be a
            `matplotlib.ticker.Formatter`.  Optional.
        textcolors
            A pair of colors.  The first is used for values below a threshold,
            the second for those above.  Optional.
        threshold
            Value in data units according to which the colors from textcolors are
            applied.  If None (the default) uses the middle of the colormap as
            separation.  Optional.
        **kwargs
            All other arguments are forwarded to each call to `text` used to create
            the text labels.
        """

        if not isinstance(data, (list, np.ndarray)):
            data = im.get_array()

        # Normalize the threshold to the images color range.
        if threshold is not None:
            threshold = im.norm(threshold)
        else:
            threshold = im.norm(data.max())/2.

        # Set default alignment to center, but allow it to be
        # overwritten by textkw.
        kw = dict(horizontalalignment="center",
                verticalalignment="center")
        kw.update(textkw)

        # Get the formatter in case a string is supplied
        if isinstance(valfmt, str):
            valfmt = mp.ticker.StrMethodFormatter(valfmt)

        # Loop over the data and create a `Text` for each "pixel".
        # Change the text's color depending on the data.
        texts = []
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
                text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
                texts.append(text)

        return texts

    @staticmethod
    def generate_heatmap_matrix(df):
        df = df.copy()

        cols = ['CpxA C-term linkage', 'CpxA N-term linkage']
        df_NC_crossed =  df.groupby(cols)['count'].sum().unstack()
        cols_ordered = sorted(df_NC_crossed.columns.tolist(), key=lambda x: x[1:])
        df_NC_crossed = df_NC_crossed[cols_ordered]
        rows_ordered = sorted((df_NC_crossed.index), key=lambda x: x[1:])
        df_NC_crossed = df_NC_crossed.reindex(rows_ordered)
        df_NC_crossed.fillna(0, inplace=True)
        df_NC_crossed = df_NC_crossed.astype('int')

        df_NC_crossed_unique =  df.groupby(cols)['count'].count().unstack()
        cols_ordered = sorted(df_NC_crossed_unique.columns.tolist(), key=lambda x: x[1:])
        df_NC_crossed_unique = df_NC_crossed_unique[cols_ordered]
        rows_ordered = sorted((df_NC_crossed_unique.index), key=lambda x: x[1:])
        df_NC_crossed_unique = df_NC_crossed_unique.reindex(rows_ordered)
        df_NC_crossed_unique.fillna(0, inplace=True)
        df_NC_crossed_unique = df_NC_crossed_unique.astype('int')

        df_NC_mut = df_NC_crossed_unique / df_NC_crossed

        return df_NC_crossed, df_NC_crossed_unique, df_NC_mut
    
    @staticmethod
    def heatmap_NC_abundance(
            df,
            title='Abundance N x C Linkages Lib',
            figsize=(8,5),
            cmap="GnBu",
            alpha=0.7,
            valfmt="{x:.0f}",
            cbarlabel='count',
            fontsize=16,
        ):
        """
        """
        df = df.copy()
        fig, ax = plt.subplots(figsize=figsize,)
        im, cbar = Amplicon.heatmap(
            data=df,
            row_labels=list(df.index),
            col_labels=list(df.columns),
            ax=ax,
            cmap=cmap, #"YlGn"
            alpha=alpha,
            cbar_kw=None,
            cbarlabel=cbarlabel,
        )
        texts = Amplicon.annotate_heatmap(
            im,
            valfmt=valfmt,
        )
        ax.set_title(title)
        fig.tight_layout()

        plt.grid(False)
        plt.axis("on")
        plt.show()