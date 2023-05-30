import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import motifs

class io(object):
 
    def __init__(self, size:int=20):
        """
        Parameters
        ----------
        size : float
            radiu
        """
        super(io, self).__init__()

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
                    record = process_fastq_lines(lines)
                    #sys.stderr.write("Record: %s\n" % (str(record)))
                    lines = []
                    records.append(record)
        fh.close()
        
        df = pd.DataFrame(records)
        df['len_nt'] = df['sequence'].apply(len)
        df['P_error'] = df['quality'].apply(convert_quality_to_probabilityoferror)
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

