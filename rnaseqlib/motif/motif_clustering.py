##
## Utilities for clustering motifs
##
import os
import sys
import time

import pandas

import numpy as np

import cogent
from cogent.align.algorithm import sw_align

import rnaseqlib
import rnaseqlib.stats.stats_utils as stats_utils
import rnaseqlib.stats.clustering as clustering
import rnaseqlib.utils as utils


def get_pwm_from_clustalw(clustalw_fname):
    """
    Get PWM from CLUSTALW alignments file.

    Return PWM and motif object.
    """
    from Bio import Motif
    from Bio.Alphabet import IUPAC
    import Bio.Seq as bio_seq
    import Bio.AlignIO as align_io
    # Load CLUSTALW file
    if not os.path.isfile(clustalw_fname):
        raise Exception, "CLUSTALW file %s does not exist" %(clustalw_fname)
    clustalw_input = align_io.read(clustalw_fname, "clustal")
    motif_obj = Motif.Motif(alphabet=IUPAC.unambiguous_dna)
    # Add sequences from CLUSTALW alignment to motif object
    for clustalw_seq in clustalw_input.get_all_seqs():
        curr_seq = bio_seq.Seq(str(clustalw_seq.seq), IUPAC.unambiguous_dna)
        motif_obj.add_instance(curr_seq)
    # Compute PWM
    pwm = motif_obj.pwm()
    return pwm, motif_obj


def pwm_to_csv(pwm, output_fname,
               delim="\t",
               alphabet_order=["A", "C", "G", "T"]):
    """
    Output PWM object (BioPython PWM object from
    Motif class) to a text file in tab-delimited format.

    Each row is an alphabet letter and each column a position.
    """
    matrix_pwm = []
    for letter in alphabet_order:
        # Frequencies for all positions for this letter
        freqs = [float(row[letter]) for row in pwm]
        matrix_pwm.append(freqs)
    matrix_pwm = np.array(matrix_pwm)
    # Renormalize to 1 with lower precision    
    matrix_pwm = matrix_pwm / sum(matrix_pwm)
    output_file = open(output_fname, "w")
    for row in matrix_pwm:
        row_line = "%s\n" %(delim.join([str(f) for f in row]))
        output_file.write(row_line)
    output_file.close()


class MotifCluster:
    def __init__(self, kmers):
        """
        Store a list of kmers.
        """
        self.kmers = kmers
        # Compute kmer len
        self.kmer_len = None
        # Levenshtein distance function
        self.edit_dist_func = \
            lambda x, y: stats_utils.leven_dist(x[0], y[0])
        #self.edit_dist_func = \
        #    lambda x, y: stats_utils.leven_dist(x, y)


    def cluster_by_seq(self, method="sw"):
        """
        Cluster the motifs by sequence.
        """
        result = None
        if method == "sw":
            # Smith-Waterman clustering
            result = self.cluster_by_sw()
        return result

            
    def cluster_by_sw(self):
        """
        Cluster sequences pairwise by Smith-Waterman alignment.

        Returns distance matrix.
        """
        # Make pdist matrix with ij entry corresponding
        # to alignment between sequence i and sequence j
        score_matrix = []
        for kmer_i in self.kmers:
            score_row = []
            for kmer_j in self.kmers:
                alignment = sw_align(kmer_i, kmer_j,
                                     return_score=True)
                sw_score = alignment[1]
                score_row.append(sw_score)
            score_matrix.append(score_row)
        score_matrix = np.array(score_matrix)
        return score_matrix


    def cluster_by_edit(self, kmers, linkage_method):
        """
        Cluster sequences by edit distances.

        Parameters:
        -----------
        kmers : flat list of kmers
        linkage_method : determines linkage function for hierarchical
        clustering ('average', 'single', ...).

        """
        # Nest kmers to create matrix for clustering purposes.
        kmers = [kmers]
        data = np.array(kmers)
        hclust = clustering.hierarchical_clust(data,
                                               self.edit_dist_func,
                                               linkage_method)
        return hclust


def output_global_alignment(kmers_fname, output_dir,
                            alignment_program="mafft",
                            overwrite_alignment=False,
                            params={}):
    """
    Output a global alignment (*.aln) for a set of kmers.

    Using clustawl for now.

    Parameters:
    -----------
    kmers_fname : filename of FASTA file containing kmers
    output_dir : output directory
    """
    def get_clustal_cmd(fasta_fname, output_fname,
                        params={}):
        """
        Return clustal alignment command.
        """
        clustalw_cmd = \
            "clustalw -INFILE=%s -OUTFILE=%s -PIM" %(fasta_fname,
                                                     output_fname)
        for p in params:
            clustalw_cmd += " -%s=%s" %(p, params[p])
        return clustalw_cmd
    def get_mafft_cmd(fasta_fname, output_fname,
                      params={}):
        """
        Return mafft alignment command.
        """
        mafft_cmd = \
            "mafft --auto %s > %s" %(fasta_fname, output_fname)
        return mafft_cmd
    utils.make_dir(output_dir)
    output_fname = \
        os.path.join(output_dir,
                     "%s.aln" %(os.path.basename(kmers_fname)))
    if os.path.isfile(output_fname) and (not overwrite_alignment):
        print "Alignment filename %s exists. Skipping..." \
              %(output_fname)
    t1 = time.time()
    if alignment_program == "clustal":
        cmd = get_clustal_cmd(kmers_fname, output_fname,
                              params=params)
    elif alignment_program == "mafft":
        cmd = get_mafft_cmd(kmers_fname, output_fname,
                            params=params)
    else:
        raise Exception, "Unknown alignment program %s" \
              %(alignment_program)
    ret_val = os.system(cmd)
    if ret_val != 0:
        print "Call to %s failed!" %(alignment_program)
        print "Params were: ", str(params)
        raise Exception, "Alignment failed."
    t2 = time.time()
    print "Global alignment took %.2f minutes." %((t2 - t1)/60.)
    return output_fname
    

def main():
    kmers = ["TGTAT", "CGTAT", "TTAGT", "TCTAT", "TCTAC"]
    motif_clust = MotifCluster(kmers)
    #result = motif_clust.cluster_by_seq()
    #print "Result: ", result
    data = np.transpose(np.array(kmers))
    print "DATA: "
    print data
    dist_func = lambda x, y: sw_align(x, y, return_score=True)
    linkage_method = "average"
    hclust = clustering.hierarchical_clust(np.array(kmers),
                                           dist_func,
                                           linkage_method)


if __name__ == "__main__":
    main()
