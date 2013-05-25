##
## SpliceGraph class.
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.events.parseTables as parseTables

import gffutils

import collections
from collections import OrderedDict, namedtuple

class OrderedDefaultdict(collections.OrderedDict):
    """
    An ordered, default dictionary.
    """
    def __init__(self, *args, **kwargs):
        if not args:
            self.default_factory = None
        else:
            if not (args[0] is None or callable(args[0])):
                raise TypeError('first argument must be callable or None')
            self.default_factory = args[0]
            args = args[1:]
        super(self.__class__, self).__init__(*args, **kwargs)

    def __missing__ (self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = default = self.default_factory()
        return default

    def __reduce__(self):  # optional, for pickle support
        args = self.default_factory if self.default_factory else tuple()
        return type(self), args, None, None, self.items()


class SpliceEdges:
    """
    Splice site edges.
    """
    def __init__(self, strand):
        self.strand = strand
        self.edges = OrderedDefaultdict(list)


    def count_from(self, node):
        """
        Count edges from node.
        """
        pass


    def count_to(self, node):
        """
        Count edges to node.
        """
        pass


    def get_edges_from(self, node):
        """
        Get edges from a given node.
        """
        return self.edges[node]


    def add_edge(self, source, sink):
        self.edges[source].append(sink)


    def __getitem__(self, node):
        return self.get_edges_from(node)

        

class Acceptors(SpliceEdges):
    """
    Acceptors to donors class.
    """
    def __init__(self, strand):
        SpliceEdges.__init__(self, strand)


    def __str__(self):
        return "Acceptors(%s)" %(self.edges)


    def __repr__(self):
        return self.__str__()
    


class Donors(SpliceEdges):
    """
    Donors class.
    """
    def __init__(self, strand):
        SpliceEdges.__init__(self, strand)
        self.strand = None


    def __str__(self):
        return "Donors(%s)" %(self.edges)


    def __repr__(self):
        return self.__str__()
    
    


def define_RI(sg,
              gff_out,
              min_intron_len=10,
              multi_iso=False):
    """
    Define retained introns.

    [ A ]---[ B ]
    [     C     ]

    Look at B's donors (in this case A) and see if there is
    a known exonic unit (A's start spliced to B's end) that
    that matches that, in this case C.
    

    Parameters:
    -----------

    sg : SpliceGraph
    min_intron_len : minimum intron length
    """
    def get_intron_len(donor, acceptor):
        """
        Calculate intron length.
        """
        if donor.strand == "+":
            intron_len = (acceptor.start_coord - 1) - \
                         (donor.end_coord + 1) + 1
        else:
            # Minus strand -- reconsider length given that
            # donor is "first" in transcript space compared to
            # acceptor, and that start > end
            intron_len = (donor.end_coord - 1) - \
                         (acceptor.start_coord + 1) + 1
        return intron_len
    if multi_iso:
        raise Exception, "Multiple isoforms not supported."
    # Keep track of observed retained introns
    retained_introns = {}
    for strand in sg.acceptors_to_donors:
        for acceptor in sg.acceptors_to_donors[strand].edges:
            # If any of the donor units to this acceptor
            # have the acceptor end as their end coordinate, it's a retained
            # intron event
            donors = sg.acceptors_to_donors[strand].edges[acceptor]
            for donor in donors:
                # If there's a node that has this acceptor end as its end coordinate
                # and the donor's start as the start coordinate, then it's a
                # retained intron
                intron_as_exon = Unit(donor.start, acceptor.end)
                if sg.has_node(intron_as_exon):
                    # Get length of intron
                    intron_len = get_intron_len(donor, acceptor)
                    if intron_len < min_intron_len:
                        continue
                    ri = (donor.coords_str, acceptor.coords_str)
                    if ri in retained_introns:
                        continue
                    # Output retained intron to gff file
                    output_RI(gff_out, donor, acceptor, intron_len)
#                    print "It's a retained intron between %s and %s: %d" \
#                          %(donor, acceptor, intron_len)
                    # Output retained intron
                    retained_introns[ri] = True


def output_RI(gff_out, donor, acceptor, intron_len,
              source="RI"):
    """
    Output a retained intron event.
    """
    chrom = donor.chrom
    strand = donor.strand
    donor_start = donor.start_coord
    donor_end = donor.end_coord
    acceptor_start = acceptor.start_coord
    acceptor_end = acceptor.end_coord
    ri_name = "%s@%s" %(donor.coords_str,
                        acceptor.coords_str)
    gene_start = donor_start
    gene_end = acceptor_end
    # For GFF record purposes, ensure start < end always
    if gene_start > gene_end:
        gene_start, gene_end = gene_end, gene_start
    gene_rec = gffutils.Feature(seqid=chrom,
                                source=source,
                                featuretype="gene",
                                start=gene_start,
                                end=gene_end,
                                strand=strand,
                                attributes={"ID": [ri_name],
                                            "Name": [ri_name]})
    # Output mRNA containing the retained intron and then output its exons
    # First output retained intron using "withRI" suffix
    long_mRNA_name = "%s.A" %(ri_name)
    # Long mRNA record has same start/end as gene record
    long_mRNA_rec = gffutils.Feature(seqid=chrom,
                                     source=source,
                                     featuretype="mRNA",
                                     start=gene_start,
                                     end=gene_end,
                                     strand=strand,
                                     attributes={"ID": [long_mRNA_name],
                                                 "Parent": [ri_name]})
    # Retained intron belongs to long mRNA
    ri_exon_name = "%s.withRI" %(long_mRNA_name)
    # Retained intron record has same start/end as gene record as well
    ri_exon_rec = gffutils.Feature(seqid=chrom,
                                   source=source,
                                   featuretype="exon",
                                   start=gene_start,
                                   end=gene_end,
                                   strand=strand,
                                   attributes={"ID": [ri_exon_name],
                                               "Parent": [long_mRNA_name]})
    # Output mRNA splicing out the intron and then output its exons
    short_mRNA_name = "%s.B" %(ri_name)
    # Short mRNA has same start/end as gene record
    short_mRNA_rec = gffutils.Feature(seqid=chrom,
                                      source=source,
                                      featuretype="mRNA",
                                      start=gene_start,
                                      end=gene_end,
                                      strand=strand,
                                      attributes={"ID": [short_mRNA_name],
                                                  "Parent": [ri_name]})
    up_exon_name = "%s.up" %(short_mRNA_name)
    up_exon_rec = gffutils.Feature(seqid=chrom,
                                   source=source,
                                   featuretype="exon",
                                   start=donor.gff_start,
                                   end=donor.gff_end,
                                   strand=strand,
                                   attributes={"ID": [up_exon_name],
                                               "Parent": [short_mRNA_name]})
    dn_exon_name = "%s.dn" %(short_mRNA_name)
    dn_exon_rec = gffutils.Feature(seqid=chrom,
                                   source=source,
                                   featuretype="exon",
                                   start=acceptor.gff_start,
                                   end=acceptor.gff_end,
                                   strand=strand,
                                   attributes={"ID": [dn_exon_name],
                                               "Parent": [short_mRNA_name]})
    # Serialize records to GFF
    # gene
    gff_out.write_rec(gene_rec)
    # long mRNA
    gff_out.write_rec(long_mRNA_rec)
    # retained intron
    gff_out.write_rec(ri_exon_rec)
    # short mRNA
    gff_out.write_rec(short_mRNA_rec)
    gff_out.write_rec(up_exon_rec)
    gff_out.write_rec(dn_exon_rec)



class Unit(namedtuple("Unit", ["start", "end"])):
    """
    Unit is an exon like part of a transcript.
    """
    __slots__ = ()
    @property
    def start_coord(self):
        return int(self.start[1])

    @property
    def end_coord(self):
        return int(self.end[1])

    @property
    def gff_start(self):
        return min((int(self.start[1]), int(self.end[1])))

    @property
    def gff_end(self):
        return max((int(self.start[1]), int(self.end[1])))        

    @property
    def chrom(self):
        return self.start[0]

    @property
    def strand(self):
        return self.start[-1]

    @property 
    def coords_str(self):
        return "%s:%s-%s:%s" %(self.chrom,
                               str(self.start_coord),
                               str(self.end_coord),
                               self.strand)

    @property
    def gff_coords_str(self):
        """
        Return coordinates string in a way that start < end
        always (for GFF purposes.)
        """
        s = self.start_coord
        e = self.end_coord
        if s > e:
            s, e = e, s
        return "%s:%s-%s:%s" %(self.chrom,
                               s,
                               e,
                               self.strand)
    
    @property
    def minus_coords_str(self):
        """
        Version of coords_str() that returns minus strand
        coordinates convention, where start > end (contrary
        to GFF) so reflect transcript order.
        """
        if self.strand == "+":
            raise Exception, "Why call minus_coords_str on a plus strand unit?"
        return "%s:%s-%s:%s" %(self.chrom,
                               str(self.end_coord),
                               str(self.start_coord),
                               self.strand)
        
    def __str__(self):
        return "Unit(%s)" %(self.coords_str)

    
class SpliceGraph:
    """
    Represent the possible splicing graph transitions.
    """
    def __init__(self, table_fnames):
        self.table_fnames = table_fnames
        # Mapping from table name to table
        self.tables = {}
        self.acceptors_to_donors = {"+": Acceptors(strand="+"),
                                    "-": Acceptors(strand="-")}
        self.donors_to_acceptors = {"+": Donors(strand="+"),
                                    "-": Donors(strand="-")}
        self.all_nodes = OrderedDict()
        # Load the UCSC tables
        self.load_tables()
        # Populate the splice graph
        self.populate_graph()
        

    def add_edge(self, donor_unit, acceptor_unit, strand):
        """
        Add edge.
        """
        # Record donor -> acceptor edge
        self.donors_to_acceptors[strand].add_edge(donor_unit, acceptor_unit)
        # Record acceptor <- donor edge
        self.acceptors_to_donors[strand].add_edge(acceptor_unit, donor_unit)
        # Record nodes non-redundantly
        if donor_unit not in self.all_nodes:
            self.all_nodes[donor_unit] = True
        if acceptor_unit not in self.all_nodes:
            self.all_nodes[acceptor_unit] = True


    def has_node(self, node):
        """
        Return True if the node (i.e. exon) exists in the splice graph,
        meaning it occurs in some transcript.
        """
        return (node in self.all_nodes)


    def count_donor_to_acceptor(donor, acceptor):
        """
        Return frequency of donor to acceptor, i.e. how often
        they are spliced to each other.
        """
        pass
        

    def load_tables(self):
        """
        Load tables.
        """
        print "Loading tables..."
        for table_fname in self.table_fnames:
            table_label = os.path.basename(table_fname)
            self.tables[table_label] = parseTables.readTable(table_fname)
                    

    def populate_graph(self):
        """
        Add edges from acceptors to donors, donors to acceptors,
        on distinct strands.
        """
        print "Populating graph..."
        t1 = time.time()
        for table_name in self.tables:
            print "Adding splice edges from table %s" %(table_name)
            for item in self.tables[table_name]:
                chrom, startvals, endvals, strand, gene = item
                startvals = map(int, startvals.split(",")[:-1])
                # Adds +1 since downloaded UCSC tables are 0-based start!
                startvals = map(str, [x + 1 for x in startvals])
                endvals = endvals.split(",")[:-1]
                indices = range(len(startvals))
                if strand == "-":
                    # If it's a minus strand event, walk the transcript from
                    # end (in order of transcription)
                    startvals = startvals[::-1]
                    endvals = endvals[::-1]
                for curr_i, next_i in utils.iter_by_pair(indices, step=1):
                    # Splice from end of current exon to start of next exonp
                    donor_unit = Unit((chrom, startvals[curr_i], strand),
                                      (chrom, endvals[curr_i], strand))
                    acceptor_unit = Unit((chrom, startvals[next_i], strand),
                                         (chrom, endvals[next_i], strand))
                    if strand == "-":
                        # Reverse start/end of donor and acceptor units
                        # if it's a minus strand event
                        donor_unit = Unit((chrom, endvals[curr_i], strand),
                                          (chrom, startvals[curr_i], strand))
                        acceptor_unit = Unit((chrom, endvals[next_i], strand),
                                             (chrom, startvals[next_i], strand))
#                        donor_unit, acceptor_unit = acceptor_unit, donor_unit
                    # Record splice site as edge
                    self.add_edge(donor_unit, acceptor_unit,
                                  strand=strand)
        t2 = time.time()
        print "Populating graph took %.2f seconds" %(t2 - t1)


    def __str__(self):
        F = "\n".join(["%s->%s" %(donor, self.donors["+"].edges[donor]) \
                      for donor in self.donors["+"].edges])
        R = "\n".join(["%s->%s" %(donor, self.donors["-"].edges[donor]) \
                      for donor in self.donors["-"].edges])
        s = "SpliceGraph(\nF=%s\nR=%s\n)" %(F,R)
        return s
 

    def output_gff(self):
        pass



def main():
    table_fname = \
        os.path.expanduser("/home/yarden/jaen/rnaseqlib/rnaseqlib/test/test-data/ri-test/ensGene.txt")
    tables = [table_fname]
    sg = SpliceGraph(tables)
    gff_out = gffutils.gffwriter.GFFWriter("./ri_test.gff3")
    define_RI(sg, gff_out)
    gff_out.close()

    
        

if __name__ == "__main__":
    main()


