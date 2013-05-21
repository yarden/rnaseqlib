##
## SpliceGraph class.
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.events.parseTables as parseTables

from collections import OrderedDict, namedtuple


import collections

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
              min_intron_len=10,
              multi_iso=False):
    """
    Define retained introns.

    [ A ]---[ B ]
    [     C     ]

    Look at B's donors (in this case A) and see if
    they have a start coordinate.

    Parameters:
    -----------

    sg : SpliceGraph
    min_intron_len : minimum intron length
    """
    if multi_iso:
        raise Exception, "Multiple isoforms not supported."
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
                    intron_len = (acceptor.start_coord - 1) - \
                                 (donor.end_coord + 1) + 1
                    if intron_len < min_intron_len:
                        continue
                    print "It's a retained intron between %s and %s" \
                          %(donor, acceptor)
            

def RI(DtoA_F, AtoD_F, DtoA_R, AtoD_R, gff3_f,
       multi_iso=False):
    """
    Define retained introns.
    """
    print "Generating retained introns (RI)"
#    if os.path.isfile(gff3_f):
#        print "  - Found file, skipping..."
#        return
    out = open(gff3_f, 'w')

    print "ATOD_F: ", AtoD_F

    for acceptor in AtoD_F:                                 # iterate through acceptors
        chrom, acceptorcoord, strand = acceptor.split(":")
        donors = [x for x in list(AtoD_F[acceptor]) if x in DtoA_R]
                                                           # get the 5'ss in same exon for these acceptors

        if len(donors) > 1:
            for donor in donors:                           # iterate through these 5'ss
                rilist = []
                riAcceptors = [x for x in list(DtoA_R[donor]) if x in AtoD_R]
                print "riAcceptors: ", riAcceptors
                                                           # get the upstream 3'ss for this 5'ss
                for riAcceptor in riAcceptors:             # iterate through these 3'ss and get their 5'ss
                    riDonors = [x for x in list(AtoD_R[riAcceptor]) if x in DtoA_R]
                    for riDonor in riDonors:               # if the 3'ss upstream of these 5'ss is the upstream
                        upAcceptors = list(DtoA_R[riDonor])# acceptor, this is a retained intron 
                        if acceptor in upAcceptors:
                            rilist.append([riDonor, riAcceptor])
                if len(rilist) > 0:
                    ss1 = acceptor.split(":")[1]
                    donorlist = list(set([x[0].split(":")[1] for x in rilist]))
                    acceptorlist = list(set([x[1].split(":")[1] for x in rilist]))
                    ss4 = donor.split(":")[1]

                    if multi_iso:
                        if strand == '+':

                            upexon = ":".join([chrom, ss1, "|".join(donorlist), strand])
                            dnexon = ":".join([chrom, "|".join(acceptorlist), ss4, strand])
                            name = "@".join([upexon, dnexon])
                            out.write("\t".join([chrom, 'RI', 'gene', ss1, ss4,\
                                '.', strand, '.', "ID=" + name + ";Name=" + name]) + "\n")
                            out.write("\t".join([chrom, 'RI', 'mRNA', ss1, ss4,\
                                '.', strand, '.', "ID=" + name + ".A;Parent=" + name]) + "\n")
                            for i in range(len(rilist)):
                                out.write("\t".join([chrom, 'RI', 'mRNA', ss1, ss4,\
                                    '.', strand, '.', "ID=" + name + "." + LETTERS[i + 1] + ";Parent=" + name]) + "\n")
                            out.write("\t".join([chrom, 'RI', 'exon', ss1, ss4,\
                                '.', strand, '.', "ID=" + name + ".ri;Parent=" + name + ".A"]) + "\n")
                            for i in range(len(rilist)):
                                out.write("\t".join([chrom, 'RI', 'exon', ss1, donorlist[i],\
                                    '.', strand, '.', "ID=" + name + ".up;Parent=" + name + "." + LETTERS[i + 1]]) + "\n")
                                out.write("\t".join([chrom, 'RI', 'exon', acceptorlist[i], ss4,\
                                    '.', strand, '.', "ID=" + name + ".dn;Parent=" + name + "." + LETTERS[i + 1]]) + "\n")

                        else:

                            upexon = ":".join([chrom, ss1, "|".join(donorlist), strand])
                            dnexon = ":".join([chrom, "|".join(acceptorlist), ss4, strand])
                            name = "@".join([upexon, dnexon])
                            out.write("\t".join([chrom, 'RI', 'gene', ss4, ss1,\
                                '.', strand, '.', "ID=" + name + ";Name=" + name]) + "\n")
                            out.write("\t".join([chrom, 'RI', 'mRNA', ss4, ss1,\
                                '.', strand, '.', "ID=" + name + ".A;Parent=" + name]) + "\n")
                            for i in range(len(rilist)):
                                out.write("\t".join([chrom, 'RI', 'mRNA', ss4, ss1,\
                                    '.', strand, '.', "ID=" + name + "." + LETTERS[i + 1] + ";Parent=" + name]) + "\n")
                            out.write("\t".join([chrom, 'RI', 'exon', ss4, ss1,\
                                '.', strand, '.', "ID=" + name + ".up;Parent=" + name + ".A"]) + "\n")
                            for i in range(len(rilist)):
                                out.write("\t".join([chrom, 'RI', 'exon', donorlist[i], ss1,\
                                    '.', strand, '.', "ID=" + name + ".up;Parent=" + name + "." + LETTERS[i + 1]]) + "\n")
                                out.write("\t".join([chrom, 'RI', 'exon', ss4, acceptorlist[i],\
                                    '.', strand, '.', "ID=" + name + ".dn;Parent=" + name + "." + LETTERS[i + 1]]) + "\n")

                    # 2-isoform case
                    else:
                        for iso1 in range(len(donorlist)):
                            for iso2 in range(len(acceptorlist)):
                                if strand == '+':
                                    upexon = ":".join([chrom, ss1, donorlist[iso1], strand])
                                    dnexon = ":".join([chrom, acceptorlist[iso2], ss4, strand])
                                    name = "@".join([upexon, dnexon])

                                    out.write("\t".join([chrom, 'RI', 'gene', ss1, ss4,\
                                        '.', strand, '.', "ID=" + name + ";Name=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'mRNA', ss1, ss4,\
                                        '.', strand, '.', "ID=" + name + ".A;Parent=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'mRNA', ss1, ss4,\
                                        '.', strand, '.', "ID=" + name + ".B;Parent=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'exon', ss1, ss4,\
                                        '.', strand, '.', "ID=" + name + ".ri;Parent=" + name + ".A"]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'exon', ss1, donorlist[iso1],\
                                        '.', strand, '.', "ID=" + name + ".up;Parent=" + name + ".B"]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'exon', acceptorlist[iso2], ss4,\
                                        '.', strand, '.', "ID=" + name + ".dn;Parent=" + name + ".B"]) + "\n")

                                else:
                                    upexon = ":".join([chrom, ss1, donorlist[iso1], strand])
                                    dnexon = ":".join([chrom, acceptorlist[iso2], ss4, strand])
                                    name = "@".join([upexon, dnexon])
                                    
                                    out.write("\t".join([chrom, 'RI', 'gene', ss4, ss1,\
                                        '.', strand, '.', "ID=" + name + ";Name=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'mRNA', ss4, ss1,\
                                        '.', strand, '.', "ID=" + name + ".A;Parent=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'mRNA', ss4, ss1,\
                                        '.', strand, '.', "ID=" + name + ".B;Parent=" + name]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'exon', ss4, ss1,\
                                        '.', strand, '.', "ID=" + name + ".ri;Parent=" + name + ".A"]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'exon', donorlist[iso1], ss1,\
                                        '.', strand, '.', "ID=" + name + ".up;Parent=" + name + ".B"]) + "\n")
                                    out.write("\t".join([chrom, 'RI', 'exon', ss4, acceptorlist[iso2],\
                                        '.', strand, '.', "ID=" + name + ".dn;Parent=" + name + ".B"]) + "\n")
                                # Record seen RI
                                #seen_RIs[name] = True
    out.close()


class Unit(namedtuple("Unit", ["start", "end"])):
    """
    Unit is an exon like part of a transcript.
    """
    __slots__ = ()
    def __str__(self):
        return "Unit(%s:%s-%s:%s)" %(self.start[0],
                                     self.start[1],
                                     self.end[1],
                                     self.start[2])

    @property
    def start_coord(self):
        return int(self.start[1])

    @property
    def end_coord(self):
        return int(self.end[1])

    @property
    def chrom(self):
        return self.start[0]

    @property
    def strand(self):
        return self.start[-1]
        
    
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
                    startvals = startvals[::-1]
                    endvals = endvals[::-1]
                for curr_i, next_i in utils.iter_by_pair(indices, step=1):
                    # Splice from end of current exon to start of next exonp
                    donor_unit = Unit((chrom, startvals[curr_i], strand),
                                      (chrom, endvals[curr_i], strand))
                    acceptor_unit = Unit((chrom, startvals[next_i], strand),
                                         (chrom, endvals[next_i], strand))
                    if strand == "-":
                        # Reverse donor, acceptor if on minus strand
                        donor_unit, acceptor_unit = acceptor_unit, donor_unit
                    # Record splice site as edge
                    self.add_edge(donor_unit, acceptor_unit,
                                  strand=strand)


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
    define_RI(sg)

    
        

if __name__ == "__main__":
    main()


