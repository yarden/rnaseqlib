##
## SpliceGraph class.
##
import os
import sys
import time

import rnaseqlib
import rnaseqlib.utils as utils
import rnaseqlib.events.parseTables as parseTables

from collections import OrderedDict


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


    def get_donors(self, acceptor):
        """
        Return donors that link up to a particular acceptor.
        """
        pass


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


    def get_acceptors(self, donor):
        """
        Return acceptors that link up to a particular donor.
        """
        pass


    def __str__(self):
        return "Donors(%s)" %(self.edges)


    def __repr__(self):
        return self.__str__()
    
    


def define_RI(sg, multi_iso=False):
    """
    Define retained introns.

    [ A ]---[ B ]
    [     C     ]

    [ A ]---[ B ]
    ----[ C ]----

    Look at B's donors (in this case A) and see if
    they have a start coordinate
    """
    if multi_iso:
        raise Exception, "Multiple isoforms not supported."
    for strand in sg.acceptors:
        for acceptor in sg.acceptors[strand].edges:
            # If any of the donor units to this acceptor
            # have the acceptor end as their end coordinate, it's a retained
            # intron event
            donors_of_acceptor = sg.acceptors[strand].get_edges_from(acceptor)
            for donor in donors_of_acceptor:
                # If there's a node that has this acceptor end as its end coordinate
                # and the donor's start as the start coordinate, then it's a
                # retained intron
                intron_as_exon = Unit(donor.start, acceptor.end)
                if sg.has_node(intron_as_exon):
                    print "RETAINED INTRON!"

            

def RI(DtoA_F, AtoD_F, DtoA_R, AtoD_R, gff3_f,
       multi_iso=False):
    """
    Define retained introns.

    Arguments are:
    splice site dictionaries, followed by output file and a
    keyword to denote method of selecting flanking exon
    coordinates:
    
      - shortest
      - longest
      - commonshortest
      - commonlongest 
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


class Unit:
    """
    An exon-like unit. Start and end.
    """
    def __init__(self, start, end):
        self.start = start
        self.end = end
        # Get chrom/strand from start (should be identical
        # to end!)
        self.chrom = start[0]
        self.strand = start[-1]
        if self.chrom != end[0]:
            raise Exception, "Unit has mismatching chromosomes for start/end."
        

    def __repr__(self):
        return self.__str__()
    

    def __str__(self):
        return "Unit(Start=%s, End=%s)" %(self.start, self.end)
    

class SpliceGraph:
    """
    Represent the possible splicing graph transitions.
    """
    def __init__(self, table_fnames):
        self.table_fnames = table_fnames
        # Mapping from table name to table
        self.tables = {}
        self.acceptors = {"+": Acceptors(strand="+"),
                          "-": Acceptors(strand="-")}
        self.donors = {"+": Donors(strand="+"),
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
        self.donors[strand].add_edge(donor_unit, acceptor_unit)
        # Record acceptor <- donor edge
        self.acceptors[strand].add_edge(acceptor_unit, donor_unit)
        # Record nodes non-redundantly
        if donor_unit not in self.all_nodes:
            self.all_nodes[donor_unit] = True
        if acceptor_unit not in self.all_nodes:
            self.all_nodes[acceptor_unit] = True


    def has_node(self, node):
        if node in self.all_nodes:
            return True
        else:
            return False


    def get_unit_with_start(start_unit, strand):
        """
        Get all units starting with the given start unit
        on the requested strand.
        """
        

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
        first_table_name = self.tables.keys()[0]
        for item in self.tables[first_table_name]: 
            chromval, startvals, endvals, strandval, gene = item
            startvals = map(int, startvals.split(",")[:-1])
            # Adds +1 since downloaded UCSC tables are 0-based start!
            startvals = map(str, [x + 1 for x in startvals])
            endvals = endvals.split(",")[:-1]
            if strandval == '+':
                ##
                ## + strand
                ##
                for i in range(len(startvals) - 1):
                    prevacceptor = ":".join([chromval, startvals[i], strandval])
                    donor = ":".join([chromval, endvals[i], strandval])
                    acceptor = ":".join([chromval, startvals[i + 1], strandval])
                    nextdonor = ":".join([chromval, endvals[i + 1], strandval])

                    # if donor not in ss5_ss3_F:
                    #     ss5_ss3_F[donor] = [] 
                    # ss5_ss3_F[donor].append(acceptor)
                    # if acceptor not in ss3_ss5_F:
                    #     ss3_ss5_F[acceptor] = []
                    # ss3_ss5_F[acceptor].append(nextdonor)

                    # if donor not in ss5_ss3_R:
                    #     ss5_ss3_R[donor] = []
                    # ss5_ss3_R[donor].append(prevacceptor)
                    # if acceptor not in ss3_ss5_R:
                    #     ss3_ss5_R[acceptor] = []
                    # ss3_ss5_R[acceptor].append(donor)
            else:
                ##
                ## - strand
                ##
                # Walk the start and end coordinates with end first
                startvals = startvals[::-1]
                endvals = endvals[::-1]
            for i in range(len(startvals) - 1):
                # Splice from end of current exon to start of next exon
                donor_unit = Unit((chromval, startvals[i], strandval),
                                  (chromval, endvals[i], strandval))
                acceptor_unit = Unit((chromval, startvals[i + 1], strandval),
                                     (chromval, endvals[i + 1], strandval))
                # Record splice site as edge
                self.add_edge(donor_unit, acceptor_unit,
                              strand=strandval)
                    
                    #### OLD CODE
                    # prevacceptor = ":".join([chromval, endvals[i], strandval])
                    # donor = ":".join([chromval, startvals[i], strandval])
                    # acceptor = ":".join([chromval, endvals[i + 1], strandval])
                    # nextdonor = ":".join([chromval, startvals[i + 1], strandval])

#                    if donor not in ss5_ss3_F:
#                        ss5_ss3_F[donor] = [] 
#                    ss5_ss3_F[donor].append(acceptor)
#                    if acceptor not in ss3_ss5_F:
#                        ss3_ss5_F[acceptor] = []
#                    ss3_ss5_F[acceptor].append(nextdonor)

                    
#                    if donor not in ss5_ss3_R:
#                        ss5_ss3_R[donor] = []
#                    ss5_ss3_R[donor].append(prevacceptor)
#                    if acceptor not in ss3_ss5_R:
#                        ss3_ss5_R[acceptor] = []
#                    ss3_ss5_R[acceptor].append(donor)


    def __str__(self):
        print "SELF.DONORS: ", self.donors["+"]
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
    print "SG: "
    print sg.donors["+"]["a"]
    print sg
    define_RI(sg)
    
        

if __name__ == "__main__":
    main()


