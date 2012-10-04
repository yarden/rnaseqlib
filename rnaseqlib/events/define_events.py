##
## Main script for defining AS events from a series of databases
##

#import misopy
#import misopy.gff_utils as gff_utils
#import misopy.Gene as gene_utils

##
## Human:
## Refseq: 
##    ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz 
## UCSC AltEvents: 
##    ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownAlt.txt.gz
##    ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz
## Ensembl events: ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ensGene.txt.gz 
## Swiss Institute of Bioinformatics events: 
##    ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/sibGene.txt.gz 
##    ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/sibTxGraph.txt.gz 
##
## Mouse:
## Refseq: 
##    ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/refGene.txt.gz 
## UCSC AltEvents: ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/knownAlt.txt.gz
## Ensembl events: ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/ensGene.txt.gz 
## Swiss Institute of Bioinformatics events: 
##    ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/sibGene.txt.gz 
##    ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/sibTxGraph.txt.gz 
##


def define_events(event_type):
    """
    Define alternative splicing events.

    Output a GFF3 file with event representations and
    related gene annotations/domain in attributes field.
    """
    pass


# Build a transcript graph using refGene and ensGene. 
# Columns for refGene and ensGene are: bin, name, chrom, strand, 
# txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, 
# score, name2, cdsStartStat, cdsEndStat, exonFrames
# 
def defineAltFromTable(table_f):

    # First, iterate through table and associate exons with transcripts
    # and transcripts with genes.
    geneToTx = {}
    txToExons = {}
    geneToExons = {}
    for line in open(table_f):
        vals = line.strip().split("\t")
        bin_, name_, chrom_, strand_, txStart_, txEnd_, cdsStart_,\
            cdsEnd_, exonCount_, exonStarts_, exonEnds_, score_,\
            name2_, cdsStartStat_, cdsEndStat_, exonFrames = vals

        # May need to check if name2 is not defined
        try:
            geneToTx[name2_].append(name_)
        except:
            geneToTx[name2_] = [name_]
        exonStarts = map(int, exonStarts_.split(",")[:-1])
        exonStarts = map(str, [s + 1 for s in exonStarts])
        exonEnds = exonEnds_.split(",")[:-1]
        for i in range(len(exonStarts)):
            exon = ":".join([chrom_, exonStarts[i], exonEnds[i], strand_])
            try:
                txToExons[name_].append(exon)
            except:
                txToExons[name_] = [exon]
            if name2_ not in geneToExons:
                geneToExons[name2_] = {}
            try:
                geneToExons[name2_][exon] += 1
            except:
                geneToExons[name2_][exon] = 1

    print len(geneToTx), 'genes'
    print len(txToExons), 'transcripts'

    # Now iterate through each gene and identify exons that are not in all
    # transcripts for the gene.
    nAltExons = 0
    for gene in geneToTx:
        txs = geneToTx[gene] 
        exons = geneToExons[gene]
        altexons = [exon for exon in geneToExons[gene] \
            if geneToExons[gene][exon] < len(txs)]
        nAltExons += len(altexons) 
   
    print nAltExons 
    # Get flanking exons for each event and define event type.
    # TODO 
    

# Get constitutive regions which flank alternative events for the
# knownAlt table.  Find the transcripts that contain the 
# altevents exons, and then get constitutive exons.
def defineConstitutiveForAltEvents(altevents_f, kg_f):

    exonToTx = {}
    txToExons = {}
    for line in open(kg_f):
        vals = line.strip().split("\t")
        name_, chrom_, strand_, txStart_, txEnd_, cdsStart_, cdsEnd_,\
            exonCount_, exonStarts_, exonEnds_, proteinID_, alignID_ = vals
        exonStarts = map(int, exonStarts_.split(",")[:-1])
        exonStarts = map(str, [s + 1 for s in exonStarts])
        exonEnds = exonEnds_.split(",")[:-1]
        exons = []
        for i in range(len(exonStarts)):
            exon = ":".join([chrom_, exonStarts[i], exonEnds[i], strand_])
            try:
                exonToTx[exon].append(name_)
            except:
                exonToTx[exon] = [name_]
            exons.append(exon)
        txToExons[name_] = exons

    nAltEvents = 0
    missing = 0
    for line in open(altevents_f):
        vals = line.strip().split("\t")
        bin_, chrom_, start_, end_, type_, score_, strand_ = vals
        start = str(int(start_) + 1)
        exon = ":".join([chrom_, start, end_, strand_])
        if exon in exonToTx:
            txs = exonToTx[exon]
            nAltEvents += 1
        else:
            print "Missing transcripts for", exon        
            missing += 1
       
    print nAltEvents, 'alternative events found'
    print missing, 'missing'

# Parse sibAlt events.



def define_se():
    """
    Define skipped exons
    """
    pass

def define_multi_se():
    """
    Define multiply skipped exons.
    """
    pass

def define_a3ss():
    """
    A3SS
    """
    pass

def define_a5ss():
    """
    A5SS
    """
    pass

def define_ri():
    """
    Retained introns.
    """
    pass

def define_afe():
    """
    AFE
    """
    pass

def define_ale():
    """
    ALE
    """
    pass

def define_mxe():
    """
    MXE
    """
    pass

def define_tandem_utr():
    """
    Tandem UTRs.
    """
    pass
    

def main():
    pass

if __name__ == "__main__":
    main()
