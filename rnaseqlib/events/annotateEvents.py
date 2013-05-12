import os, sys, operator, shelve
import parseTables


def getLookup(gff3_f, table_f, xref_f, db_f, out_f):

    ssToGene = parseTables.populateGenelist(table_f)
    print len(ssToGene), 'splice sites'
    geneToInfo = parseTables.readXref(xref_f)

    eventToInfo = shelve.open(db_f, 'c')
    for line in open(gff3_f):
        vals = line.strip().split("\t")
        if vals[2] == 'gene':
            # Get coordinates from gene record
            eventtype = vals[1]
            event = vals[8].split(";")[0].split("ID=")[1]
            sslist = eval(eventtype + '("' + event + '")')
            genes = []
            for ss in sslist:
                if ss in ssToGene:
                    genes.append(ssToGene[ss])
 
            if len(genes) > 0 and genes[0] in geneToInfo:
                symb, desc = geneToInfo[genes[0]]
                eventToInfo[event] = [genes[0], symb, desc]

    out = open(out_f, 'w')
    out.write("#Event\tGene\tSymb\tDesc\n")
    for event in eventToInfo:
        item = [event]
        item.extend(eventToInfo[event])
        out.write("\t".join(item) + "\n")
    out.close()
    print len(eventToInfo), 'events matched to genes'
    eventToInfo.close()


def SE(event):
    exons = event.split("@")
    ss = []
    for exon in exons:
        chrom, start, end, strand = exon.split(":")
        ss1 = ":".join([chrom, start, strand])
        ss2 = ":".join([chrom, end, strand])
        ss.append(ss1)
        ss.append(ss2)
    return ss

def MXE(event):
    return SE(event)

def A5SS(event):
    alt, dn = event.split("@")
    ss = []
    chrom, start, endcoords, strand = alt.split(":")
    ss1 = ":".join([chrom, start, strand])
    ss.append(ss1)
    for end in endcoords.split("|"):
        ss2 = ":".join([chrom, end, strand])
        ss.append(ss2)
    chrom, start, end, strand = dn.split(":")
    ss1 = ":".join([chrom, start, strand])
    ss2 = ":".join([chrom, end, strand])
    ss.append(ss1)
    ss.append(ss2)
    return ss

def A3SS(event):
    up, alt = event.split("@")
    ss = []
    chrom, start, end, strand = up.split(":")
    ss1 = ":".join([chrom, start, strand])
    ss2 = ":".join([chrom, end, strand])
    ss.append(ss1)
    ss.append(ss2)
    chrom, startcoords, end, strand = alt.split(":")
    for start in startcoords.split("|"):
        ss1 = ":".join([chrom, start, strand])
        ss.append(ss1)
    ss2 = ":".join([chrom, end, strand])
    ss.append(ss2)
    return ss

def RI(event):
    alt1, alt2 = event.split("@")
    ss = []
    chrom, start, endcoords, strand = alt1.split(":")
    ss1 = ":".join([chrom, start, strand])
    ss.append(ss1)
    for end in endcoords.split("|"):
        ss2 = ":".join([chrom, end, strand])
        ss.append(ss2)
    chrom, startcoords, end, strand = alt2.split(":")
    ss1 = ":".join([chrom, end, strand])
    ss.append(ss1)
    for end in startcoords.split("|"):
        ss2 = ":".join([chrom, start, strand])
        ss.append(ss2)
    return ss

