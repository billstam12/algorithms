import gzip
import pysam
import gzip
import pysam


def read_fasta(data):
    reads = 0
    bases = 0
    with open(data, 'rb') as read:
        for id in read:
            seq = next(read)
            reads += 1
            bases += len(seq.strip())
            next(read)
            next(read)

    print("reads", reads)
    print("bases", bases)

def set_exon_length(length, exon_length):
    if(length <= 10):
        exon_length["1-10"] += 1
    elif(length <= 20):
        exon_length["11-20"] += 1
    elif(length <= 50):
        exon_length["21-50"] += 1
    elif(length <= 100):
        exon_length["51-100"] += 1
    elif(length <= 200):
        exon_length["101-200"] += 1
    elif(length <= 500):
        exon_length["201-500"] += 1
    elif(length <= 1000):
        exon_length["501-1000"] += 1
    elif(length <= 5000):
        exon_length["1001-5000"] += 1
    else:
        exon_length["5001-25000"] += 1

        
#FSM: Full Splice Match
#ISM: Incomplete Splice Match
#NIC: Novel in catalog
#NNC: Novel not in catalog
#NO_SPLICE
#unaligned
def transform_categories(xc):
    print(xc)
    
    return {"FSM": xc["FSM"], 
            "ISM": xc["ISM/NIC_known"] if "ISM/NIC_known" in xc.keys() else 0, 
            "NIC": xc["NIC_novel"] if "NIC_novel" in xc.keys() else 0,
            "NNC": xc["Insufficient_junction_coverage_unclassified"] if "Insufficient_junction_coverage_unclassified" in xc.keys() else 0,
            "NO_SPLICE": xc["NO_SPLICE"] if "NO_SPLICE" in xc.keys() else 0,
            "unaligned": xc["unaligned"] if "unaligned" in xc.keys() else 0,
            #"unindexed": xc["uLTRA_unindexed"] if "uLTRA_unindexed" in xc.keys() else 0
    }
                        
def get_read_stats(data):
    samfile = pysam.AlignmentFile(data, "rb")
    xc = {}
    correct_exon_length = {"1-10" : 0, "11-20" : 0, "21-50" : 0, "51-100" : 0, "101-200" : 0,
                   "201-500" : 0, "201-500" : 0, "501-1000" : 0, "1001-5000" : 0, "5001-25000" : 0, }
    
    for read in samfile.fetch():
        length = read.query_length
        if(length != None):
            set_exon_length(length, correct_exon_length)
        if(read.has_tag("XC")):
            a = read.get_tag("XC")
            if a in xc.keys():
                xc[a] += 1
            else:
                xc[a] = 1       
    samfile.close()

    return correct_exon_length, transform_categories(xc)