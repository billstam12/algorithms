import pickle
import re
from collections import defaultdict
import math
import os
import pysam
import gffutils
import intervaltree
import json

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].split()[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break

                
def decide_primary_locations(sam_file): # maybe this function is not needed if only one primary alignment from minimap2
    SAM_file = pysam.AlignmentFile(sam_file, "r", check_sq=False)
    reads_primary = {}
    reads_multiple_primary = set()
    # reads_tmp = {}

    for read in SAM_file.fetch(until_eof=True):
        if read.flag == 0 or read.flag == 16:
            # ins = sum([length for type_, length in read.cigartuples if type_ == 1])
            # del_ = sum([length for type_, length in read.cigartuples if type_ == 2 and length < min_intron ])
            # subs = sum([length for type_, length in read.cigartuples if type_ == 8])
            # matches = sum([length for type_, length in read.cigartuples if type_ == 7])
            # tot_align = ins + del_ + subs + matches
            # identity = matches/float(tot_align)

            if read.query_name in reads_primary:
                reads_multiple_primary.add(read.query_name)

                # if identity >= reads_tmp[read.query_name][0] and  matches >= reads_tmp[read.query_name][1]:
                #     reads_primary[read.query_name] = read
                #     reads_tmp[read.query_name] = (identity, matches)
                # elif identity <= reads_tmp[read.query_name][0] and  matches <= reads_tmp[read.query_name][1]:
                #     continue
                # else:
                #     if identity * matches > reads_tmp[read.query_name][0] * reads_tmp[read.query_name][1]:
                #         reads_primary[read.query_name] = read
                #         reads_tmp[read.query_name] = (identity, matches)
                #     else: 
                #         continue


            else:
                reads_primary[read.query_name] = read
                # reads_tmp[read.query_name] = (identity, matches)
    print("TOTAL READS FLAGGED WITH MULTIPLE PRIMARY:", len(reads_multiple_primary))
    return reads_primary    

def pickle_dump(data, filename):
    with open(filename, "wb") as f:
        # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)


def pickle_load(filename):
    with open(filename, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        data = pickle.load(f)
    return data

def modify_ref_header_for_alignment(header):
    if header.isdigit() or header == 'X' or header == 'Y':
        return 'chr'+ header
    elif header == 'MT':
        return 'chrM'
    else:
        return header
    

def get_exon_sites(cigar_tuples, first_exon_start, last_exon_end, annotated_chr_coordinate_pairs):
    ref_pos = first_exon_start
    exon_sites = [ref_pos]
    
    for i, (l,t) in enumerate(cigar_tuples):
        if t == "D":
                if (ref_pos, ref_pos + l) in annotated_chr_coordinate_pairs:
                    exon_sites.append( ref_pos )
                    exon_sites.append( ref_pos + l )
                    # print("HEERE")
                ref_pos += l

        elif t == "=" or t== "M" or t == "X":
            ref_pos += l
        elif t == "N":
            splice_start = ref_pos
            ref_pos += l
            splice_stop = ref_pos
            exon_sites.append( splice_start)
            exon_sites.append( splice_stop )

        elif t == "I" or t == "S" or t == "H": # insertion or softclip
            ref_pos += 0

        else: # reference skip or soft/hardclip "~", or match =
            print("UNEXPECTED!", t)
            sys.exit()

    exon_sites.append(last_exon_end)
    exon_sites = [(exon_sites[i],exon_sites[i+1]) for i in range(0, len(exon_sites), 2) ]  
    return exon_sites

def get_read_alignment_exon_sites(reads_primary_locations, annotated_splice_coordinates_pairs):
    read_exon_sites = {}
    for acc in reads_primary_locations:
        read = reads_primary_locations[acc]
        assert read.flag == 0 or read.flag == 16
        q_cigar = read.cigarstring
        q_start = read.reference_start
        q_end = read.reference_end
        read_cigar_tuples = []
        result = re.split(r'[=DXSMINH]+', q_cigar)
        i = 0
        for length in result[:-1]:
            i += len(length)
            type_ = q_cigar[i]
            i += 1
            read_cigar_tuples.append((int(length), type_ ))  

        mod_ref = modify_ref_header_for_alignment(read.reference_name) # convert 1,..,22,X,Y,MT to chr1, .., chr22, chrX, chrY, chrM
        if mod_ref in annotated_splice_coordinates_pairs:
            annotated_chr_coordinate_pairs = annotated_splice_coordinates_pairs[mod_ref]
        else:
            annotated_chr_coordinate_pairs = set()

        read_exon_sites[read.query_name] = {}  
        read_exon_sites[read.query_name][mod_ref] = get_exon_sites(read_cigar_tuples, q_start, q_end, annotated_chr_coordinate_pairs)

    return read_exon_sites


def get_alignment_classifications(true_exon_sites, aligned_exon_sites, is_ultra = False):
    read_annotations = {}
    total_count_exon_sizes = defaultdict(int)
    correct_count_exon_sizes = defaultdict(int)

    for acc in true_exon_sites:
        exons_true = true_exon_sites[acc][list(true_exon_sites[acc].keys())[0]]
        # if is_ultra:
        #     acc_mod = mod_long_acc(acc)
        # else:
        #     acc_mod = acc
        if acc not in aligned_exon_sites:
            read_annotations[acc] = ('unaligned', len(exons_true), 0, 0)
        else:
            chr_, exons_aligned = list(aligned_exon_sites[acc].items())[0]
            # exons_aligned = list(aln_dict.values())[0]
            aln_start, aln_stop = exons_aligned[0][0], exons_aligned[-1][1]
            true_start, true_stop = exons_true[0][0], exons_true[-1][1]
            if aln_start > true_stop or  aln_stop < true_start:
                read_annotations[acc] = ("diff_location", len(exons_true), len(exons_aligned), 0)

            elif len(exons_true) != len(exons_aligned):
                read_annotations[acc] = ("diff_exon_count", len(exons_true), len(exons_aligned), 0)
            else:
                max_diff = 0
                for i in range(len(exons_true)):
                    (true_start, true_stop) = exons_true[i]
                    (aln_start, aln_stop) = exons_aligned[i]
                    tmp_diff = max(math.fabs(true_start - aln_start), math.fabs(true_stop - aln_stop) )
                    
                    if tmp_diff <= 15:
                        correct_count_exon_sizes[true_stop-true_start + 1] += 1

                    if tmp_diff > max_diff:
                        max_diff = tmp_diff

                if max_diff > 15:
                    is_correct = False
                else:
                    is_correct = True

                if is_correct:
                    read_annotations[acc] = ("correct", len(exons_true), len(exons_aligned), max_diff)
                else:
                    read_annotations[acc] = ("site_diff", len(exons_true), len(exons_aligned), max_diff)

        for (true_start, true_stop) in exons_true:
            total_count_exon_sizes[true_stop-true_start + 1] += 1

    return read_annotations, total_count_exon_sizes, correct_count_exon_sizes

def print_correctness_per_exon_size(correctness_per_exon_size_outfile, total_count_exon_sizes, correct_count_exon_sizes, alignment_method):
    for e_size in sorted(list(total_count_exon_sizes.keys())):
        nr_total = total_count_exon_sizes[e_size]
        if e_size in correct_count_exon_sizes:
            nr_corr = correct_count_exon_sizes[e_size]
        else:
            nr_corr = 0
        correctness_per_exon_size_outfile.write("{0},{1},{2},{3},{4}\n".format(e_size, nr_total, nr_corr, float(nr_corr)/nr_total, alignment_method))
        

def get_read_candidate_splice_sites(reads_primary_locations, annotated_splice_coordinates_pairs):
    read_splice_sites = {}
    for acc in reads_primary_locations:
        read = reads_primary_locations[acc]
        if read.flag == 0 or read.flag == 16:
            # compare cs tag at intron sites
            q_cigar = read.cigarstring
            q_start = read.reference_start
            q_end = read.reference_end
            read_cigar_tuples = []
            result = re.split(r'[=DXSMINH]+', q_cigar)
            i = 0
            for length in result[:-1]:
                i += len(length)
                type_ = q_cigar[i]
                i += 1
                read_cigar_tuples.append((int(length), type_ ))  

            mod_ref = read.reference_name #modify_ref_header_for_alignment(read.reference_name) # convert 1,..,22,X,Y,MT to chr1, .., chr22, chrX, chrY, chrM
            if mod_ref in annotated_splice_coordinates_pairs:
                annotated_chr_coordinate_pairs = annotated_splice_coordinates_pairs[mod_ref]
            else:
                annotated_chr_coordinate_pairs = set()

            read_splice_sites[read.query_name] = {}  
            read_splice_sites[read.query_name][mod_ref] = get_splice_sites(read_cigar_tuples, q_start, annotated_chr_coordinate_pairs)

    return read_splice_sites


def get_annotated_splicesites(ref_gff_file, infer_genes, outfolder):
    db_name = os.path.join(outfolder, 'database.db')
    if infer_genes:
        fn = gffutils.example_filename(os.path.abspath(ref_gff_file))
        db = gffutils.create_db(fn, dbfn=db_name, force=True, keep_order=True, merge_strategy='merge', 
                                sort_attribute_values=True)
        db = gffutils.FeatureDB(db_name, keep_order=True)
    else:
        fn = gffutils.example_filename(os.path.abspath(ref_gff_file))
        db = gffutils.create_db(fn, dbfn=db_name, force=True, keep_order=True, merge_strategy='merge', 
                                sort_attribute_values=True, disable_infer_genes=True, disable_infer_transcripts=True)
        db = gffutils.FeatureDB(db_name, keep_order=True)


    exon_intervals = defaultdict(intervaltree.IntervalTree)

    splice_coordinates = {} # to calc individual fraction of correct sites and NIC
    splice_coordinates_pairs = {} 
    ref_isoforms = {} # To calculate Full splice matches
    # minimum_annotated_intron = 1000000000
    for gene in db.features_of_type('gene'):
        chromosome = str(gene.seqid)
        if chromosome not in ref_isoforms:
            ref_isoforms[chromosome] = {}
        if chromosome not in splice_coordinates:
            splice_coordinates[chromosome] = set()
            splice_coordinates_pairs[chromosome] = set()

        #add annotated transcripts
        for transcript in db.children(gene, featuretype='transcript', order_by='start'):
            # print(dir(transcript))
            consecutive_exons = [exon for exon in db.children(transcript, featuretype='exon', order_by='start')]
            # print('transcript', transcript.id, transcript.start, transcript.stop, [ (exon.start, exon.stop) for exon in db.children(transcript, featuretype='exon', order_by='start')])
            # for individual sites
            for e in consecutive_exons:
                splice_coordinates[chromosome].add(e.stop)
                splice_coordinates[chromosome].add(e.start -1 )
                exon_intervals[chromosome].addi(e.start -1, e.stop, None)


            # for splice pairs
            tmp_splice_sites = []
            for e1,e2 in zip(consecutive_exons[:-1], consecutive_exons[1:]):
                tmp_splice_sites.append( (e1.stop, e2.start -1 ))           
                splice_coordinates_pairs[chromosome].add( (e1.stop, e2.start -1 ) )
            
            ref_isoforms[chromosome][tuple(tmp_splice_sites)] = transcript.id

    return ref_isoforms, splice_coordinates, splice_coordinates_pairs, exon_intervals


def get_splice_sites(cigar_tuples, first_exon_start, annotated_chr_coordinate_pairs):
    splice_sites = []
    ref_pos = first_exon_start
    
    for i, (l,t) in enumerate(cigar_tuples):
        if t == "D":
            # if l >= minimum_annotated_intron -1:
            #     # print("long del", l)
            #     splice_start = ref_pos
            #     ref_pos += l
            #     splice_stop = ref_pos
            #     splice_sites.append( (splice_start, splice_stop) )
                
            # else:
                if (ref_pos, ref_pos + l) in annotated_chr_coordinate_pairs:
                    splice_sites.append( (ref_pos, ref_pos + l) )
                    #print("HEERE")
                ref_pos += l

                # if l > 15:
                #     print("Large deletion!!!", l)


        elif t == "=" or t== "M" or t == "X":
            ref_pos += l
        elif t == "N":
            splice_start = ref_pos
            ref_pos += l
            splice_stop = ref_pos
            splice_sites.append( (splice_start, splice_stop) )

        elif t == "I" or t == "S" or t == "H": # insertion or softclip
            ref_pos += 0

        else: # reference skip or soft/hardclip "~", or match =
            print("UNEXPECTED!", t)
            sys.exit()

    return splice_sites

def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)


from collections import namedtuple
def get_splice_classifications(annotated_ref_isoforms, annotated_splice_coordinates, annotated_splice_coordinates_pairs,  all_reads_splice_sites, ref_seqs, reads_primary_locations, total_no_of_reads):

    total_reads = 0
    
    #individual splice sites
    total_individual_true = 0
    total_individual_in_data = 0

    #pairs of splice sites
    total_pairs_in_data = 0
    tot_nic_pairs = 0
    tot_sm_pairs = 0

    #whole structure
    total_transcript_fsm = 0
    total_transcript_nic = 0
    total_transcript_ism = 0
    total_transcript_nnc = 0
    total_transcript_no_splices = 0

    read_annotations = {}
    alignments_not_matching_annotated_sites = 0
    # print(annotated_ref_isoforms["SIRV5"])
    # for tr in annotated_ref_isoforms["SIRV5"]:
    #     print(tr)
    canonical_splice = 0
    all_splice = 0
    for read_acc in all_reads_splice_sites:
        read_annotations[read_acc] = {}
        total_reads += 1
        assert len(all_reads_splice_sites[read_acc]) == 1
        chr_id, read_splice_sites = list(all_reads_splice_sites[read_acc].items())[0]
        
        # print(reads_primary_locations[read_acc].flag, read_splice_sites)
        # print(chr_id)
        # print(annotated_splice_coordinates)
        if chr_id not in annotated_splice_coordinates:
            alignments_not_matching_annotated_sites += 1
            annotated_sites = set()
            annotated_pairs = set()
            annotated_isoforms = set()
        else:
            annotated_sites = annotated_splice_coordinates[chr_id]
            annotated_pairs = annotated_splice_coordinates_pairs[chr_id]
            annotated_isoforms = annotated_ref_isoforms[chr_id]
        
        # print(annotated_pairs)
        # print( chr_id, read_splice_sites)
        # print(annotated_ref_isoforms[chr_id])
        read_sm_junctions = 0
        read_nic_junctions = 0
        read_splice_letters = []
        read_splice_choords = []
        total_individual_in_data += len(read_splice_sites)*2
        total_pairs_in_data += len(read_splice_sites)
        if chr_id not in ref_seqs.keys(): continue;
        for splice_site in read_splice_sites:
            start_sp, stop_sp = splice_site
            if reads_primary_locations[read_acc].flag == 0:
                donor = ref_seqs[chr_id][start_sp: start_sp + 2] 
                acceptor = ref_seqs[chr_id][stop_sp - 2: stop_sp]
            else:
                acceptor = reverse_complement(ref_seqs[chr_id][start_sp: start_sp + 2])
                donor = reverse_complement(ref_seqs[chr_id][stop_sp - 2: stop_sp])

            if donor == "GT" and acceptor == "AG":
                canonical_splice += 1
            all_splice += 1

            read_splice_letters.append( donor + str("-") + acceptor )
            read_splice_choords.append( str(start_sp) + str("-") + str(stop_sp) )
            # print(splice_site)
            if (start_sp, stop_sp) in annotated_pairs:
                tot_sm_pairs += 1 
                total_individual_true += 2
                read_sm_junctions += 1

            elif start_sp in annotated_sites and stop_sp in annotated_sites:
                # print((start_sp, stop_sp), annotated_pairs )
                tot_nic_pairs += 1
                total_individual_true += 2
                read_nic_junctions += 1

            elif start_sp in annotated_sites:
                total_individual_true += 1  

            elif stop_sp in annotated_sites:
                total_individual_true += 1  

            # for sp in splice_site:
            #     # print(sp)
            #     # print(annotated_sites)
            #     total += 1
            #     if sp in annotated_sites:
            #         total_true += 1

        # check set intersection between read splice sites and annotated splice sites
        read_nic, read_ism, read_nnc, read_no_splices = 0,0,0,0
        read_type = ''
        transcript_fsm_id = "NA"
        # print(tuple(read_splice_sites))
        # print(annotated_ref_isoforms[chr_id])
        # print(annotated_ref_isoforms)

        if len(read_splice_sites) > 0:

            if tuple(read_splice_sites) in annotated_isoforms:
                total_transcript_fsm += 1
                read_type = 'FSM'
                transcript_fsm_id = annotated_isoforms[ tuple(read_splice_sites) ]

            elif len(read_splice_sites) == read_sm_junctions + read_nic_junctions:
                if read_nic_junctions >= 1:
                    total_transcript_nic += 1
                    read_type = 'NIC'
                    # print('NIC', read_acc)

                else:
                    total_transcript_ism += 1
                    read_type = 'ISM'

                    # print(read_acc)
                    # print(tuple(read_splice_sites))
                    # for an in annotated_ref_isoforms[chr_id]:
                    #     print(an)
                    # print()
                    # print(annotated_ref_isoforms[chr_id])
            else:
                total_transcript_nnc += 1
                read_nnc = 1
                read_type = 'NNC'
        else:
            total_transcript_no_splices += 1                
            read_type = 'NO_SPLICE'

        read_annotation = namedtuple('Annotation', ['tot_splices', 'read_sm_junctions', 'read_nic_junctions', 'annotation', "donor_acceptors", "donor_acceptors_choords", "transcript_fsm_id" ])
        if read_splice_letters:
            donor_acceptors = ":".join([str(item) for item in read_splice_letters])
            donor_acceptors_choords = ":".join([str(item) for item in read_splice_choords])

        else: 
            donor_acceptors = "NA"
            donor_acceptors_choords = "NA"
        read_annotations[read_acc] = read_annotation( len(read_splice_sites), read_sm_junctions, read_nic_junctions, read_type, donor_acceptors, donor_acceptors_choords, transcript_fsm_id )
                # print("FSM!!")
    # print(annotated_ref_isoforms[chr_id])
    # print( tuple(read_splice_sites))
    print("Total splice sizes found in cigar in reads (individual, pairs):", total_individual_in_data, total_pairs_in_data, "total matching annotations (individual):", total_individual_true,
             "total annotated junctions (splice match pairs):", tot_sm_pairs,  "total NIC junctions in (pairs):", tot_nic_pairs, "total reads aligned:", total_reads)
    print("total transcripts FSM:", total_transcript_fsm)
    print("total transcripts NIC:", total_transcript_nic)
    print("total transcripts ISM:", total_transcript_ism)
    print("total transcripts NNC:", total_transcript_nnc)
    print("total transcripts no splice sites:", total_transcript_no_splices)
    print("total splice sites:", all_splice)
    print("GT-AG splice sites:", canonical_splice)

    counts = {"FSM": total_transcript_fsm, "NIC": total_transcript_nic, "ISM": total_transcript_ism,
              "NNC": total_transcript_nnc, "NO_SPLICE":  total_transcript_no_splices, "unaligned": total_no_of_reads - total_reads}
    return counts, read_annotations


def print_detailed_values_to_file(error_rates, annotations_dict, reads, outfile, read_type, read_alignments, exon_intervals):
    for acc in reads:
        if acc in error_rates:
            err_rate = error_rates[acc]
        else:
            err_rate = "-"

        if acc not in read_alignments:
            reference_name = 'unaligned'
            reference_start = '-'  #read.reference_name, read.reference_start, read.reference_end + 1, read.flag,
            reference_end = '-'
            flag = '-'
            read_class = ("-","-","-","unaligned","-","-","-") # namedtuple('Annotation', ['tot_splices', 'read_sm_junctions', 'read_nic_junctions', 'annotation', "donor_acceptors", "donor_acceptors_choords", "transcript_fsm_id" ])
            is_exonic = '-'
        else:
            read = read_alignments[acc]
            reference_name, reference_start, reference_end, flag = read.reference_name, read.reference_start, read.reference_end + 1, read.flag
            read_class = annotations_dict[acc] 
            is_exonic = 1 if exon_intervals[reference_name].overlaps(reference_start, reference_end) else 0
        read_length = len(reads[acc])
        # is_unaligned_in_other_method = 1 if acc in reads_unaligned_in_other_method else 0
        info_tuple = (acc, read_type, err_rate, read_length, *read_class, reference_name, reference_start, reference_end, flag, is_exonic) # 'tot_splices', 'read_sm_junctions', 'read_nic_junctions', 'fsm', 'nic', 'ism', 'nnc', 'no_splices'  )
        outfile.write( ",".join( [str(item) for item in info_tuple] ) + "\n")



        
def evaluate_splice_annotations(args, sam_file, alg):
    counts_path = os.path.join(args.outfolder, "counts.json")
    splice_annotations_path = os.path.join(args.outfolder, "splice_annotations.csv")
    if (os.path.exists(counts_path) and os.path.exists(splice_annotations_path)):
        return json.load(open(counts_path))
    annotated_ref_isoforms, annotated_splice_coordinates, annotated_splice_coordinates_pairs, exon_intervals = get_annotated_splicesites(args.gtf, True, args.outfolder)
    pickle_dump(annotated_ref_isoforms, os.path.join( args.outfolder, 'annotated_ref_isoforms.pickle') )
    pickle_dump(annotated_splice_coordinates, os.path.join( args.outfolder, 'annotated_splice_coordinates.pickle') )
    pickle_dump(annotated_splice_coordinates_pairs, os.path.join( args.outfolder, 'annotated_splice_coordinates_pairs.pickle') )
    pickle_dump(exon_intervals, os.path.join( args.outfolder, 'exon_intervals.pickle') )
    if args.simulated:
        error_rates = get_error_rates(reads)
    else:
        error_rates = {}
        
    refs = { acc.split()[0] : seq for i, (acc, (seq, _)) in enumerate(readfq(open(args.ref, 'r')))}
    reads = { acc.split()[0] : seq for i, (acc, (seq, qual)) in enumerate(readfq(open(args.reads, 'r')))}
    total_reads = len(reads)
    # modify_reference_headers(refs)
    # print("SHORTEST INTRON:", minimum_annotated_intron)
    # minimum_annotated_intron = max(minimum_annotated_intron,  min_intron)

    detailed_results_outfile = open(os.path.join(args.outfolder, "splice_annotations.csv"), "w")
    detailed_results_outfile.write("acc,algorithm,error_rate,read_length,tot_splices,read_sm_junctions,read_nic_junctions,annotation,donor_acceptors,donor_acceptors_choords,transcript_fsm_id,chr_id,reference_start,reference_end,sam_flag,is_exonic\n")

    
    primary_locations = decide_primary_locations(sam_file)
    aligned_splice_sites = get_read_candidate_splice_sites(primary_locations, annotated_splice_coordinates_pairs)
    counts, read_annotations = get_splice_classifications(annotated_ref_isoforms, annotated_splice_coordinates, annotated_splice_coordinates_pairs, aligned_splice_sites, refs, primary_locations, total_reads)
    json.dump(counts, open(counts_path, "w"  ))
    print_detailed_values_to_file(error_rates, read_annotations, reads, detailed_results_outfile, alg, primary_locations, exon_intervals)

    return counts  

