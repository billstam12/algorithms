import os
import subprocess
from sys import stdout, exit

from collections import defaultdict
from collections import namedtuple

mem = namedtuple('Mem', ['x', 'y', 'c', 'd', 'val', 'j', "exon_part_id"])
globals()[mem.__name__] = mem # Global needed for multiprocessing

def find_mems_mummer(mummer_path, outfolder, read_path, refs_path, mummer_out_path, min_mem):
    with open(mummer_out_path, "w") as output_file:
        # print('Running spoa...', end=' ')
        stdout.flush()
        null = open(os.path.join(outfolder, "mummer_errors.1") , "w")
        subprocess.check_call([ mummer_path,   '-maxmatch', '-l' , str(min_mem),  refs_path, read_path], stdout=output_file, stderr=null)
        # print('Done.')
        stdout.flush()
    output_file.close()


def find_mems_slamem(slamem_path, mummer_path, outfolder, read_path, refs_path, out_path, min_mem):
    # time slaMEM -l 14 /Users/kxs624/tmp/ULTRA/human_test/refs_sequences.fa /Users/kxs624/tmp/ULTRA/human_test_new_flanking_strat/reads_tmp.fq -o /Users/kxs624/tmp/ULTRA/human_test/slamem_test.tx
    # with open(out_path, "w") as output_file:
    tmp = out_path.split("seeds_batch_")[1]
    batch_id =  tmp.split('.')[0]
    stdout.flush()
    stderr_file = open(os.path.join(outfolder, "slamem_stderr_{0}.1".format(batch_id)) , "w")
    stdout_file = open(os.path.join(outfolder, "slamem_stdout_{0}.1".format(batch_id)) , "w")

    try: # slaMEM throws error if no MEMs are found in any of the sequences
        subprocess.check_call([ slamem_path, '-l' , str(min_mem),  refs_path, read_path, '-o', out_path ], stdout=stdout_file, stderr=stderr_file, env = os.environ)
        print("Using SLAMEM")
    except:
        find_mems_mummer(mummer_path, outfolder, read_path, refs_path, out_path, min_mem)
        print("Using MUMMER")

    print('Done.')
    stdout.flush()
    # output_file.close()


def find_nams_strobemap(slamem_path, outfolder, read_path, refs_path, out_path, nr_cores, min_mem):
    # /usr/bin/time -l ./StrobeMap -n 2 -k 9 -w 30 -t 3 -s  -o /Users/kxs624/tmp/ULTRA/human_test/refs_sequences.fa /Users/kxs624/tmp/ULTRA/human_test_new_flanking_strat/reads_tmp.fq
    # StrobeMap -n 2 -k 9 -w 30 -t 3 -s?? -o ~/tmp/STROBEMERS/multithreading/ /Users/kxs624/tmp/ULTRA/dros_tmp/refs_sequences.fa /Users/kxs624/tmp/ULTRA/dros_test/reads_16xrep.fa
    # with open(out_path, "w") as output_file:
    stdout.flush()
    stderr_file = open(os.path.join(outfolder, "strobemap_stderr.1") , "w")
    stdout_file = open(os.path.join(outfolder, "strobemap_stdout.1") , "w")
    try: # slaMEM throws error if no MEMs are found in any of the sequences
        subprocess.check_call([ 'StrobeMap', '-n' , "2", '-k' , '10', '-v' , "11", '-w' , "35", '-C' , '500', '-L' , '1000', '-S', '-t', str(nr_cores), '-s', '-o', outfolder, refs_path, read_path ], stdout=stdout_file, stderr=stderr_file)
        print("Using StrobeMap")
        print([ 'StrobeMap', '-n' , "2", '-k' , '10', '-v' , "11", '-w' , "35", '-C' , '500', '-L' , '1000', '-S', '-t', str(nr_cores), '-s', '-o', outfolder, refs_path, read_path ])
    except:
        find_mems_slamem(slamem_path, outfolder, read_path, refs_path, out_path, min_mem)
        print("An unexpected error happend in StrobeMap, check error log at:", stderr_file)
        print("If you beileive this is a bug in StrobeMap, report an issue at: https://github.com/ksahlin/strobemers")
        print("You can always sidestep this issue by providing another seed finder to uLTRA, i.e., remove option --use_NAM_seeds.")
        sys.exit()
    stdout.flush()

# def parse_results(mems_path):
#     # file = open(os.path.join(mems_folder, "mummer_mems.txt"), 'r')
#     mems_db = {}
#     # read_mems_tmp = {}

#     for i, line in enumerate(open(mems_path, 'r')):
#         if line[0] == '>':
#             if i == 0:
#                 read_acc = line.split()[1].strip()  # mems_db[line.split()[1].strip)()] = [] 
#             else:
#                 mems_db[read_acc] = read_mems_tmp 
#                 read_acc = line.split()[1].strip() 
            
#             read_mems_tmp = defaultdict(list)

#         else:
#             vals =  line.split() #11404_11606           1     11405       202
#             exon_part_id = vals[0]
#             chr_id, ref_coord_start, ref_coord_end = exon_part_id.split('^')
#             mem_len = int(vals[3])
#             mem_ref_exon_part_start = int(vals[1])
#             mem_read_start = int(vals[2])
#             # convert to 0-indexed as python, however last coordinate is inclusive of the hit, not as in python end-indexing
#             mem_tuple = mem(int(ref_coord_start) + mem_ref_exon_part_start - 1, int(ref_coord_start) + mem_ref_exon_part_start -1 + mem_len - 1,
#                             mem_read_start-1, mem_read_start-1 + mem_len - 1, 
#                             mem_len, None, exon_part_id)
            
#             read_mems_tmp[chr_id].append( mem_tuple )
#         # print(line)
#     # add last read
#     mems_db[read_acc] = read_mems_tmp 


#     return mems_db


def get_mem_records(mems_path, reads):
    '''
        Reads contains all the relevant reads in the batch to read mems from 
    '''
    relevant = False
    relevant_read_cnt = 0
    for i, line in enumerate(open(mems_path, 'r')):
        if line[0] == '>':
            acc = line[1:].strip()
            if acc not in reads:
                relevant = False
                continue
            else:
                relevant = True
                relevant_read_cnt +=1

            if relevant_read_cnt == 1:
                read_acc = acc  
            else:

                for chr_id in list(read_mems_tmp.keys()):
                    coordinate_sorted_tuples = sorted(read_mems_tmp[chr_id], key = lambda x: x[1])
                    sorted_mems = [ mem(x,y,c,d,val,j,e_id) for j, (x, y, c, d, val, e_id) in enumerate(coordinate_sorted_tuples) ]
                    read_mems_tmp[chr_id] = sorted_mems

                yield read_acc, read_mems_tmp
                read_acc = acc 
            
            read_mems_tmp = defaultdict(list)

        elif relevant:
                vals =  line.split() #11404_11606           1     11405       202
                exon_part_id = vals[0]
                chr_id, ref_coord_start, ref_coord_end = exon_part_id.split('^')
                chr_id = int(chr_id)
                mem_len = int(vals[3])
                # convert to 0-indexed reference as in python
                # however, for MEM length last coordinate is inclusive of the hit in MEM solvers, not as in python end-indexing
                mem_ref_exon_part_start = int(vals[1]) - 1
                mem_read_start = int(vals[2]) - 1
                ref_coord_start = int(ref_coord_start) # has already been 0-indexed when constructing parts
                mem_genome_start = ref_coord_start + mem_ref_exon_part_start
                
                
                # mem_tuple = mem(int(ref_coord_start) - 1 + mem_ref_exon_part_start - 1, int(ref_coord_start) - 1 + mem_ref_exon_part_start -1 + mem_len - 1,
                #                 mem_read_start-1, mem_read_start-1 + mem_len - 1, 
                #                 mem_len, None, exon_part_id)
                # read_mems_tmp[chr_id].append( mem_tuple )
                info_tuple = ( mem_genome_start, mem_genome_start + mem_len - 1,
                                mem_read_start, mem_read_start + mem_len - 1, 
                                mem_len, exon_part_id)
                read_mems_tmp[chr_id].append( info_tuple )


    for chr_id in list(read_mems_tmp.keys()):
        coordinate_sorted_tuples = sorted(read_mems_tmp[chr_id], key = lambda x: x[1])
        sorted_mems = [ mem(x,y,c,d,val,j,e_id) for j, (x, y, c, d, val, e_id) in enumerate(coordinate_sorted_tuples) ]
        read_mems_tmp[chr_id] = sorted_mems
    print("READ {0} RECORDS.".format(relevant_read_cnt))
    yield read_acc, read_mems_tmp


# def find_file_positions(read_pos_list, mummer_out_path, mummer_out_path_rc):
#     print(read_pos_list)



