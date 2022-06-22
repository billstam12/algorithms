import os 

class arguments:
    def __init__(self, ref, gtf, bed, desalt_gtf, reads, ont, isoseq, disable_infer, name, tool_name, mm2, desalt_index):
        self.index = None
        self.outfolder = os.path.join("output/", name, tool_name)
        self.ont = ont # Set parameters suitable for ONT (Currently sets: --min_mem 17, --min_acc 0.6 --alignment_threshold 0.5).
        self.isoseq = isoseq  # Set parameters suitable for IsoSeq (Currently sets: --min_mem 20, --min_acc 0.8 --alignment_threshold 0.5).
        self.simulated = not (ont or isoseq)
        self.flank_size = 1000 # Size of genomic regions surrounding genes.
        self.use_NAM_seeds = False # Uses StrobeMap to generate NAM seeds instead of MEMs. Uses StrobeMap parameters ./StrobeMap -n 2 -k 15 -v 16 -w 45 -t [nr_cores] -s.
        self.min_mem = 17 # Threshold for minimum mem size considered.
        self.min_segm = 25 # Threshold for minimum segment size considered.
        self.min_acc = 0.5 # Minimum accuracy of MAM to be considered in mam chaining.
        self.max_intron = 1200000 #Set global maximum size between mems considered in chaining solution. This is otherwise inferred from GTF file per chromosome.')
        self.max_loc = 5 # Limit read to be aligned to at most max_loc places.
        self.mask_threshold = 200 #Abundance occurance threshold. Masks more abundant k-mers than this threshold before MEM finding.
        self.alignment_threshold = 0.5 #Lower threshold for considering an alignment. 
        self.small_exon_threshold = 200 # Considered in MAM solution even if they cont contain MEMs.
        self.non_covered_cutoff = 15 #Threshold for what is counted as varation/intron in alignment as opposed to deletion.
        self.dropoff = 0.95 #Ignore alignment to hits with read coverage of this fraction less than the best hit.
        self.nr_cores = 4  # number of cores to be used
        self.mm2_ksize = 15
        self.disable_mm2 = mm2 # Disables use of minimap2
        self.genomic_frac = 0.1 # If parameter prefilter_genomic is set, this is the threshild for fraction of aligned portion of read that is outside uLTRA indexed regions to be considered genomic (default 0.1).
        self.reads_tmp = ""
        self.reduce_read_ployA = 8 # Reduces polyA tails longer than X bases (default 10) in reads to 5bp before MEM finding. This helps MEM matching to spurios regions but does not affect final read alignments.'
        self.reads_rc = "" 
        self.prefix = "reads" # Outfile prefix [default=reads]. "--prefix sample_X" will output a file sample_X.sam.
        self.keep_temporary_files = False # Keeps all intermediate files used for the alignment. This parameter is manily good for bugfixing and development.'

        self.slamem_path = "uLTRA/slaMEM/slaMEM"
        self.mummer_path = "uLTRA/mummer-4.0.0rc1/mummer"
        self.minimap_path = "uLTRA/minimap2/minimap2"
        self.desalt_index = desalt_index
        self.ref = ref # path to reference genome
        self.bed = bed # path to bed annotations file
        self.gtf = gtf # path to annotations file
        self.desalt_gtf = desalt_gtf # path to desalt annotations file
        self.reads = reads # path to reads file
        self.disable_infer = disable_infer
        self.verbose = False
        
