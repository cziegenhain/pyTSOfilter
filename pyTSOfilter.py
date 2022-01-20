import os
import pysam
import argparse
import multiprocessing as mp
import regex
#from Bio import SeqIO

#make the lookup dictionary for reverse complementing
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

def reverse_complement(seq):
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def complementseq(seq):
    bases = list(seq)
    bases = [complement.get(base,base) for base in bases]
    bases = ''.join(bases)
    return bases

def idx_bam(bam, threads):
    threads = str(threads)
    try:
        pysam.index("-@"+threads,bam)
    except:
        outcome = 'idxerror'
    else:
        outcome = 'idxsuccess'
    if outcome == 'idxerror':
        print("indexing failed, trying to sort bam file...")
        inbam = bam
        bam = bam+".sorted.bam"
        pysam.sort("-@"+threads,"-o", bam, inbam)
        print("indexing bam file...")
        pysam.index(bam)
    return bam

def makeBAMheader(args, v):
    bam = pysam.AlignmentFile(args.bam, 'rb')
    hdr = bam.header.to_dict()
    bam.close()
    cmdlinecall = 'pyTSOfilter --bam '+args.bam+' --out '+args.out+' --fa '+args.fa+' --p '+str(args.p)+' --n_mismatch '+str(args.n_mismatch)
    if args.ggg:
        cmdlinecall = cmdlinecall+' --ggg'
    pg = {'ID': 'pyTSOfilter', 'PN': 'pyTSOfilter',
          'CL': cmdlinecall, 'VN': v}
    if 'PG' in hdr:
        pglines = hdr['PG']
        pglines.append(pg)
    else:
        pglines = [pg]
    hdr['PG'] = pglines
    return hdr

def is_fiveprime_read(read):
    fivep = False
    if read.has_tag('UB'):
        UMIseq = read.get_tag('UB')
        if UMIseq != "":
            fivep = True
    return fivep

def get_gene(read):
    if read.has_tag("GE"):
        gene = read.get_tag("GE")
    else:
        if read.has_tag("GI"):
            gene = read.get_tag("GI")
        else:
            gene = "none"
    return gene

def collect_bam_chunks(inpath, chrs, outpath):
    allpaths = [inpath+".tmp."+c+".bam" for c in chrs[:-1]]
    allpaths.append(inpath+".tmp."+"unmapped"+".bam")
    cat_args = ['-o', outpath]+allpaths
    pysam.cat(*cat_args)
    x = [os.remove(f) for f in allpaths]

def filterTSO(inpath, threads, fastapath, chr, mm, ggg, outheader):
    fa = pysam.FastaFile(fastapath)
    if chr == '*':
        chrlabel = 'unmapped'
    else:
        chrlabel = chr
        #reflen = fa.get_reference_length(chr)
    #open in/out files
    outpath = inpath+".tmp."+chrlabel+".bam"
    inp = pysam.AlignmentFile(inpath, 'rb', threads = threads)
    out = pysam.AlignmentFile(outpath, 'wb', header = outheader, threads = threads)
    for read in inp.fetch(chr):
        isfivep = is_fiveprime_read(read)
        gene = get_gene(read)
        if chr in fa.references and isfivep and gene != "none":
            if read.is_reverse:
                strand = "-"
                rloffset = read.query_length
            else:
                strand = "+"
            if read.is_paired:
                readtype = "PE"
                if read.is_read1:
                    startpos = read.reference_start
                else:
                    if read.mate_is_unmapped:
                        continue
                    startpos = read.next_reference_start
                    if read.mate_is_reverse:
                        strand = "-"
                        rloffset = read.query_length-22
                    else:
                        strand = "+"
            else:
                startpos = read.reference_start
                readtype = "SE"
            UMIseq = read.get_tag('UB')
            if ggg:
                UMIseq = UMIseq+"GGG"
            if strand == "-":
                #if(startpos+rloffset+21 > reflen):
                #    upstream_seq = fa.fetch(chr, startpos+rloffset, reflen)
                #else:
                #    upstream_seq = fa.fetch(chr, startpos+rloffset, startpos+rloffset+21)
                #if the fetch is longer than reflen, it just returns empty string but no error
                upstream_seq = fa.fetch(chr, startpos+rloffset, startpos+rloffset+21)
                upstream_seq = reverse_complement(upstream_seq)
            else:
                if startpos-20 < 0:
                    startpos = startpos + 20-startpos
                upstream_seq = fa.fetch(chr, startpos-20, startpos+1)
            pattern = r"(" + UMIseq+ "){s<="+ str(mm) + "}"
            m = regex.search(pattern, upstream_seq)
            if m is not None:
                continue
        out.write(read)
    inp.close()
    out.close()



def main():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('--bam', type=str, metavar='FILENAME',
                        help='Path to input BAM file', required = True)
    parser.add_argument('--out', type=str, metavar='FILENAME',
                        help='Path to output bam file', required = True)
    parser.add_argument('--fa', type=str, metavar='FILENAME',
                        help='Path to genome reference fasta', required = True)
    parser.add_argument('--p', type=int, default = 10,
                        help='Number of processes to use')
    parser.add_argument('--n_mismatch', type=int, default = 1,
                        help='Number of mismatches allowed when matching TSO')
    parser.add_argument('--ggg', action = 'store_true',
                        help='Require template-switching Gs as part of sequence match')



    args = parser.parse_args()
    v = '0.1.2'
    print("pyTSOfilter.py v"+v)
    print("WARNING: This script only works correctly when using stranded gene assignment via zUMIs.")
    bampath = args.bam

    print("Loading reference sequence...")
    #reffa = '/home/chris/projects/pyTSOfilter/hg38.primary_assembly.fa'
    reffa = args.fa
    #fa_dict = SeqIO.to_dict(SeqIO.parse(reffa, "fasta"))

    #check if fasta file is indexed
    try:
        fa = pysam.FastaFile(reffa)
    except ValueError:
        print("Error: Reference fasta file is not indexed!")
        print("Please run: samtools faidx "+args.fa)
        quit()
    fa.close()

    if not os.path.exists(bampath+'.bai'):
        print("input bam index not found, indexing...")
        bampath = idx_bam(bampath,args.p)

    chrs = pysam.idxstats(bampath).split('\n')
    chrs = [c.split('\t')[0] for c in chrs[:-1]]
    #chrs.append('*')

    if args.p > 20:
        pysam_workers = 2
        n_jobs = int(args.p/2)
    else:
        pysam_workers = 1
        n_jobs = args.p

    #Construct the new bam header to work with
    bamheader = makeBAMheader(args, v)
    print("Working...")

    pool = mp.Pool(n_jobs)
    results = [pool.apply_async(filterTSO, (args.bam,pysam_workers,args.fa,chr,args.n_mismatch,args.ggg,bamheader,  )) for chr in chrs]
    x = [r.get() for r in results]

    print("Creating final output .bam file...")
    collect_bam_chunks(inpath = bampath, chrs = chrs, outpath = args.out)
    print("Indexing final output .bam file...")
    y = idx_bam(args.out,args.p)

    print("Done!")

if __name__ == "__main__":
    main()
