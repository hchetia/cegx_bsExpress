#!/usr/bin/env python

import argparse
import sys
import os
import subprocess
import tempfile
import re
import gzip
import shlex
import shutil
import oxbs_qc 
import oxbs_qc_func
import validate_args

VERSION= '0.5cegx' ## MEMO: Change version in setup.py accordingly

CONTROL_FASTA= 'oxBS_controls-v1.0.fa'

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Execute QC of BS/oxBS libraries.

SEE ALSO:
    http://code.google.com/p/oxbs-sequencing-qc/wiki/bsExpressDoc

    """, formatter_class= argparse.RawDescriptionHelpFormatter)

# -----------------------------------------------------------------------------
input_args= parser.add_argument_group('Input options', '')

input_args.add_argument('--input', '-i',
                   required= True,
                   nargs= '+',
                   help='''Input to QC: One or two fastq files or one aligned
bam. If two fastq files are given they are taken as paired-end reads.
                   ''')

input_args.add_argument('--prefix', '-p',
                   required= True,
                   help='''Prefix assigned to all the output files. Default "oxbs_qc."
                   ''')

input_args.add_argument('--outdir', '-o',
                   required= False,
                   default= None,
                   help='''Directory where output will be sent. Created if it doesn't
exist. Default to a directory named as the --prefix.
                   ''')

input_args.add_argument('--ref', '-r',
                   required= True,
                   help='''The path to the reference fasta file . Inside the same
containing the fasta file there should be also: 
- Tab delimited file of base modifications <ref>.txt
- subdir `Bisulfite_Genome` prepared by `bismark_genome_preparation --bowtie2`
                   ''')

input_args.add_argument('--encoding', '-e',
                   required= False,
                   choices= ['Sanger', 'Illumina', None],
                   default= None,
                   help='''Quality encoding of FASTQ files.
Options are 'Sanger' or 'Illumina'. If not set, encoding will be inferred.
                   ''')

input_args.add_argument('--maxlen', '-ml',
                   type= int,
                   default= None,
                   help='''Maximum length of a read before it is shortened.
Ignored if flag --skip_shorten is used. Shortening is skipped if -ml is None and
the reference is not the control sequences (%s).
If -ml is not set and the control sequences are used, ml is set to 58.
                   ''' %(CONTROL_FASTA))

input_args.add_argument('--listpos', '-l',
                   default= None,
                   help='''BED file of positions where methylation should be called.
If none (default), methylation will be called at *each* cytosine in the reference genome.
This argument passed to `mpileup -l`.
                   ''')

input_args.add_argument('--verbose', '-V',
                   action= 'store_true',
                   help='''Produce more verbose progress report. Currently,
--verbose opt only shows the executed commands in addition to standard output.
                   ''')

input_args.add_argument('--version', '-v',
                   action= 'version',
                   version= VERSION,
                   help='''Version
                   ''')

# -----------------------------------------------------------------------------
wflow_opts= parser.add_argument_group('Workflow options', '''
Steps of the pipeline to execute. The order printed here reflect the pipeline
order.''')

wflow_opts.add_argument('--check',
                   action= 'store_true',
                   help='''Perform a check that all the settings and requirements are
satisfied and exit.
                   ''')

wflow_opts.add_argument('--skip_fastqc', '-Sf',
                   action= 'store_true',
                   help='''Do not perform FastQC on input file(s)
                   ''')

wflow_opts.add_argument('--skip_trim', '-St',
                   action= 'store_true',
                   help='''Do not perform trim reads by removing 3'-adapters and
low quality ends. 
                   ''')

wflow_opts.add_argument('--skip_shorten', '-Ss',
                   action= 'store_true',
                   help='''Do not perform shortening of reads. Use this option
if your reads are 3 or more bases shorter than the shortest reference sequence.
                   ''')

wflow_opts.add_argument('--skip_aln', '-Sa',
                   action= 'store_true',
                   help='''Do not align reads. Input is already aligned SAM or BAM.
                   ''')

wflow_opts.add_argument('--skip_sam2bam', '-Sb',
                   action= 'store_true',
                   help='''Do not perform conversion from SAM to BAM, sorting and indexing.
This flag implies that input is a sorted and indexed bam file.
                   ''')

wflow_opts.add_argument('--mark_duplicates', '-Md',
                   action= 'store_true',
                   help='''Mark duplicates in bam file using picard MarkDuplicates.
NB: The mpileup engine skips reads marked as duplicate therefore is not adviceable
to set this option on control DNA where most of the reads are duplicates.
                   ''')

wflow_opts.add_argument('--skip_clip', '-Sc',
                   action= 'store_true',
                   help='''Do not clip paired-end reads. Use this option if your
paired-end reads are not going to overlap after alignment. Clipping is automatically skipped
if the input is a single FASTQ file.
                   ''')

wflow_opts.add_argument('--skip_mcall', '-Sm',
                   action= 'store_true',
                   help='''Do not call methylation on the aligned BAM file.
                   ''')

wflow_opts.add_argument('--skip_report', '-Sr',
                   action= 'store_true',
                   help='''Do not call produce QC report.
                   ''')

wflow_opts.add_argument('--rrbs',
                   action= 'store_true',
                   help='''Trim reads from RRBS libraries prepared with MspI. This option passed to
trim_galore. From trim_galore: Specifies that the input file was an MspI digested RRBS sample (recognition
site: CCGG). Sequences which were adapter-trimmed will have a further 2 bp
removed from their 3' end. This is to avoid that the filled-in C close to the
second MspI site in a sequence is used for methylation calls. Sequences which
were merely trimmed because of poor quality will not be shortened further.
                   ''')

# -----------------------------------------------------------------------------

prog_opts= parser.add_argument_group('Program options', '''
Further options passed to the programs in the pipeline''')

prog_opts.add_argument('--trim_galore_opts',
                   default= '',
                   help='''String of options to pass to trim_galore. Do not include
`-o/--output_dir` since output directory will be `--outdir`.
                   ''')

prog_opts.add_argument('--bismark_opts',
                   default= '',
                   help='''String of options to pass to bismark. Do not include
`-o/--output_dir` since output directory will be `--outdir`.
                   ''')

prog_opts.add_argument('--mpileup_opts',
                   default= '',
                   help='''String of options to pass to samtools mpileup.
Do not include options -d (set to 1000000), -Q (set to 0), -B, -f (use --ref as
reference fasta) and -l (use --listpos instead).
                   ''')


# -----------------------------------------------------------------------------

prog_paths= parser.add_argument_group('Program paths', '''
Paths to executables.''')

#prog_paths.add_argument('--samtools_path',
#                   default= '',
#                   help='''Path to samtools. Default is '' (samtools is
#on PATH)
#                   ''')

prog_paths.add_argument('--trim_galore_path',
                   default= None,
                   help='''Path to trim_galore, default is to use trim_galore
shipped with the package. Use '' to chose the one on your path.
                   ''')

prog_paths.add_argument('--bismark_path',
                   default= None,
                   help='''Path to bismark, default is to use trim_galore
shipped with the package. Use '' to chose the one on your PATH.
                   ''')

prog_paths.add_argument('--clipoverlap_path',
                   default= '',
                   help='''Path to bam clipOverlap. Default is '' (clipOverlap is
on PATH)
                   ''')

prog_paths.add_argument('--rscript_path',
                   default= '',
                   help='''Path to Rscript. Default is '' (Rscript is
on PATH)
                   ''')


# -----------------------------------------------------------------------------

def main():
    """_MEMO_: Each step takes as file input the variable `workingFileList` (list of
    one or two input file names) and returns `workingFileList` updated with its output.
    
    Output names are <prefix>.<prog>.<ext>[.gz]

    The output dict from each step is appended to logBuilder
    """
    args = parser.parse_args()
    print('\n\nbsExpress version %s\n' %(VERSION))
    if not args.outdir:
        args.outdir= args.prefix
    if not os.path.exists(args.outdir):
        try:
            os.makedirs(args.outdir)
        except OSError:
            sys.exit('\n\nI cannot create output dir "%s"\n' %(args.outdir))
    
    checkList= validate_args.check_settings(args)
    if args.check:
        sys.exit()
    wdir= oxbs_qc_func.get_wdir(args.outdir)
    
    validate_args.validate_prefix(args.prefix)
    validate_args.validate_trim_galore_opts(args.trim_galore_opts)
    validate_args.validate_bismark_opts(args.bismark_opts)
    validate_args.validate_mpileup_opts(args.mpileup_opts)
    validate_args.validate_listpos(args.listpos)
    inputType= validate_args.validate_input(args.input)
    logBuilder= []
    print('PIPELINE:')
    workingFileList= [x for x in args.input]
    
    if not args.skip_fastqc:
        print('\n-- Executing FastQC')
        p= oxbs_qc_func.task_fastqc(infile= args.input, outdir= args.outdir, opts= '', verbose= args.verbose)
        logBuilder.append(p)    

    if inputType == 'raw':
        if len(workingFileList) == 1:
            "Do not clip reads if only one fastq file is given"
            args.skip_clip= True
        if not args.skip_trim:
            """trim_galore.
            _MEMO_: Output fastq file must not include PATH. Use opts= '-o' to assign output dir.
            This should be changed!
            """
            print('\n-- Executing trim_galore')
            outFiles= []
            opts= shlex.split(args.trim_galore_opts) + ['-o', args.outdir] + [' --stringency 13'] +  [' --quality 20'] + [' --three_prime_clip_R1 50']
            opts= ' '.join(opts)
            n= 1
            for x in args.input:
                outFiles.append(args.prefix + '.R' + str(n) + '.trim.fq.gz')
                n += 1

            if args.trim_galore_path is None:
                tgpath= oxbs_qc_func.mypath
            else:
                tgpath= args.trim_galore_path

            p= oxbs_qc_func.task_trim_galore(inFastq= args.input, outFastq= outFiles, rrbs= args.rrbs, tmpdir= None, opts= opts,
                path= tgpath, encoding= args.encoding, verbose= args.verbose)
            logBuilder.append(p)
            workingFileList= [os.path.join(args.outdir, x) for x in outFiles]

        ## With --skip_shorten skip this block altogether
        if not args.skip_shorten:
            if args.maxlen > 0 or os.path.split(args.ref)[1] == CONTROL_FASTA:
                if os.path.split(args.ref)[1] == CONTROL_FASTA and args.maxlen is None:
                    shortLen= 60
                else:
                    shortLen= args.maxlen
                print('\n-- Shortening reads')
                outFiles= []
                n= 1
                for x in workingFileList:
                    shortF= os.path.join(args.outdir, args.prefix + '.R' + str(n) + '.short.fq.gz')
                    xout= oxbs_qc_func.task_shorten_fastq(x, shortF, shortLen, verbose= args.verbose)
                    outFiles.append(shortF)
                    logBuilder.append({'cmd': xout['cmd']})
                    n += 1
                workingFileList= [x for x in outFiles]

        if not args.skip_aln:
            print("\n-- Aligning reads")
            outSam= os.path.join(args.outdir, args.prefix + '.sam')
            if not args.ref:
                sys.exit('I cannot aling reads without a reference file and directory!')
            if not os.path.isfile(args.ref):
                sys.exit('''Reference fasta file "%s" not found.''' %(args.ref))

            if args.bismark_path is None:
                bspath= oxbs_qc_func.mypath
            else:
                bspath= args.bismark_path
            opts= shlex.split(args.bismark_opts) + ['-o', args.outdir] + [' --sam']
            opts= ' '.join(opts)
#            opts= ' '.join(['-o', args.outdir])
            p= oxbs_qc_func.task_bismark_aln(inFastq= workingFileList, outSam= os.path.split(outSam)[1], ref= args.ref,
                tmpdir= None, opts= opts, path= bspath, encoding= args.encoding, verbose= args.verbose)
            logBuilder.append(p)
            workingFileList= [outSam]            
                
    ## From now on you have BAM files.
    ## The input file is always a list of length 1 (workingFileList)
    if not args.skip_sam2bam:
        ## Convert SAM to BAM
        print('\n-- Converting SAM to BAM, sorting and indexing')
        p= oxbs_qc_func.task_sam2bam(workingFileList[0], rmSam= True, verbose= args.verbose)
        logBuilder.append(p)
        workingFileList= [p['bam']]

    if args.mark_duplicates:
        print('\n-- Marking duplicates')
        outBam= outBam= os.path.join(args.outdir, args.prefix + '.mdup.bam') ## re.sub('\.bam$', '.mdup.bam', workingFileList[0])
        p= oxbs_qc_func.task_markDuplicates(workingFileList[0], outBam, verbose= args.verbose) ## outBam= None, opts= None, path= mypath, java_opts= ''
        logBuilder.append(p)
        workingFileList= [p['bam']]
        
    if not args.skip_clip:
        print("\n-- Clipping read overlap")
        outBam= os.path.join(args.outdir, args.prefix + '.clip.bam')
        p= oxbs_qc_func.task_clipOverlap(inBam= workingFileList[0], clippedBam= outBam, clipOverlapPath= args.clipoverlap_path, verbose= args.verbose)
        logBuilder.append(p)
        workingFileList= [outBam]
    
    if not args.skip_mcall:
        print("\n-- Calling methylation")

        if len(workingFileList) == 1:
           print "Single End as just one fastq given"

        outBdg= os.path.join(args.outdir, args.prefix + '.mcall.bdg.gz')
        if args.listpos is not None:
            l= ' -l ' + args.listpos
        else:
            l= ''
        p= oxbs_qc_func.task_callMethylation(inBam= workingFileList[0], outBdg= outBdg,
            opt_ref= args.ref, opt_l= args.listpos, verbose= args.verbose, opt_relation= len(workingFileList))
        logBuilder.append(p)
        workingFileList= [outBdg]

    if not args.skip_report:
        print("\n-- Reporting results")
        reftxt= os.path.splitext(args.ref)[0] + '.txt'
        if not os.path.isfile(reftxt):
            sys.exit('Reference txt file %s not found' %(ref))
        p= oxbs_qc_func.task_oxbs_report(inBdg= workingFileList[0], ref= reftxt, output_dir= args.outdir, prefix= args.prefix, rpath= args.rscript_path, verbose= args.verbose)
        logBuilder.append(p)
    
#    print('-' * 80)
#    print('EXECUTED COMMANDS:\n')
#    for x in logBuilder:
#        print(x['cmd'] + '\n')
    print('-' * 80)
    print('\n\n-- Successfully completed --\n\n')
           
#    print(logBuilder)    
if __name__ == '__main__':
    main()
    sys.exit()
    
