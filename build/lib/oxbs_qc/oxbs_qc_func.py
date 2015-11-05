import sys
import os
import subprocess
import tempfile
import re
import gzip
import shlex
import shutil
import oxbs_qc_func
from distutils import spawn

"""TODO:

* Function to check executable (file) exists either on PATH if no path is given
or on the specified path.

* Get report files from trim_galore and bismark

* Allow custom scirpts (e.g .R) to be found on path other than mypath

* Change direcotry assignment in trim-galore & bismark

* Auto detect shortest refernce length ?

* Add verbose opt to all task_* so that cmd is printed before it is executed.

"""

mypath= os.path.abspath(os.path.split(oxbs_qc_func.__file__)[0])

class FastqcException(Exception):
    pass

class TrimGaloreException(Exception):
    """Base class for exceptions raised by trim_galore."""
    pass

class ShortenFastqException(Exception):
    pass

class BismarkException(Exception):
    pass

class GetWdirException(Exception):
    pass

class Sam2BamException(Exception):
    pass

class MarkDuplicatesException(Exception):
    pass

class ClipOverlapException(Exception):
    pass

class CallMethylationException(Exception):
    pass

class RoxbsReport(Exception):
    pass

def get_wdir(wdir):
    """Set a working directory.
    Possibly TODO: use tempfile.mkdir if wdir is None 
    """
    if wdir == '' or type(wdir) != str:
        raise GetWdirException('\n\nwdir must be a non-empty str. Got "%s"\n' %(wdir))
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    return(wdir)

def task_fastqc(infile, outdir, opts, path= mypath, verbose= False):
    """Execute FastQC
    infile:
        List of input fastq of bam. Can be of length 2 for paired fastq files
    outdir:
        Output dir for fastqc report
    opts:
        String of further options to pass to fastqc (not implemented).
    path:
        Path to fastqc
    """
    if type(infile) != list:
        raise FastqcException('Input file(s) must be a list e.g. ["myfile.fq"].')
    for x in infile:
        if not os.path.isfile(x):
            raise FastqcException('Cannot find input file %s' %(x))
    if not os.path.isdir(outdir):
        raise FastqcException('%s is not a directory' %(outdir))

    fastqc= os.path.join(path, 'FastQC/fastqc')
    if not os.path.isfile(fastqc):
        raise FastqcException('\n\nCannot find %s\n' %(fastqc))
    cmd= ' '.join(['perl', fastqc, '-t', str(len(infile)), '-o', outdir,' '.join(infile)])
    if verbose:
        print(cmd)
    p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
    stdout, stderr= p.communicate()
    if p.returncode != 0:
        print(stderr)
        print('Exit code %s' %(p.returncode))
        raise FastqcException('\n\nError produced by fastqc\n')

    return({'cmd':cmd, 'stderr': stderr, 'stdout': stdout, 'returncode': p.returncode})
    
    
def task_trim_galore(inFastq, outFastq, rrbs= False, tmpdir= None, opts= '', path= mypath, encoding= None, verbose= False):
    """Execute trim_galore on the input FastqKit.
    inFastq:
        List of length one or two. If length 2 it is taken as paired-end
    tmpdir:
        A working directory where to dump output files. Inside this dir python will
        make a tempdirectory and the ouput of trim_galore put there
    outFastq:
        List of length 1 or 2 for the name of the output fastq file(s). DO NOT include
        path (it will be stripped anyway).
    rrbs:
        True/False to enable --rrbs opt
    opts:
        String of options passed to trim_galore
    path:
        path to trim_galore
    echo:
        Only print the command that would be executed.
    return:
        Dictionary with:
        {'cmd': 'command-executed-by-subprocess',
         'fastq': ['fastq1', 'fastq2'],
         'report': 'trim_galore-report'}
    """
    if type(inFastq) != list:
        raise TrimGaloreException('\n\nInput fastq files must be a *list*. Got %s\n' %(type(inFastq)))
    if type(outFastq) != list:
        raise TrimGaloreException('\n\nOutput fastq files must be a *list*. Got %s\n' %(type(outFastq)))
    if len(inFastq) < 1 or len(inFastq) > 2:
        raise TrimGaloreException('\n\nExpected a list of 1 or 2 fastq files. Got %s.\n' %(len(inFastq)))
    for x in inFastq:
        if not os.path.exists(x):
            raise TrimGaloreException('\n\nFile %s not found\n' %(x))
    for i in range(0, len(outFastq)):
        fq= outFastq[i]
        p,f= os.path.split(fq)
        if p != '':
            raise TrimGaloreException('\n\nOutput fastq file(s) must not have directory path. Got "%s"\n' %(fq))
    if len(inFastq) != len(outFastq):
        raise TrimGaloreException('\n\n%s fastq files in input but %s found in output arg\n' %(len(inFastq), len(outFastq)))
    if len(set(inFastq)) != len(inFastq):
        raise TrimGaloreException('\n\nInvalid input (same file given twice?): %s\n' %(inFastq))
    if len(set(outFastq)) != len(outFastq):
        raise TrimGaloreException('\n\nInvalid input (same file given twice?): %s\n' %(outFastq))
    ## Set working dir:
    if tmpdir is not None:
        wdir= get_wdir(tmpdir)
    else:
        wdir= tempfile.mkdtemp(prefix= 'trim_galore_wdir')
    tg_outdir= tempfile.mkdtemp(prefix= 'tmp_trim_galore_', dir= wdir)
    opts= shlex.split(opts)
    
    if '--dont_gzip' in opts:
        sys.exit('Sorry, trim_galore option --dont_gzip is not supported.')
    
    if '--gzip' not in opts:
        opts.append('--gzip')
    
    ## Quality Encoding
    if encoding == 'Sanger':
        opts.append('--phred33')
    elif encoding == 'Illumina':
        opts.append('--phred64')
    elif '--phred33' not in opts and '--phred64' not in opts:
        encoding= get_fastq_encoding(inFastq[0])
        if encoding == 'Sanger':
            opts.append('--phred33')
        elif encoding == 'Illumina 1.5+':
            opts.append('--phred64')
        else:
            sys.exit('Unable to determine quality encoding for %s. Consider using option --encoding' %(inFastq))
    else:
        sys.exit("Unable to set encoding")

    ## This otuput dir will NOT be passed to trim_galore. Instead the output files
    ## from trim_galore will be moved there afterwards.
    if '-o' in opts:
        i= opts.index('-o')
        final_outdir= opts[i+1]
        del opts[i:(i+2)]
    elif '--output_dir' in opts:
        i= opts.index('--output_dir')
        final_outdir= opts[i+1]
        del opts[i:(i+2)]
    else:
        final_outdir= '.'
    ## Handling gzip
    ## If gzip was given to trim_galore add *.gz to output name.
    ## If --dont_gzip id given remove *.gz from output names.
    ## -----------------
    if '--gzip' in opts:
        outFastq= [x.rstrip('.gz') + '.gz' for x in outFastq]
    if '--dont_gzip' in opts:
        outFastq= [x.rstrip('.gz') for x in outFastq]
    
    outFastqFull= [os.path.join(final_outdir, x) for x in outFastq]

    if '--paired' not in opts and len(inFastq) > 1:
        opts.append('--paired')
    else:
        pass

    if rrbs and '--rrbs' not in opts:
        opts.append('--rrbs')
        
    opts= ' '.join(opts)

    trim_galore= os.path.join(path, 'trim_galore')

    inFastq= [os.path.abspath(x) for x in inFastq]
    fastqName= [os.path.split(x)[1] for x in inFastq]
    ## Make a symlink from input fastq to tmp dir to avoid parallel runs of trim_galore
    ## to overwrite each other
    cmd_sym= []
    for fq in inFastq:
        cmd_sym.append('ln -s ' + fq + ' . ')
    cmd_sym= ' ; '.join(cmd_sym)
    cmd_rm_sym= []

    ## Remove symlinks at the end of execution
    for fq in fastqName:
        cmd_rm_sym.append('rm ' + fq)
    cmd_rm_sym= ' ; '.join(cmd_rm_sym)

    ## You need to be in the output dir or trim_galore will throw an error with --rrbs
    cmd_cd= ' '.join(['cd', tg_outdir])

    ## trim_galore command line    
    cmd_trim_galore= ' '.join(['perl', trim_galore, opts, ' '.join(fastqName)])
    
    cmd= ';\n'.join(['set -e', cmd_cd, cmd_sym, cmd_trim_galore, cmd_rm_sym])
    if verbose:
        print(cmd)
    p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
    stdout, stderr= p.communicate()
    if p.returncode != 0:
        print(stderr)
        print('Exit code %s' %(p.returncode))
        raise TrimGaloreException('\n\nError produced by trim_galore\n')
    ## Now get the names of the output fastq file(s) and rename them to outFastq
    ## -------------------------------------------------------------------------
    tg_output_files= sorted(os.listdir(tg_outdir))
    fqz= [x for x in tg_output_files if x.endswith('.fq.gz')]
    fq= [x for x in tg_output_files if x.endswith('.fq')]
    ## Make sure you got the right number of files.
    if len(fqz) > 2 or len(fq) > 2:
        raise TrimGaloreException('\n\nToo many fastq files found in output: %s, %s\n' %(fqz, fq))
    if len(fqz) > 0 and len(fq) > 0:
        raise TrimGaloreException('\n\nFound both gzipped and unzipped fastq files: %s, %s\n' %(fqz, fq))
    if len(fqz) != len(inFastq) and len(fq) != len(inFastq):
        if len(fqz) == 0:
            nout= len(fq)
        else:
            nout= len(fqz)
        raise TrimGaloreException('\n\nInconsistent number of output fastq files found: %s input, %s output\n' %(len(inFastq), len(nout)))
    tg_fq= [os.path.join(tg_outdir, x) for x in fqz + fq]
    for old, new in zip(tg_fq, outFastqFull):
        shutil.move(old, new)
    
    ## Get report:
    ## -----------
    tg_report=       [x for x in tg_output_files if x.endswith('_trimming_report.txt')]
    tg_report_tmp=   [os.path.join(tg_outdir, x) for x in tg_report] 
    tg_report_final= [os.path.join(final_outdir, x) for x in tg_report]
    if len(tg_report) != len(inFastq):
        raise TrimGaloreException('\n\n%s txt files found in trim_galore output dir. Expected %s (one per input).\n' %(len(tg_report), len(inFastq)))
    for src,dest in zip(tg_report_tmp, tg_report_final):
        shutil.move(src, dest)

    ## Exit.
    ## -----
    if len(os.listdir(tg_outdir)) != 0:
        raise TrimGaloreException('\n\nTemporary output dir from trim_galore should be empty now! Found: %s\n' %(sorted(os.listdir(tg_outdir))))
    if tmpdir is None:
        shutil.rmtree(wdir)
    return({'fastq': outFastqFull, 'cmd': cmd, 'report': tg_report_final, 'stderr': stderr, 'stdout': stdout, 'returncode': p.returncode})

def task_shorten_fastq(fq_in, fq_out, mlen, path= mypath, verbose= False):
    """Replace fastx_trimmer with ShortenFastq.jar: Trim reasd from start to given
    length. See help from ShortenFastq.jar
    fq_in, fq_out:
        Name of fastq files for input and output.
    mlen:
        Int
    path:
        Path to ShortenFastq.jar
    Return:
        Dict with executed command.
    """
    shortenFastq= os.path.join(path, 'ShortenFastq.jar')
    if not os.path.isfile(shortenFastq):
        raise ShortenFastqException('\n\n%s not found!\n' %(shortenFastq))
    if not os.path.isfile(fq_in):
        raise ShortenFastqException('\n\n%s not found!\n' %(fq_in))
    if mlen < 0:
        raise ShortenFastqException('\n\nRead length must be >0. Got %s\n' %(mlen))
    try:
        x= open(fq_out, 'w').close()
    except IOError:
        raise ShortenFastqException('\n\nI cannot write to file %s!\n' %(fq_out)) 
    os.remove(fq_out)
    
    cmd= ' '.join(['java -jar', shortenFastq, fq_in, fq_out, str(mlen)])
    if verbose:
        print(cmd)
    p= subprocess.Popen(cmd, shell= True, stderr= subprocess.PIPE, stdout= subprocess.PIPE)
    stdout, stderr= p.communicate()
    if p.returncode != 0:
        print(stderr)
        print('Exit code %s' %(p.returncode))
        print(cmd)
        raise ShortenFastqException('\n\nError while executing task_shorten_fastq\n')
        
    return({'cmd':cmd, 'out': fq_out, 'stderr': stderr, 'stdout': stdout, 'returncode': p.returncode})

def task_bismark_aln(inFastq, outSam, ref, tmpdir= None, opts= '', path= mypath, verbose= False, encoding= None):
    """Execute bismark alignment.
    inFastq:
        List of len 1 or 2 of fastq files. If 2 they are taken as paired.
    outSam:
        Output name for output sam file. Do not inlcude path (it will be stripped).
        Must end in *.sam
    ref:
        Path to reference file (not just dir).
    tmpdir:
        Working dir for bismark. It is NOT where you find the output. Output is given
        by `bismark --output_dir`
    opt:
        String of options to be passed to bismark
    path:
        Path to bismark
    """
    opts= shlex.split(opts)
    
    inFastq= [os.path.abspath(x) for x in inFastq]
    ref= os.path.abspath(ref)

    if type(inFastq) != list:
        raise BismarkException('\n\nInput fastq files must be a *list*. Got %s\n' %(type(inFastq)))
    if len(inFastq) < 1 or len(inFastq) > 2:
        raise BismarkException('\n\nExpected a list of 1 or 2 fastq files. Got %s.\n' %(len(inFastq)))
    if len(set(inFastq)) != len(inFastq):
        raise BismarkException('\n\nInput fastq files are the same! Got %s.\n' %(inFastq))
    for x in inFastq:
        if not os.path.exists(x):
            raise BismarkException('\n\nFile %s not found\n' %(x))
    if not os.path.isfile(ref):
        raise BismarkException('\n\nPath to reference FASTA %s not found\n' %(ref))
    if not outSam.endswith('.sam'):
        raise BismarkException('\n\nOutput sam file must have extension .sam. Got %s\n' %(outSam))
    if os.path.split(outSam)[0] != '':
        raise BismarkException('\n\nOutput sam name must not inlcude directory (use opts= "bismark --output_dir ..."). Got %s\n' %(outSam))
    
    ref= os.path.split(os.path.abspath(ref))[0] ## Get only fll path to ref FASTA dir
    
    ## Set working dir:
    if tmpdir is None:
        wdir= tempfile.mkdtemp(prefix= 'bismark_wdir')
    else:
        wdir= get_wdir(tmpdir)
    bis_outdir= tempfile.mkdtemp(prefix= 'tmp_bismark_wdir_', dir= wdir)
    if '--bowtie2' not in opts:
        opts.append('--bowtie2')

    bismark= os.path.join(path, 'bismark')

    ## Get encoding, use only first FASTQ if a pair:

    if encoding == 'Sanger':
        opts.append('--phred33-quals')
    elif encoding == 'Illumina':
        opts.append('--phred64-quals')
    elif '--phred33-quals' not in opts and '--phred64-quals' not in opts and '--solexa-quals' not in opts and '--solexa1.3-quals' not in opts:
        encoding= get_fastq_encoding(inFastq[0])
        if encoding == 'Sanger':
            opts.append('--phred33-quals')
        elif encoding == 'Illumina 1.5+':
            opts.append('--phred64-quals')
        else:
            sys.exit('Unable to determine quality encoding for %s. Consider setting option --encoding' %(inFastq))
    else:
        sys.exit('Unable to detect encoding');
        
    ## This otuput dir will NOT be passed to bismark. Instead the output files
    ## from bismark will be moved there afterwards.
    if '-o' in opts:
        i= opts.index('-o')
        final_outdir= opts[i+1]
        del opts[i:(i+2)]
    elif '--output_dir' in opts:
        i= opts.index('--output_dir')
        final_outdir= opts[i+1]
        del opts[i:(i+2)]
    else:
        final_outdir= '.'

    ## Single or paired end
    if len(inFastq) == 1:
        fq= [inFastq[0]]
    else:
        fq= ['-1', inFastq[0], '-2', inFastq[1]]

    cmd_cd= ' '.join(['set -e', '\n', 'cd', bis_outdir]) ## cd to bismark wdir
    cmd_bismark= ['perl'] + [bismark] + ['-o .'] + opts + [ref] + fq
    cmd_bismark= ' '.join(cmd_bismark)
    cmd= cmd_cd + '\n' + cmd_bismark
    if verbose:
        print(cmd)
    p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
    stdout, stderr= p.communicate()
    if p.returncode != 0:
        print(stderr)
        print('Exit code %s' %(p.returncode))
        raise BismarkException('\n\nError produced by bismark\n')

    ## Bismark OUTPUT files
    output_files= sorted(os.listdir(bis_outdir))
    output_files_path= [os.path.join(bis_outdir, x) for x in output_files]
    sam_output= [x for x in output_files if x.endswith('.sam') or x.endswith('.bam')]
    if len(sam_output) != 1:
        raise BismarkException('\n\nExpected 1 sam output. Found %s (%s)\n' %(len(sam_output), sam_output))
    sam_output= sam_output[0]
    
## Not necessary with bismark 0.10+
    ## If PE clean read names to remove /1 /2 suffix:
    ## ----------------------------------------------
#    if len(inFastq) == 2:
#        print("\nRemoving /1 and /2 from read names...")
#        CleanReadNames= os.path.join(path, 'CleanReadNames.jar')
#    
#        if not os.path.isfile(CleanReadNames):
#            raise BismarkException('\n\n%s not found!\n' %(CleanReadNames))
#    
#        input= os.path.join(bis_outdir, sam_output)
#        output= os.path.join(bis_outdir, "CleanReadNames." + sam_output)
#    
#        cmd_clean= ' '.join(["java -jar", CleanReadNames, input, output])
#        if verbose:
#            print(cmd_clean)
#        cmd+=cmd_clean
#        p= subprocess.Popen(cmd_clean, shell= True, stderr= subprocess.PIPE, stdout= subprocess.PIPE)
#        stdout, stderr= p.communicate()
#        if p.returncode != 0:
#            print(stderr)
#            print('Exit code %s' %(p.returncode))
#            print(cmd_clean)
#            raise BismarkException('\n\nError while executing CleanReadNames\n')
#        if verbose:
#            print("\n%s\nrenamed to\n%s" %(output, input))
#        os.rename(output, input)


    ## Move output. Sam will be renamed. Everything else stays the same
    ## ----------------------------------------------------------------
    sam= os.path.join(final_outdir, outSam)
    shutil.move(os.path.join(bis_outdir, sam_output), sam)
    all_out= [sam]
    for x in output_files_path:
        if x.endswith('.sam'):
            continue
        dst= os.path.join(final_outdir, os.path.split(x)[1])
        shutil.move(x, dst)
        all_out.append(dst)
    if tmpdir is None:
        shutil.rmtree(wdir)
    return({'cmd':cmd, 'stderr': stderr, 'stdout': stdout, 'returncode': p.returncode, 'sam': sam, 'all_out': all_out})

def task_sam2bam(inSam, rmSam= True, verbose= False):
    """With sam input convert sam to bam, sort and index. For bam input sort and index only.
    inSam:
        sam/bam file to convert
    rmSam:
        Remove original SAM? NB: If input is BAM the unsorted BAM is going to be overwritten!
    return:
        Dictionary with list of cammands executed ('cmd_list'), output bam file ('bam')
    """
    if not os.path.exists(inSam):
        raise Sam2BamException('\n\nCould not find input file "%s"\n' %(inSam))
    if inSam.endswith('.sam'):
        isSam= True
    elif inSam.endswith('.bam'):
        isSam= False
    else:
        raise Sam2BamException('\n\nInvalid extension in %s. Expected .sam or .bam\n' %(inSam))
    
    if isSam:
        outBamUnsorted= re.sub('\.sam$', '.unsorted.bam', inSam)
        ## Check sam has no alignment
        fin= open(inSam)
        hasReads= False
        for line in fin:
            if not line.startswith('@'):
                hasReads= True
                break
        fin.close()
        if not hasReads:
            sys.exit('\n\nNo reads in file "%s"\n. Are the reads longer than reference sequence(s)?\n' %(inSam))        
    else:
        outBamUnsorted= inSam
    outBamPrefix= re.sub('\.bam$|\.sam$', '.sorted', inSam)
    outBam= re.sub('\.bam$|\.sam$', '', inSam) + '.bam'
    samtools= 'samtools'
    cmds= []
    if isSam:
        cmds.append(' '.join([samtools, 'view', '-h', '-S', '-b', '-o', outBamUnsorted, inSam]))
    cmds.append(' '.join([samtools, 'sort',  outBamUnsorted, outBamPrefix]))
    cmds.append(' '.join(['mv', outBamPrefix + '.bam', outBam]))
    cmds.append(' '.join([samtools, 'index', outBam]))
    if isSam:
        cmds.append(' '.join(['rm', outBamUnsorted]))
    if rmSam and isSam:
        cmds.append(' '.join(['rm', inSam]))
    if verbose:
        print('\n'.join(cmds))    
    for cmd in cmds:
        ## Execute commands one by one
        if verbose:
            print("Executing:\n" + cmd)
        p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
        stdout, stderr= p.communicate()
        if p.returncode != 0:
            print(stderr)
            print('Exit code %s' %(p.returncode))
            print(cmd)
            raise Sam2BamException('\n\nError produced while converting sam to bam\n')
    return({'cmd': ';\n'.join(cmds), 'stderr': stderr, 'stdout': stdout, 'returncode': p.returncode, 'bam': outBam})

def task_markDuplicates(inBam, outBam= None, opts= None, path= mypath, java_opts= '', verbose= False):
    """Execute picard MarkDuplicates on input *sorted* bam files
    inBam, outBam:
        Name of input and output bam files. Output name defaults to input with
        .bam replaced by .mdup.bam
    opts:
        Further options passed to MarkDuplicates. Defualt None sets some sensible
        settings
    path:
        Path to MarkDuplicates.jar.
    java_opts:
        Options passed to java jvm.
    """
    if not os.path.isfile(inBam):
        raise MarkDuplicatesException('\n\nFile "%s" not found\n' %(inBam))
    if outBam is not None and not outBam.endswith('.bam'):
        raise MarkDuplicatesException('\n\nOutput file will be in BAM format and must have .bam extension. Got "%s"\n' %(outBam))
    
    if not outBam:
        outBam= re.sub('\.bam$', '.mdup.bam', inBam)

    metrics_file= re.sub('\.bam$', '.metrics.txt', outBam)
    
    markDuplicates= os.path.join(path, 'MarkDuplicates.jar')
    if not os.path.isfile(markDuplicates):
        raise MarkDuplicatesException('\n\n%s not found!\n' %(markDuplicates))
    
    if not opts:
        opts= 'VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true'
    
    input= 'INPUT=' + inBam
    output= 'OUTPUT=' + outBam
    metrics= 'METRICS_FILE=' + metrics_file
    opts= opts ## 
    
    cmd= ' '.join(['java -Xmx1g', java_opts, '-jar', markDuplicates, input, output, metrics, opts])
    if verbose:
        print(cmd)
    p= subprocess.Popen(cmd, shell= True, stderr= subprocess.PIPE, stdout= subprocess.PIPE)
    stdout, stderr= p.communicate()
    if p.returncode != 0:
        print(stderr)
        print('Exit code %s' %(p.returncode))
        print(cmd)
        raise MarkDuplicatesException('\n\nError while executing task_markDuplicates\n')
        
    return({'cmd':cmd, 'bam': outBam, 'metrics_file': metrics_file, 'stderr': stderr, 'stdout': stdout})


def task_clipOverlap(inBam, clippedBam, clipOverlapPath= '', path= mypath, verbose= False):
    """Excute clipOverlap. Input file is passed thorugh samtools to convert
    bam -> sam, then cleanReadNames to remove trailing /1 and /2 the clipOverlap
    inBam:
        Input bam file to clip. Must be sorted.
        See http://genome.sph.umich.edu/wiki/BamUtil:_clipOverlap
    clippedBam:
        Output clipped bam
    path:
        Path `bam` utils
    path:
        Path to cleanReadNames.py (should be in site-packages with the others)
    return:
        dict with {'cmd': executed-command}
        Side effect: Clipped bam file, indexed and 
    
    TODO: cleanReadNames shouldn;t be necessary with bismark 0.10+
    """
    if not os.path.isfile(inBam):
        raise ClipOverlapException('\n\nFile "%s" not found\n' %(inBam))
    if not clippedBam.endswith('.bam'):
        raise ClipOverlapException('\n\nOutput file is going to be BAM and must have .bam extension. Got "%s"\n' %(clippedBam))
    
    bamutils= os.path.join(clipOverlapPath, 'bam')
    samtools= 'samtools'
    cleanReadNames= os.path.join(path, 'cleanReadNames.py')
    
    if inBam.endswith('.bam'):
        view= 'view -h'
    elif inBam.endswith('.sam'):
        view= 'view -S -h'
    else:
        raise ClipOverlapException('\n\nInvalid input extension. Expected .sam or .bam got %s\n' %(inBam))
    
    cmd_clip = getBash() + " -c '"
    cmd_clip += ' '.join(['set -o pipefail; set -e;\n', samtools, view, inBam, '| python', cleanReadNames, '-i - |', bamutils, 'clipOverlap  --stats --in -.sam',  '--out', clippedBam])
    cmd_clip += "'" ## Close single quotes
    
    if verbose:
        print(cmd_clip)
    p= subprocess.Popen(cmd_clip, shell= True, stderr= subprocess.PIPE, stdout= subprocess.PIPE)
    stdout, stderr= p.communicate()
    if p.returncode != 0:
        print(stderr)
        print('Exit code %s' %(p.returncode))
        print(cmd_clip)
        raise ClipOverlapException('\n\nError while executing task_clipOverlap\n')
    statsFile= open(clippedBam + '.clipStats', 'w')
    for line in stderr:
        statsFile.write(line)
    statsFile.close()
    
    cmd_idx=  ' '.join([samtools, 'index', clippedBam])
    if verbose:
        print(cmd_idx)
    p= subprocess.Popen(cmd_idx, shell= True, stderr= subprocess.PIPE, stdout= subprocess.PIPE)
    stdout, stderr= p.communicate()
    if p.returncode != 0:
        print(stderr)
        print('Exit code %s' %(p.returncode))
        print(cmd_idx)
        raise ClipOverlapException('\n\nError while executing %s\n' %(cmd_idx))
        
    return({'cmd': cmd_clip + ';\n' + cmd_idx, 'stderr': stderr, 'stdout': stdout, 'returncode': p.returncode})

def task_callMethylation(inBam, outBdg,
        opt_ref, opt_l= None, opt_A= '', opt_samargs= '', opt_region= '', 
        bam2methylation_path= mypath, verbose= False, opt_relation= ''):
    """
    inBam:
        Input bam sorted
    outBdg:
        Output file (bedgraph-like format)
    opt_*:
        Options passed to bam2methyaltion.py
    *_path;
        Path to scirpts & binaries
    Return:
        Dict: {'cmd': executed-command}
    """
    bam2met= os.path.join(bam2methylation_path, 'bam2methylation.py')
    if not os.path.exists(inBam):
        raise CallMethylationException('\n\nFile %s not found\n' %inBam())
    pathToOut= os.path.split(outBdg)[0]
    if not os.path.exists(pathToOut) and pathToOut != '':
        raise CallMethylationException('\n\nPath to output file %s does not exists!\n' %(outBdg))

    ## Get mandatory args:
    opt_input= '-i ' + inBam
    opt_ref= '-r ' + opt_ref

    ## Get optional args to bam2met:
    if opt_l != '' and opt_l is not None:
        opt_l= '-l ' + opt_l
    else:
        opt_l= ''
    if opt_A != '':
        opt_A= '-A'
    if opt_samargs != '':
        opt_samargs= '-s ' + opt_samargs
    if opt_region != '':
        opt_region= '--region ' + opt_region
    if opt_relation == 1:
        opt_relation = '--se '
    else:
       opt_relation = ''

         
    cmd= ' '.join(['set -e;\n', 'python', bam2met, opt_input, opt_ref, opt_l, opt_A, opt_samargs, opt_region, opt_relation, '| gzip >', outBdg])
    
    if verbose:
        print(cmd)
    p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
    stdout, stderr= p.communicate()
    if p.returncode != 0:
        print(stderr)
        print('Exit code %s' %(p.returncode))
        print(cmd)
        raise CallMethylationException('\n\nError while calling methylation\n')    
    return({'cmd':cmd, 'stderr': stderr, 'stdout': stdout, 'returncode': p.returncode})


def task_oxbs_report(inBdg, ref, output_dir, prefix, path= mypath, rpath= '', verbose= False):
    """
    inBdg:
        Input bedgraph-like file as returned by mpileup2methylation.py
    ref:
        Reference file of modification. See R code for expected format.
    path:
        Path to oxbs_report.R
    rpath:
        Path tp Rscript

    MEMO: Commands args to oxbs_report.R:
    input<- opts[1]
    reference<- opts[2]
    output_dir<- opts[3]
    output_prefix<- opts[4]
    """
    if not os.path.isfile(inBdg):
        raise RoxbsReport('\n\nInput file %s does no exist\n' %(inBdg))
    if not os.path.isfile(ref):
        raise RoxbsReport('\n\nReference file %s does no exist\n' %(ref))

    oxbs_report= os.path.join(path, 'oxbs_report.R')
    if not os.path.isfile(oxbs_report):
        raise RoxbsReport('\n\nCannot find %s does no exist\n' %(oxbs_report))
    
    if rpath == '':
        Rscript= spawn.find_executable('Rscript')
    else:
        Rscript= os.path.join(rpath, 'Rscript')
    if Rscript is None or not os.path.isfile(Rscript):
        raise RoxbsReport('\n\nRscript not found in %s\n' %(rpath))
    
    cmd= [Rscript, oxbs_report, inBdg, ref, output_dir, prefix]
#    if output_dir != '':
#        cmd.append(output_dir)
#    if prefix != '':
#        cmd.append(prefix)
    cmd= ' '.join(cmd)
    if verbose:
        print(cmd)
    p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
    stdout, stderr= p.communicate()
    if p.returncode != 0:
        print(stderr)
        print('Exit code %s' %(p.returncode))
        print(cmd)
        raise RoxbsReport('\n\nError executing oxbs_report\n')    
    return({'cmd':cmd, 'stderr': stderr, 'stdout': stdout, 'returncode': p.returncode})

def get_fastq_encoding(fastq):
    """Read fastq file to determine the quality encoding. Fastq file can be either
    plain text or gzipped.
   
    **Returns**:
        String ``Sanger`` or ``Illumina 1.5+`` or ``Undetermined``.
        There is no distinction within Illumina variations and no distictions between Sanger and `Illumina 1.8+ Phred+33`.
        If only ambiguos codes are found (e.g. 'BCDEFGHI') the string ``Undertermined`` is returned

    .. seealso::

        A sample `perl <https://www.uppnex.uu.se/content/check-fastq-quality-score-format>`_ script using the same approach and
        `Wikipedia fastq format <http://en.wikipedia.org/wiki/FASTQ_format#Encoding>`_
       
        SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................
        ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................
        ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................
        .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................
        LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................
        !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
        |                         |    |        |                              |                     |
       33                        59   64       73                            104                   126
        0........................26...31.......40                                
                                 -5....0........9.............................40
                                       0........9.............................40
                                          3.....9.............................40
        0........................26...31........41                              
     
       S - Sanger        Phred+33,  raw reads typically (0, 40)
       X - Solexa        Solexa+64, raw reads typically (-5, 40)
       I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
       J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
           with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
           (Note: See discussion above).
       L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)

    .. note::
       
        * A file is ``Sanger`` if a quality ascii <59 is found and ``Illumina`` if ascii >74 is found. Files where neither of these ranges are found will throw an exception.
        * NB: There is no check of whether the codes are within the known encodings
          (e.g. no error or warning if you the fastq has qualities like 'lmnopqrstuvwxy')    
    """
    if fastq.endswith('.gz'):
        fin= gzip.open(fastq)
    else:
        fin= open(fastq)
    while True:
        fin.readline()
        fin.readline()
        fin.readline()
        qual= fin.readline().strip()
        if qual == '':
            fin.close()
            return('Undetermined')
            ## sys.exit(inspect.stack()[0][3] + ': Encoding of %s could not be determined after having read the whole file' %(fastq))
        encoding= [ord(x) for x in qual]
        for x in encoding:
            if x < 59:
                fin.close()
                return('Sanger')
            if x > 74:
                fin.close()
                return('Illumina 1.5+')

def getBash():
    """This is a hacky function to get the path to bash.
    """
    cmd= 'which bash'
    p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
    stdout, stderr= p.communicate()
    if p.returncode != 0:
        sys.exit("Error excuting %s" %(cmd))
    return(stdout.strip())
