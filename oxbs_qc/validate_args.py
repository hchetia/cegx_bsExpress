"""Functions to validate command line arguments
"""

import os
import shlex
import sys
from distutils import spawn
import oxbs_qc
import re

class ValidArgException(Exception):
    pass

def check_settings(args):
    """args is a parser obejct with all the input options. Check they are all ok.
    check_only:
        If True, sys.exit() after checking
    """
    print('')
    print('-'*80)
    print('CHECK SETTINGS:\n')

    pypath= os.path.abspath(os.path.split(oxbs_qc.__file__)[0])
    checkList= []
    sys.stdout.write('Check output directory "%s" is writable...        ' %(args.outdir))
    passed= os.access(args.outdir, os.W_OK)
    if passed:
        sys.stdout.write('OK\n')
    else:
        sys.stdout.write('FAILED\n')
    checkList.append(('outdir', passed))

    ## samtools
    ## -----------
    tg= spawn.find_executable('samtools')

    if tg is None:
        passed= False
    else:
       passed= os.path.isfile(tg)

    sys.stdout.write('Check samtools "%s"...     ' %(tg))
    if passed:
        sys.stdout.write('OK\n')
    else:
        sys.stdout.write('FAILED\n')
    checkList.append(('samtools', passed))

    ## samtools
    ## -----------
    tg= spawn.find_executable('bedtools')

    if tg is None:
        passed= False
    else:
       passed= os.path.isfile(tg)

    sys.stdout.write('Check bedtools "%s"...     ' %(tg))
    if passed:
        sys.stdout.write('OK\n')
    else:
        sys.stdout.write('FAILED (methyaltion can\'t be called)\n')
    checkList.append(('bedtools', passed))


    ## cutadapt
    ## trim_glore doesn't have a cutadapt path options. So it must be on the PATH
    ## -----------
    tg= spawn.find_executable('cutadapt')
    if tg is None:
        passed= False
    else:
       passed= os.path.isfile(tg)

    sys.stdout.write('Check cutadapt "%s"...     ' %(tg))
    if passed:
        sys.stdout.write('OK\n')
    else:
        sys.stdout.write('FAILED\n')
    checkList.append(('cutadapt', passed))

    ## Trim galore
    ## -----------
    if args.trim_galore_path is None:
        tgpath= os.path.abspath(os.path.split(oxbs_qc.__file__)[0])
    else:
        tgpath= args.trim_galore_path    
    tg= os.path.join(tgpath, 'trim_galore')
    passed= os.path.isfile(tg)

    sys.stdout.write('Check trim_galore "%s"...     ' %(tg))
    if passed:
        sys.stdout.write('OK\n')
    else:
        sys.stdout.write('FAILED\n')
    checkList.append(('trim_galore', passed))

    ## bismark
    ## -------
    if args.bismark_path is None:
        tgpath= os.path.abspath(os.path.split(oxbs_qc.__file__)[0])
    else:
        tgpath= args.bismark_path    
    tg= os.path.join(tgpath, 'bismark')
    passed= os.path.isfile(tg)

    sys.stdout.write('Check bismark "%s"...     ' %(tg))
    if passed:
        sys.stdout.write('OK\n')
    else:
        sys.stdout.write('FAILED\n')
    checkList.append(('bismark', passed))

    ## clipOverlap
    ## -----------
    if args.clipoverlap_path == '':
        tg= spawn.find_executable('bam')
    else:
        tg= os.path.join(args.clipoverlap_path, 'bam')
    
    if tg is None:
        passed= False
    else:
       passed= os.path.isfile(tg)

    sys.stdout.write('Check bam clipOverlap "%s"...     ' %(tg))
    if passed:
        sys.stdout.write('OK\n')
    else:
        sys.stdout.write('FAILED\n')
    checkList.append(('clipOverlap', passed))

    ## R
    ## -----------
    if args.rscript_path == '':
        tg= spawn.find_executable('Rscript')
    else:
        tg= os.path.join(args.rscript_path, 'Rscript')

    if tg is None:
        passed= False
    else:
       passed= os.path.isfile(tg)

    sys.stdout.write('Check R/Rscript "%s"...     ' %(tg))
    if passed:
        sys.stdout.write('OK\n')
    else:
        sys.stdout.write('FAILED\n')
    checkList.append(('R/Rscript', passed))


    ## Custom scripts
    ## --------------
    tg= os.path.join(pypath, 'FastQC/fastqc')
    passed= os.path.isfile(tg)
    sys.stdout.write('Check fastqc "%s"...     ' %(tg))
    if passed:
        sys.stdout.write('OK\n')
    else:
        sys.stdout.write('FAILED\n')
    checkList.append(('fastqc', passed))

    tg= os.path.join(pypath, 'ShortenFastq.jar')
    passed= os.path.isfile(tg)
    sys.stdout.write('Check ShortenFastq.jar "%s"...     ' %(tg))
    if passed:
        sys.stdout.write('OK\n')
    else:
        sys.stdout.write('FAILED\n')
    checkList.append(('ShortenFastq.jar', passed))

    tg= os.path.join(pypath, 'MarkDuplicates.jar')
    passed= os.path.isfile(tg)
    sys.stdout.write('Check MarkDuplicates.jar "%s"...     ' %(tg))
    if passed:
        sys.stdout.write('OK\n')
    else:
        sys.stdout.write('FAILED\n')
    checkList.append(('MarkDuplicates.jar', passed))

    tg= os.path.join(pypath, 'cleanReadNames.py')
    passed= os.path.isfile(tg)
    sys.stdout.write('Check cleanReadNames.py "%s"...     ' %(tg))
    if passed:
        sys.stdout.write('OK\n')
    else:
        sys.stdout.write('FAILED\n')
    checkList.append(('cleanReadNames.py.py', passed))

    #tg= os.path.join(pypath, 'mpileup2methylation.py')
    #passed= os.path.isfile(tg)
    #sys.stdout.write('Check mpileup2methylation.py "%s"...     ' %(tg))
    #if passed:
    #    sys.stdout.write('OK\n')
    #else:
    #    sys.stdout.write('FAILED\n')
    #checkList.append(('mpileup2methylation.py', passed))
   
    tg= os.path.join(pypath, 'bam2methylation.py')
    passed= os.path.isfile(tg)
    sys.stdout.write('Check bam2methylation.py "%s"...     ' %(tg))
    if passed:
        sys.stdout.write('OK\n')
    else:
        sys.stdout.write('FAILED\n')
    checkList.append(('bam2methylation.py', passed))
 
    tg= os.path.join(pypath, 'oxbs_report.R')
    passed= os.path.isfile(tg)
    sys.stdout.write('Check oxbs_report.R "%s"...     ' %(tg))
    if passed:
        sys.stdout.write('OK\n')
    else:
        sys.stdout.write('FAILED\n')
    checkList.append(('oxbs_report.R', passed))

    ## Reference FASTA
    sys.stdout.write('Check reference FASTA "%s"...     ' %(args.ref))
    if args.ref is None:
        tg= None
        passed= False
        sys.stdout.write('FAILED\n')
    elif os.path.isfile(args.ref):
        sys.stdout.write('OK\n')
        passed= True
    else:
        sys.stdout.write('FAILED\n')
        passed= False
    checkList.append(('Ref. FASTA', passed))
    
    ## Reference TXT
    sys.stdout.write('Check reference TXT "%s"...     ' %(args.ref))
    if args.ref is None:
        tg= None
        passed= False
        sys.stdout.write('FAILED\n')
    else:
        txt= os.path.splitext(args.ref)[0] + '.txt'
        if os.path.isfile(txt):
            sys.stdout.write('OK\n')
            passed= True
        else:
            sys.stdout.write('FAILED\n')
            passed= False
    checkList.append(('Ref. TXT', passed))
    
    ## List of BED positions:
    if args.listpos is not None:
        sys.stdout.write('Check bed file of positions "%s"...     ' %(args.listpos))
        if os.path.isfile(args.listpos):
            sys.stdout.write('OK\n')
            passed= True
        else:
            sys.stdout.write('FAILED\n')
            passed= False
    checkList.append(('Bed file of positions', passed))
    
    ## Check bowtie2 indexes
    ##refbt2= os.path.split(args.ref)[0] + 'Bisulfite_Genome'
    ## ...
    
    sys.stdout.write('Check prefix "%s"...     ' %(args.prefix))
    try:
        passed= validate_prefix(args.prefix)
    except ValidArgException:
        passed= False
    if passed:
        sys.stdout.write('OK\n')
    else:
        sys.stdout.write('FAILED\n')
    checkList.append(('prefix', passed))
        
    print('')
    print('-'*80)
    for x in checkList:
        if not x[1]:
            print(x[0] + '    ' + 'FAILED')
    print('')
    return(checkList)

def validate_listpos(arg):
    if arg is not None:
        if os.path.isfile(arg):
            return(True)
        else:
            raise ValidArgException('\n\nInput file --listpos "%s" not found\n' %(arg))
            return(False)
    else:
        return(True)

def validate_prefix(arg):
    if os.path.split(arg)[0] != '':
        raise ValidArgException('\n\nOutput prefix must not contain a directory path (use --outdir instead). Got "%s"\n' %(arg))
        return(False)
    return(True)
    
def validate_input(arg):
    """Check input arg is valid
    Return: 'bam', 'sam', 'raw' for bam, sam or raw read files.
    """
    if not type(arg) == list:
        raise ValidArgException('Input "%s" must be a list. Got %s' %(arg, type(arg)))
    
    if len(arg) != len(set(arg)):
        raise ValidArgException('\n\nDuplicate files found in input list %s\n' %(arg))
    
    bnames= [os.path.split(x)[1] for x in arg]
    bnames= [re.sub('\.gz$', '', x) for x in bnames]
    if len(bnames) == 2 and len(set(bnames)) == 1:
        raise ValidArgException('\n\nPaired fastq files must have different, unzipped names even if they are in different directories.\nGot %s\n' %(arg))
    
    for x in arg:
        if not os.path.isfile(x):
            raise ValidArgException('\n\nFile "%s" not found\n' %(x))
            
    if len(arg) == 2:
        return('raw')
    elif len(arg) == 1:
        ext= os.path.splitext(arg[0])[1]
        if ext in ['.sam', '.bam']:
            return(ext.strip('.'))
        else:
            return('raw')
    else:
        raise ValidArgException('\n\n1 or 2 item must be in input "%s". Got %s\n' %(arg, len(arg)))

def validate_trim_galore_opts(arg):
    """Check opts passed to trim_galore 
    """
    opts= shlex.split(arg)
    if '-o' in opts or '--output_dir' in opts:
        raise ValidArgException('\n\nOption -o/--output_dir should not be passed to trim_galore. Got "%s"\n' %(arg))
        return(False)
    return(True)

def validate_bismark_opts(arg):
    """Check opts passed to bismark 
    """
    opts= shlex.split(arg)
    if '-o' in opts or '--output_dir' in opts:
        raise ValidArgException('\n\nOption -o/--output_dir should not be passed to bismark. Got "%s"\n' %(arg))
        return(False)
    return(True)

def validate_mpileup_opts(arg):
    """Check opts passed to samtools mpileup
    """
    opts= shlex.split(arg)
    
    for x in opts:
        if x.startswith('-d'):
            raise ValidArgException('\n\nOption -d should not be passed to samtools mpileup. Got "%s"\n' %(arg))
            return(False)
    for x in opts:
        if x.startswith('-B'):
            raise ValidArgException('\n\nOption -B should not be passed to samtools mpileup. Got "%s"\n' %(arg))
            return(False)
    for x in opts:
        if x.startswith('-Q'):
            raise ValidArgException('\n\nOption -Q should not be passed to samtools mpileup. Got "%s"\n' %(arg))
            return(False)
    for x in opts:
        if x.startswith('-f'):
            raise ValidArgException('\n\nOption -f should not be passed to samtools mpileup. Got "%s"\n' %(arg))
            return(False)
    for x in opts:
        if x.startswith('-l'):
            raise ValidArgException('\n\nOption -l should not be passed to samtools mpileup. Got "%s"\n' %(arg))
            return(False)

    return(True)


## TODO: validate_reference_fasta()
