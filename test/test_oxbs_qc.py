#!/usr/bin/env py.test


import sys
import os
import shutil
import gzip
import subprocess

sys.path.insert(0, '../') ## To import oxbs_qc

import oxbs_qc
from oxbs_qc import oxbs_qc_func
from oxbs_qc import mpileup2methylation
from oxbs_qc import validate_args
from oxbs_qc import cleanReadNames

testDir= 'test_out'
if not os.path.exists(testDir):
    os.makedirs(testDir)

def test_get_wdir():
    passed= False
    try:
        wdir= oxbs_qc_func.get_wdir('')
    except oxbs_qc_func.GetWdirException:
            passed= True
    assert passed
    
    passed= False
    try:
        wdir= oxbs_qc_func.get_wdir(['foo'])
    except oxbs_qc_func.GetWdirException:
            passed= True
    assert passed

    wdir= oxbs_qc_func.get_wdir('.')
    assert wdir == '.'
    wdir= oxbs_qc_func.get_wdir('test_out/test_wdir')
    assert wdir == 'test_out/test_wdir'
    assert os.path.exists(wdir) == True
    os.rmdir(wdir) ## _MEMO_: os.rmdir removes dir only if it is empty

def test_task_fastqc():
    try:
        shutil.rmtree('test_out/fastqc')
    except OSError:
        pass    
    os.mkdir('test_out/fastqc')

    xout= oxbs_qc_func.task_fastqc(['test_data/mjb042_oxBS_R1.fastq.gz', 'test_data/mjb042_oxBS_R2.fastq.gz'], 'test_out/fastqc')
    assert xout['returncode'] == 0
    assert os.path.exists('test_out/fastqc/mjb042_oxBS_R1_fastqc.zip')
    assert os.path.exists('test_out/fastqc/mjb042_oxBS_R2_fastqc.zip')
    
    passed= False
    try:
        ## Input not list
        xout= oxbs_qc_func.task_fastqc('test_data/mjb042_oxBS_R1.fastq.gz', 'test_out/fastqc')
    except oxbs_qc_func.FastqcException:
        passed= True
    assert passed

    passed= False
    try:
        ## Outdir does not exist
        xout= oxbs_qc_func.task_fastqc(['test_data/mjb042_oxBS_R1.fastq.gz'], '/bla/bla')
    except oxbs_qc_func.FastqcException:
        passed= True
    assert passed
    
    
def test_task_trim_galore_SE():
    passed= False
    try:
        tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.fastq.gz'], outFastq= ['bla/fq1.fq.gz'], tmpdir= 'test_out', opts= '-o test_out/trim_galore')
    except oxbs_qc_func.TrimGaloreException:
        "outFastq has path attached"
        passed= True
    assert passed

def test_task_trim_galore_SE_1():
    tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.fastq.gz'], outFastq= ['fq1.fq.gz'], tmpdir= 'test_out', opts= '-o test_out/trim_galore_1')
    print(tg)
    assert tg['fastq'][0] == 'test_out/trim_galore_1/fq1.fq.gz'
    assert os.path.exists(tg['fastq'][0])
    assert os.path.exists('test_out/trim_galore_1/mjb042_oxBS_R1.fastq.gz_trimming_report.txt')

def test_task_trim_galore_SE_2():
    tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.txt'], outFastq= ['fq1.txt'], tmpdir= 'test_out', opts= '-o test_out/trim_galore_2')
    print(tg)
    assert tg['fastq'][0] == 'test_out/trim_galore_2/fq1.txt'
    assert os.path.exists(tg['fastq'][0])

def test_task_trim_galore_SE_3():
    tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.txt.gz'], outFastq= ['fq1.txt'], tmpdir= 'test_out', opts= '-o test_out/trim_galore_3 --gzip')
    print(tg)
    assert tg['fastq'][0] == 'test_out/trim_galore_3/fq1.txt.gz'
    assert os.path.exists(tg['fastq'][0])

def test_task_trim_galore_SE_4():
    tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.txt.gz'], outFastq= ['fq1.txt'], tmpdir= 'test_out', opts= '-o test_out/trim_galore_4 --dont_gzip')
    print(tg)
    assert tg['fastq'][0] == 'test_out/trim_galore_4/fq1.txt'
    assert os.path.exists(tg['fastq'][0])

def test_task_trim_galore_PE_exceptions():
    ## 2 inputs 1 output
    passed= False
    try:
        tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.fastq.gz', 'test_data/mjb042_oxBS_R2.fastq.gz'], outFastq= ['fq1.fq.gz'], tmpdir= 'test_out', opts= '')
    except oxbs_qc_func.TrimGaloreException:
        passed= True
    assert passed

    ## 1 inputs 2 output
    passed= False
    try:
        tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.fastq.gz'], outFastq= ['fq1.fq.gz', 'fq2.fq.gz'], tmpdir= 'test_out', opts= '')
    except oxbs_qc_func.TrimGaloreException:
        passed= True
    assert passed

    ## same input twice
    passed= False
    try:
        tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.fastq.gz', 'test_data/mjb042_oxBS_R1.fastq.gz'], outFastq= ['fq1.fq.gz', 'fq2.fq.gz'], tmpdir= 'test_out', opts= '')
    except oxbs_qc_func.TrimGaloreException:
        passed= True
    assert passed

    ## same Output twice
    passed= False
    try:
        tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.fastq.gz', 'test_data/mjb042_oxBS_R2.fastq.gz'], outFastq= ['fq1.fq.gz', 'fq1.fq.gz'], tmpdir= 'test_out', opts= '')
    except oxbs_qc_func.TrimGaloreException:
        passed= True
    assert passed

def test_task_trim_galore_PE_1():
    tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.fastq.gz', 'test_data/mjb042_oxBS_R2.fastq.gz'], outFastq= ['fq1.fq.gz', 'fq2.fq.gz'], tmpdir= 'test_out', opts= '-o test_out/trim_galore_pe1')
    print(tg)
    assert tg['fastq'][0] == 'test_out/trim_galore_pe1/fq1.fq.gz' ## Files are named as you expect
    assert tg['fastq'][1] == 'test_out/trim_galore_pe1/fq2.fq.gz'
    assert os.path.exists(tg['fastq'][0]) ## Files do exist
    assert os.path.exists(tg['fastq'][1])

def test_task_trim_galore_PE_2():
    tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.txt.gz', 'test_data/mjb042_oxBS_R2.txt.gz'], outFastq= ['fq1.fq.gz', 'fq2.fq.gz'], tmpdir= 'test_out', opts= '-o test_out/trim_galore_pe2')
    print(tg)
    assert tg['fastq'][0] == 'test_out/trim_galore_pe2/fq1.fq.gz' ## Files are named as you expect
    assert tg['fastq'][1] == 'test_out/trim_galore_pe2/fq2.fq.gz'
    for x in tg['fastq']:
        assert os.path.exists(x) ## files do exist
    for x in tg['report']:
        assert os.path.exists(x) ## Reports do exist
    assert len(os.listdir('test_out/trim_galore_pe2/')) == 4

def test_task_trim_galore_PE_3():
    tg= oxbs_qc_func.task_trim_galore(inFastq= ['test_data/mjb042_oxBS_R1.txt', 'test_data/mjb042_oxBS_R2.txt'], outFastq= ['fq1.txt', 'fq2.txt'], tmpdir= 'test_out', opts= '-o test_out/trim_galore_pe3')
    print(tg)
    assert tg['fastq'][0] == 'test_out/trim_galore_pe3/fq1.txt' ## Files are named as you expect
    assert tg['fastq'][1] == 'test_out/trim_galore_pe3/fq2.txt'
    for x in tg['fastq']:
        assert os.path.exists(x) ## files do exist
    for x in tg['report']:
        assert os.path.exists(x) ## Reports do exist
    assert len(os.listdir('test_out/trim_galore_pe3/')) == 4

def test_task_shorten_fastq_exception():
    try:
        shutil.rmtree('test_out/shorten_fastq')
    except OSError:
        pass    
    os.mkdir('test_out/shorten_fastq')

    ## Invalid input
    passed= False
    try:
        oxbs_qc_func.task_shorten_fastq('/non/sense.fq', 'test_out/shorten_fastq/test.fq', 10)
    except oxbs_qc_func.ShortenFastqException:
        passed= True
    assert passed
    
    ## Invalid output
    passed= False
    try:
        oxbs_qc_func.task_shorten_fastq('test_data/mjb042_oxBS_R1.fastq.gz', '/non/sense/shorten_fastq/test.fq', 10)
    except oxbs_qc_func.ShortenFastqException:
        passed= True
    assert passed

    ## Invalid length
    passed= False
    try:
        oxbs_qc_func.task_shorten_fastq('test_data/mjb042_oxBS_R1.fastq.gz', 'test_out/shorten_fastq/test.fq', -10)
    except oxbs_qc_func.ShortenFastqException:
        passed= True
    assert passed

    ## All ok
    passed= False
    try:
        oxbs_qc_func.task_shorten_fastq('test_data/mjb042_oxBS_R1.fastq.gz', 'test_out/shorten_fastq/test.fq', 10)
        passed= True
    except oxbs_qc_func.ShortenFastqException:
        passed= False
    assert passed

def test_task_shorten_fastq_run():
    try:
        shutil.rmtree('test_out/shorten_fastq')
    except OSError:
        pass    
    os.mkdir('test_out/shorten_fastq')

    xout= oxbs_qc_func.task_shorten_fastq('test_data/mjb042_oxBS_R1.fastq.gz', 'test_out/shorten_fastq/test.fq.gz', 10)
    assert xout['returncode'] == 0
    assert os.path.isfile(xout['out'])
    
    orifq= gzip.open('test_data/mjb042_oxBS_R1.fastq.gz').readlines()
    shortfq= gzip.open(xout['out']).readlines()
    assert len(orifq) == len(shortfq) ## Do not loose any line.
    n= -1
    for line in shortfq:
        "Reads and quals == short len"
        if n % 2 == 0:
            assert len(line.strip()) == 10
        n += 1 
#def test_trim_fastq_exception():
#    passed= False
#    try:
#        oxbs_qc_func.trim_fastq('test_data/mjb042_oxBS_R1.txt', -1, 'mjb042_oxBS_R1.fq') ### Invalid -1
#    except oxbs_qc_func.TrimFastqException:
#        passed= True
#    assert passed
#
#    passed= False
#    try:
#        x= oxbs_qc_func.trim_fastq('test_data/mjb042_oxBS_R1.txt.gz', 1, 'test_data/mjb042_oxBS_R1.txt.gz') ## input == output
#    except oxbs_qc_func.TrimFastqException:
#        passed= True
#    assert passed
#
#    passed= False
#    try:
#        x= oxbs_qc_func.trim_fastq('test_data/mjb042_oxBS_R1.txt.gz', 1, '/nonsense/mjb042_oxBS_R1.txt.gz') ## output dir does not exists
#    except oxbs_qc_func.TrimFastqException:
#        passed= True
#    assert passed
#
#def test_trim_fastq_1():
#    try:
#        os.makedirs('test_out')
#    except OSError:
#        pass
#    xout= oxbs_qc_func.trim_fastq('test_data/mjb042_oxBS_R1.txt.gz', 20, 'test_out/mjb042_oxBS_R1.fq')
#    assert xout
#    assert os.path.exists('test_out/mjb042_oxBS_R1.fq.gz') ## Note .gz included!
#
#    ## Same number of lines in input and output
#    tfq= gzip.open('test_out/mjb042_oxBS_R1.fq.gz').readlines()
#    fq= gzip.open('test_data/mjb042_oxBS_R1.txt.gz').readlines()
#    assert len(tfq) == len(fq)
#
#    ## Lenght of read AND quality == 20
#    read_quality_qlen= [len(tfq[x].strip()) for x in range(1, 100, 2)]
#    assert set(read_quality_qlen) == set([20])

def test_task_bismark_exceptions():
    passed= False
    try:
        oxbs_qc_func.task_bismark_aln(inFastq= ['NA'], ref= '../control_reference/bsseq_synthetic4', outSam= 'mjb042_oxBS_R1.sam', opts= '') ## No input fastq
    except oxbs_qc_func.BismarkException:
        passed= True
    assert passed

    passed= False
    try:
        oxbs_qc_func.task_bismark_aln(inFastq= ['test_data/mjb042_oxBS_R1.txt.gz'], ref= '../control_reference/bsseq_synthetic4', outSam= 'mjb042_oxBS_R1', opts= '') ## No exts on sam
    except oxbs_qc_func.BismarkException:
        passed= True
    assert passed

    passed= False
    try:
        x= oxbs_qc_func.task_bismark_aln(inFastq= ['test_data/mjb042_oxBS_R1.txt.gz'], ref= '../control_reference/bsseq_synthetic4/NA', outSam= 'mjb042_oxBS_R1.sam', opts= '') ## No ref
        print(x)
    except oxbs_qc_func.BismarkException:
        passed= True
    assert passed

def test_task_bismark_aln_1():
    x= oxbs_qc_func.task_bismark_aln(inFastq= ['test_data/mjb042_oxBS_R1.txt.gz'], ref= '../control_reference/bsseq_synthetic4', tmpdir= 'test_out', outSam= 'mjb042_oxBS_R1.sam', opts= '-o test_out')
    print(x)
    assert x['sam'] == 'test_out/mjb042_oxBS_R1.sam' ## Outut sam is what you expect
    assert os.path.exists(x['sam']) ## Sam exists
    assert len(x['all_out']) > 1 ## There are more files in output dir
    
def test_task_bismark_aln_2():
    x= oxbs_qc_func.task_bismark_aln(inFastq= ['test_data/mjb042_oxBS_R1.txt.gz', 'test_data/mjb042_oxBS_R2.txt.gz'], ref= '../control_reference/bsseq_synthetic4', tmpdir= 'test_out', outSam= 'mjb042.sam', opts= '-o test_out')
    print(x)
    assert x['sam'] == 'test_out/mjb042.sam' ## Outut sam is what you expect
    assert os.path.exists(x['sam']) ## Sam exists
    assert len(x['all_out']) > 1 ## There are more files in output dir

#def test_cleanSamReadNames():
#    try:
#        os.mkdir('test_out/cleanSamReadNames')
#    except OSError:
#        pass    
#    shutil.copy('test_data/cc035.pe.sam', 'test_out/cleanSamReadNames/cc035.pe.sam')
#    orifile= open('test_out/cleanSamReadNames/cc035.pe.sam').readlines()
#    x= oxbs_qc_func.cleanSamReadNames('test_out/cleanSamReadNames/cc035.pe.sam')
#    newfile= open('test_out/cleanSamReadNames/cc035.pe.sam').readlines()
#    assert x
#    assert len(orifile) == len(newfile)
#    suffix_1= [x for x in newfile if '/1' in x]
#    suffix_2= [x for x in newfile if '/2' in x]
#    print(suffix_1)
#    print(suffix_2)
#    assert suffix_1 == []
#    assert suffix_2 == []

def test_cleanReadNames():
    try:
        os.mkdir('test_out/cleanReadNames')
    except OSError:
        pass    
    line_without_trailing= """M01567:30:000000000-A4FG4:1:1101:11780:2230_1:N:0:2	0	SY009_3x3_hmc_Q68	1	255	70M	*	0	0	TTACTATCTATTATTCAATTATTTAAATATGATTAAATAACATTAATATATTATCGATTAAATAAGAATT	ABABAFFFFFFFGGGGGGGGDGGHFHFHHHGHHHHFHGGHHHHGEHGBGGFFHHFFHHHGHGHHHHGHHH	NM:i:8	XX:Z:1C2C5CC2C3CC9C40	XM:Z:.h.Hh..H..hh..hH..hh.........z..........H.............Z...............	XR:Z:CT	XG:Z:CT"""
    outline= cleanReadNames.cleanSamReadNames(line_without_trailing)
    assert outline == line_without_trailing

    line_without_trailing_with_newline= """M01567:30:000000000-A4FG4:1:1101:11780:2230_1:N:0:2	0	SY009_3x3_hmc_Q68	1	255	70M	*	0	0	TTACTATCTATTATTCAATTATTTAAATATGATTAAATAACATTAATATATTATCGATTAAATAAGAATT	ABABAFFFFFFFGGGGGGGGDGGHFHFHHHGHHHHFHGGHHHHGEHGBGGFFHHFFHHHGHGHHHHGHHH	NM:i:8	XX:Z:1C2C5CC2C3CC9C40	XM:Z:.h.Hh..H..hh..hH..hh.........z..........H.............Z...............	XR:Z:CT	XG:Z:CT\n"""
    outline= cleanReadNames.cleanSamReadNames(line_without_trailing_with_newline)
    assert outline == line_without_trailing_with_newline

    line_with_1= """M01567:30:000000000-A4FG4:1:1101:11780:2230_1:N:0:2/1	0	SY009_3x3_hmc_Q68	1	255	70M	*	0	0	TTACTATCTATTATTCAATTATTTAAATATGATTAAATAACATTAATATATTATCGATTAAATAAGAATT	ABABAFFFFFFFGGGGGGGGDGGHFHFHHHGHHHHFHGGHHHHGEHGBGGFFHHFFHHHGHGHHHHGHHH	NM:i:8	XX:Z:1C2C5CC2C3CC9C40	XM:Z:.h.Hh..H..hh..hH..hh.........z..........H.............Z...............	XR:Z:CT	XG:Z:CT\n"""
    outline= cleanReadNames.cleanSamReadNames(line_with_1)
    assert outline == line_without_trailing_with_newline ## Must equal the version w/o trailing

    line_with_2= """M01567:30:000000000-A4FG4:1:1101:11780:2230_1:N:0:2/2	0	SY009_3x3_hmc_Q68	1	255	70M	*	0	0	TTACTATCTATTATTCAATTATTTAAATATGATTAAATAACATTAATATATTATCGATTAAATAAGAATT	ABABAFFFFFFFGGGGGGGGDGGHFHFHHHGHHHHFHGGHHHHGEHGBGGFFHHFFHHHGHGHHHHGHHH	NM:i:8	XX:Z:1C2C5CC2C3CC9C40	XM:Z:.h.Hh..H..hh..hH..hh.........z..........H.............Z...............	XR:Z:CT	XG:Z:CT\n"""
    outline= cleanReadNames.cleanSamReadNames(line_with_2)
    assert outline == line_without_trailing_with_newline ## Must equal the version w/o trailing

    ## Trailers other than /1 and /2 remain:
    line_with_other_trailing= """M01567:30:000000000-A4FG4:1:1101:11780:2230_1:N:0:2/3	0	SY009_3x3_hmc_Q68	1	255	70M	*	0	0	TTACTATCTATTATTCAATTATTTAAATATGATTAAATAACATTAATATATTATCGATTAAATAAGAATT	ABABAFFFFFFFGGGGGGGGDGGHFHFHHHGHHHHFHGGHHHHGEHGBGGFFHHFFHHHGHGHHHHGHHH	NM:i:8	XX:Z:1C2C5CC2C3CC9C40	XM:Z:.h.Hh..H..hh..hH..hh.........z..........H.............Z...............	XR:Z:CT	XG:Z:CT\n"""
    outline= cleanReadNames.cleanSamReadNames(line_with_other_trailing)
    assert outline == line_with_other_trailing

    ## Header lines untouched:
    header= """@SQ     SN:SY007_1x2_hmc_Q38_indexed    LN:100\n"""
    outline= cleanReadNames.cleanSamReadNames(header)
    assert outline == header

        
    cmd= 'python ../oxbs_qc/cleanReadNames.py -h'
    p= subprocess.Popen(cmd, shell= True, stderr= subprocess.PIPE, stdout= subprocess.PIPE)    
    stdout, stderr= p.communicate()
    assert p.returncode == 0

    cmd= 'python ../oxbs_qc/cleanReadNames.py -i test_data/cc035.sam'
    p= subprocess.Popen(cmd, shell= True, stderr= subprocess.PIPE, stdout= subprocess.PIPE)    
    stdout, stderr= p.communicate()
    assert p.returncode == 0
    ## Num lines in original sam must be the same.
    ## NB: -1 to remove triling add by subprocess
    assert (len(stdout.split('\n')) - 1) == 2140 

    ## Pipe from samtools SAM:
    cmd= 'samtools view -S -h test_data/cc035.sam | python ../oxbs_qc/cleanReadNames.py -i - '
    p= subprocess.Popen(cmd, shell= True, stderr= subprocess.PIPE, stdout= subprocess.PIPE)    
    stdout, stderr= p.communicate()
    assert p.returncode == 0
    ## Num lines in original sam must be the same.
    ## NB: -1 to remove triling add by subprocess
    assert (len(stdout.split('\n')) - 1) == 2140 

    ## Pipe from samtools BAM:
    cmd= 'samtools view -h test_data/cc035.bam | python ../oxbs_qc/cleanReadNames.py -i - '
    p= subprocess.Popen(cmd, shell= True, stderr= subprocess.PIPE, stdout= subprocess.PIPE)    
    stdout, stderr= p.communicate()
    assert p.returncode == 0
    ## Num lines in original sam must be the same.
    ## NB: -1 to remove triling add by subprocess
    assert (len(stdout.split('\n')) - 1) == 2140 

    
def test_task_sam2bam_with_sam():
    try:
        shutil.rmtree('test_out/sam2bam_with_sam')
    except OSError:
        pass    
    os.mkdir('test_out/sam2bam_with_sam')
    shutil.copy('test_data/cc035.sam', 'test_out/sam2bam_with_sam/cc035.sam')
    p= oxbs_qc_func.task_sam2bam('test_out/sam2bam_with_sam/cc035.sam', rmSam= True, path= '')
    print(p['cmd'])
    assert p['returncode'] == 0
    assert os.path.exists('test_out/sam2bam_with_sam/cc035.bam')
    assert os.path.exists('test_out/sam2bam_with_sam/cc035.bam.bai')
    assert p['bam'] == 'test_out/sam2bam_with_sam/cc035.bam'


def test_task_sam2bam_with_bam():
    try:
        shutil.rmtree('test_out/sam2bam_with_bam')
    except OSError:
        pass    
    os.mkdir('test_out/sam2bam_with_bam')
    shutil.copy('test_data/cc035.bam', 'test_out/sam2bam_with_bam/cc035.bam')
    p= oxbs_qc_func.task_sam2bam('test_out/sam2bam_with_bam/cc035.bam', rmSam= True, path= '')
    print(p['cmd'])
    assert p['returncode'] == 0
    assert os.path.exists('test_out/sam2bam_with_bam/cc035.bam')
    assert os.path.exists('test_out/sam2bam_with_bam/cc035.bam.bai')
    assert p['bam'] == 'test_out/sam2bam_with_bam/cc035.bam'

def test_task_markDuplicates():
    try:
        shutil.rmtree('test_out/markDuplicates')
    except OSError:
        pass    
    os.mkdir('test_out/markDuplicates')
    shutil.copy('test_data/cc035.sorted.bam', 'test_out/markDuplicates/cc035.sorted.bam')
    passed= False
    try:
        ## input not found/not valid
        x= oxbs_qc_func.task_markDuplicates('test_data/cc035.bam.x')
    except oxbs_qc_func.MarkDuplicatesException:
        passed= True
        pass
    assert passed
    
    passed= False
    try:
        ## .jar not found
        x= oxbs_qc_func.task_markDuplicates('test_data/cc035.bam', path= '.')
    except oxbs_qc_func.MarkDuplicatesException:
        passed= True
        pass
    assert passed

    p= oxbs_qc_func.task_markDuplicates('test_out/markDuplicates/cc035.sorted.bam')
    print(p)    
    assert os.path.isfile('test_out/markDuplicates/cc035.sorted.mdup.bai')
    assert os.path.isfile(p['bam'])
    assert os.path.isfile(p['metrics_file'])
   
def test_task_clipOverlap():
    """MEMO: Need a real paired end file.
    """
    try:
        shutil.rmtree('test_out/clipOverlap')
    except OSError:
        pass    
    os.mkdir('test_out/clipOverlap')
    p= oxbs_qc_func.task_clipOverlap(inBam= 'test_data/cc035.sorted.bam', clippedBam= 'test_out/clipOverlap/cc035.sorted.bam')
    print(p['cmd'])
    assert p['returncode'] == 0
    assert os.path.exists('test_out/clipOverlap/cc035.sorted.bam')
    assert os.path.exists('test_out/clipOverlap/cc035.sorted.bam.bai')
    assert os.path.exists('test_out/clipOverlap/cc035.sorted.bam.clipStats')    
    
def test_cleanMpileup():
    """Check meyhylatio2mpileup correclty removes indels and chars including and
    following ^.
    """
    bases= '^A...+11CCCCCCCCCCC.,,.,.,..-5TTTTT.,,,.,....^k.'
    cleanString= mpileup2methylation.cleanCallString(bases)
    assert cleanString ==  '....,,.,.,...,,,.,.....'

def test_task_callMethylation():
    try:
        shutil.rmtree('test_out/callMethylation')
    except OSError:
        pass    
    os.mkdir('test_out/callMethylation')
    p= oxbs_qc_func.task_callMethylation(inBam= 'test_data/cc035.sorted.bam', outBdg= 'test_out/callMethylation/cc035.sorted.bedGraph', opts= '-f ../control_reference/bsseq_synthetic4/bsseq_synthetic4.fa')
    assert p['returncode'] == 0
    print(p['stderr'])
    assert os.path.exists('test_out/callMethylation/cc035.sorted.bedGraph')
    bdg= open('test_out/callMethylation/cc035.sorted.bedGraph').readlines()
    assert len(bdg) == 220

def test_task_oxbs_report():
    try:
        shutil.rmtree('test_out/oxbs_report')
    except OSError:
        pass    
    os.mkdir('test_out/oxbs_report')
    p= oxbs_qc_func.task_oxbs_report(inBdg= 'test_data/cc035.bedGraph', ref= '../control_reference/bsseq_synthetic4/bsseq_synthetic4.txt', output_dir= 'test_out/oxbs_report', prefix= 'pytest')
    assert p['returncode'] == 0

    ## Files exist and are not empty
    outfiles= sorted(os.listdir('test_out/oxbs_report'))
    assert outfiles == ['pytest.conversion.pdf', 'pytest.coverage.pdf', 'pytest.oxqc.txt', 'pytest.oxqc_avg.txt']
    assert set([os.stat(os.path.join('test_out/oxbs_report', x)).st_size > 0 for x in outfiles]) == set([True])
    
    ## You have the expected number of lines
    flen= open('test_out/oxbs_report/pytest.oxqc.txt').readlines()
    assert len(flen) == 224

# -----------------------------------------------------------------------------
# Testing validate args
# -----------------------------------------------------------------------------
def test_validate_input():
    x= validate_args.validate_input(['test_data/cc035.pe.sam'])        
    assert x == 'sam'
    x= validate_args.validate_input(['test_data/cc035.sorted.bam'])
    assert x == 'bam'
    x= validate_args.validate_input(['test_data/mjb042_oxBS_R1.fq.gz'])
    assert x == 'raw'
    passed= False
    try:
        x= validate_args.validate_input(['test_data/NA'])
    except validate_args.ValidArgException:
        passed= True
    assert passed

def test_settings():
    """
    cd /Users/berald01/svn_checkout/oxbs-sequencing-qc/pipeline_kit/test
    ../oxbs_qc/oxbs_qc.py -i test_data/mjb042_oxBS_R1.fq.gz -r ../control_reference/bsseq_synthetic4/bsseq_synthetic4.fa -o test_out/oxbs_qc/
    """
    cmd= "../oxbs_qc/oxbs_qc.py -i test_data/mjb042_oxBS_R1.fq.gz -r ../control_reference/bsseq_synthetic4/bsseq_synthetic4.fa -o /usr --check"
    p= subprocess.Popen(cmd, shell= True, stderr= subprocess.PIPE, stdout= subprocess.PIPE)
    stdout, stderr= p.communicate()
    print(stderr)
    print(stdout)
    assert p.returncode == 0
    assert 'FAILED' in stdout
    
    
# ------------------------------------------------------------------------------
# Test workflow
# ------------------------------------------------------------------------------

"""
cd /Users/berald01/svn_checkout/oxbs-sequencing-qc/pipeline_kit/test
rm test_out/oxbs_qc/*
../oxbs_qc/oxbs_qc.py -i test_data/mjb042_oxBS_R1.fq.gz -r ../control_reference/bsseq_synthetic4/bsseq_synthetic4.fa -o test_out/oxbs_qc/ -V -p oxbs_qc --mpileup_opts ' -l bla'

rm test_out/oxbs_qc/*
../oxbs_qc/oxbs_qc.py -i test_data/mjb042_oxBS_R1.fq.gz test_data/mjb042_oxBS_R2.fq.gz -r ../control_reference/bsseq_synthetic4/bsseq_synthetic4.fa -o test_out/oxbs_qc/ -V -p oxbs_qc -Md

rm test_out/oxbs_qc/*
bsExpress -i test_data/mjb042_oxBS_R1.fq.gz test_data/mjb042_oxBS_R2.fq.gz -r ../control_reference/bsseq_synthetic4/bsseq_synthetic4.fa -o test_out/oxbs_qc/ -V -p oxbs_qc -Md
"""

#def test_