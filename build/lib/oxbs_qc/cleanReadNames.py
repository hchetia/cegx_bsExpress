#!/usr/bin/env python

import sys
import os
import argparse

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Remove the /1 or /2 suffix from paired read names.

    bowtie2 and possibly other aligners append /1 and /2 to the read name first
    and second in pair. This causes problems to tools like bamUtils which expect
    paired read names to be identical.

INPUT:
    Must be SAM file either as file or as stream. File name must end in .sam.

EXAMPLE
    ## Read SAM note -h option to samtools and -S
    samtools view -h test.bam | cleanBamReadNames.py -i -

    ## Use with bamUtils
    samtools view -h test.bam | cleanBamReadNames.py -i - | bam clipOverlap --in -.sam --out test.clip.bam --stats

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--input', '-i',
                   required= True,
                   help='''Input sam file to clean. Use - to read from stdin.

''')

def cleanSamReadNames(samline):
    """Removes trailing /1 /2 from read names to have pairs with the same name
    samline:
        Line from sam file or stream. Lines with leading @ printed out as they are.
    """
    if samline.startswith('@'):
        return(samline)
    else:
        samline= samline.split('\t')
        readName= samline[0]
        if readName[-2:] in ('/1', '/2'):
            readName= readName[0:-2]
            samline[0]= readName
        return('\t'.join(samline))

def main():
    args = parser.parse_args()
    if args.input == '-':
        samin= sys.stdin
    elif not os.path.isfile(args.input):
        sys.exit('File "%s" not found' %(args.input))
    elif not args.input.endswith('.sam'):
        sys.exit('Input file must have extension .sam. Got %s' %(args.input))
    else:
        samin= open(args.input)
    for samline in samin:
        outline= cleanSamReadNames(samline)
        sys.stdout.write(outline)
    samin.close()

if __name__ == '__main__':
    main()
    sys.exit()
