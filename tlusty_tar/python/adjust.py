import tlusty_runs
import argparse
import re
import numpy as np
import subprocess
import astropy.table as at
from astropy.io import ascii
from types import *


def adjust(p, start, delta, file='tmp.flag'):
    #Open parameter file, read then close
    f=open(file, 'r')
    old=f.read()
    f.close()
    #Constructing pattern to find in file and its replacement
    #s0=p+"=[0-9.\-]*"  
    #Substitute new value of parameter into file
    s1=p+"="+str(start)
    pat1=re.compile(s1)
    s2=p+"="+str(start+delta)
    pat2=re.compile(s2)
    print s1,s2
    new=pat1.sub(s2, old)

    f=open(file, 'w')
    f.write(new)





##driver; parse user input
def main():
    parser=argparse.ArgumentParser(description='Wrapper to gradually adjust parameters')
    parser.add_argument('param')
    parser.add_argument('delta')
    parser.add_argument('start')
    parser.add_argument('-f', '--floating', 
        help='Flag to turn on floating point adjustment to parameters',
        action='store_true')
    #parsing the input arguments
    args=parser.parse_args()
    p=args.param
    floating=args.floating
    #if the user has specified a floating flag, set start and delta parameters to be floats
    if floating:
        delta=float(args.delta)
        start=float(args.start)
    #otherwise they should be integers    
    else:
        delta=int(args.delta)
        start=int(args.start)


    print p, delta

    adjust(p, start,delta)



if __name__ == '__main__':
    main()