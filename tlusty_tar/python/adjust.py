#!/usr/bin/env python

import tlusty_runs as tr
import argparse
import re
import numpy as np
import subprocess
import astropy.table as at
from astropy.io import ascii
from types import *

##Adjusts parameter p in 
def adjust_file(p, end, file='tmp.flag'):
    #Open parameter file, read then close
    f=open(file, 'r')
    old=f.read()
    f.close()

    #Substitute new value of parameter into file
    s1=p+"=[0-9.\-]*"
    pat1=re.compile(s1)
    s2=p+"="+str(end)
    pat2=re.compile(s2)
    print s2
    new=pat1.sub(s2, old)
    #Overwrite old file so that the parameter is adjusted like we wanted
    f=open(file, 'w')
    f.write(new)


##Adjusts parameter p in tmp.flag until we have reached target
def adjust(p, start, target, delta, model, nlte=False):
    there=np.allclose([start], [target])
    #if we are close to our target adjust delta accordingly
    if there:
        return
    else:
        if abs(target-start)<abs(delta):
            delta=target-start

        #End point of the current step
        end=start+delta
        adjust_file(p, end)
        if model:
            tr.tlusty_runs_model(model, nlte=nlte)
        else:
            tr.tlusty_runs_file(nlte=nlte)

        start=end
        adjust(p, start, target, delta, model, nlte)







##driver; parse user input
def main():
    parser=argparse.ArgumentParser(description='Wrapper to gradually adjust parameters')
    parser.add_argument('param')
    parser.add_argument('start')
    parser.add_argument('target')
    parser.add_argument('delta')
    parser.add_argument('-m','--model',
        help='Model atmosphere location')
    parser.add_argument('-fl', '--floating', 
        help='Flag to turn on floating point adjustment to parameters',
        action='store_true')
    parser.add_argument('-n','--nlte',
        help='Switch on nlte',
        action='store_true')
    #parsing the input arguments
    args=parser.parse_args()
    p=args.param
    floating=args.floating
    model=args.model
    nlte=args.nlte


    #if the user has specified a floating flag, set start and delta parameters to be floats
    if floating:
        delta=float(args.delta)
        start=float(args.start)
        target=float(args.target)
    else:
        delta=int(args.delta)
        start=int(args.start)
        target=int(args.target)    

    adjust(p, start, target, delta, model, nlte)



if __name__ == '__main__':
    main()