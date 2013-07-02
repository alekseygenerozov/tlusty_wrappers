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
    s1=p+'=[0-9.\-]*'
    s2=p+'='+'{0:.3g}'.format(end)
    pat1=re.compile(s1)

    #print old
    #print pat1.search(old)
    #If the parameter is present in the old file then substitute the new value 
    if pat1.search(old):
        pat2=re.compile(s2)
        new=pat1.sub(s2, old)
        f=open(file, 'w')
        f.write(new) 
    #If the parameter is not present then simply append it to the end of the optional parameters file
    f=open(file, 'a')
    f.write(s2+'\n')


## Adjusts parameter p in tmp.flag until we have reached target
def adjust(p, start, target, delta, model, nlte=False):
    there=np.allclose([start], [target])
    print model
    #if we are there then done
    if there:
        return
    #if close to target then adjust the step appropriately
    else:
        if abs(target-start)<abs(delta):
            delta=target-start

        #End point of the current step
        end=start+delta
        print end
        #Copy the original optional parameters file 
        if model:
            tr.bash_command('cp ' + model + '/tmp.flag ' + './')
        adjust_file(p, end)

        maxchange=1
        if model:
            out=tr.tlusty_runs_model(model, nlte=nlte, copy=False)
            maxchange=out[0]
            model=out[1]
            print maxchange
        else:
            maxchange=tr.tlusty_runs_file(nlte=nlte, copy=False)

        #If we have reached a model which doesn't converge then return 
        maxchange= np.log10(np.absolute(maxchange))
        if np.isnan(maxchange):
            maxchange=1000
        if maxchange>0:
            print 'Atmosphere does not appear to converge'
            return

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
    # parser.add_argument('-fl', '--floating', 
    #     help='Flag to turn on floating point adjustment to parameters',
    #     action='store_true')
    parser.add_argument('-n','--nlte',
        help='Switch on nlte',
        action='store_true')
    #parsing the input arguments
    args=parser.parse_args()
    p=args.param
    #floating=args.floating
    model=args.model
    nlte=args.nlte

    #if the user has specified a floating flag, set start and delta parameters to be floats
    delta=float(args.delta)
    start=float(args.start)
    target=float(args.target)
    # if floating:
    #     delta=float(args.delta)
    #     start=float(args.start)
    #     target=float(args.target)
    # else:
    #     delta=int(args.delta)
    #     start=int(args.start)
    #     target=int(args.target)    

    adjust(p, start, target, delta, model, nlte)



if __name__ == '__main__':
    main()