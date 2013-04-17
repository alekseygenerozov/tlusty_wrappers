import argparse
import numpy as np

import subprocess

import astropy.table as at
from astropy.io import ascii


##Setup input file for tlusty
def setup(log_qgrav, log_teff, log_dmtot):
        f5=open('fort.5', 'w')
        f1=open('fort.1', 'w')
        tailf=open('tail', 'r')

        f1.write('1\n')


        qgrav=10**log_qgrav
        teff=10**log_teff
        dmtot=10**log_dmtot

        #Construct fort.5 input file for tlusty
        f5.write(str(0) + ' ' + str(teff) + ' ' + str(qgrav) + ' ' +
            str(dmtot) + '\n')
        #second line of input
        f5.write('T  T\n')
        #tail 
        tail=tailf.read()
        f5.write(tail) 


def run():
    # bash_command="nice -n 19 ./t202 <fort.5 >fort.6"
    bash_command="./t202 <fort.5 >fort.6"
    process=subprocess.Popen(['/bin/bash', '-c', bash_command])
    #process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE, shell=True)
    #report any problems that occur while reading the line
    # while True:
    #     line = process.stdout.readline()
    #     if not line:
    #         break
    #     print line
    #output = process.communicate()[0]


##Construct tlusty model from scratch
def tlusty_runs(nlte=False):
    #read in list of parameters and extract the input parameters
    params=ascii.read('params.in', fill_values='')
    #number of columns for our table
    ncols=len(params.columns)

    if(nlte):
        lte='F'
    else:
        lte='T'

    for i in range(0, len(params)): 
        f=open('fort.5', 'w')
        tailf=open('tail', 'r')
        log_qgrav=params[i][0]
        log_teff=params[i][1]
        log_dmtot=params[i][2]



        setup(log_qgrav, log_teff, log_dmtot)
        run()


##driver; parse user input
def main():
    parser=argparse.ArgumentParser(
        description='Wrapper for running TLUSTY')
    # parser.add_argument('--model',
    #     help='Look for models.in in the current dir',
    #     action='store_true')
    parser.add_argument('--nlte',
        help='Switch on nlte',
        action='store_true')
    args=parser.parse_args()
    nlte=args.nlte

    # if model:
    #     tlusty_runs_model(nlte)
    # else:
    #     tlusty_runs_param(nlte)
    tlusty_runs(nlte)


if __name__ == '__main__':
    main()