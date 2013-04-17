import argparse
import numpy as np

import subprocess

import astropy.table as at
from astropy.io import ascii

##Run a command from the bash shell
def bash_command(cmd):
     process=subprocess.Popen(['/bin/bash', '-c', cmd],  stdout=subprocess.PIPE)
     process.wait()


##Setup input file for tlusty
def setup(log_qgrav, log_teff, log_dmtot, lte, ltg, inp_model):
        #files to which we will write input for tlusty
        f5=open('fort.5', 'w')
        f1=open('fort.1', 'w')
        #Tail of tlusty input file assumed to already be present in the current directory
        tailf=open('tail', 'r')

        f1.write('1\n')

        qgrav=10**log_qgrav
        teff=10**log_teff
        dmtot=10**log_dmtot

        #Construct fort.5 input file for tlusty
        f5.write(str(0) + ' ' + str(teff) + ' ' + str(qgrav) + ' ' +
            str(dmtot) + '\n')
        #Second line of input
        f5.write(str(lte) + '  ' + str(ltg) + '\n')
        #Tail 
        tail=tailf.read()
        f5.write(tail) 

        #Create a folder in which to store the tlusty output in 
        loc='t' + str(log_teff*10) + 'm' + str(log_dmtot*10) + 'q' + str(log_qgrav*10)
        bash_command("mkdir " + loc)

        #If the location of the input model is not blank then copy model atmosphere to the current directory
        if inp_model:
            bash_command('cp ' + inp_model + '/fort.7 ' './fort.8')


##Run TLUSTY from the bash shell
def run():
    cmd="./t202 <fort.5 >fort.6"
    bash_command(cmd)


##Construct tlusty model from scratch
def tlusty_runs(nlte=False):
    #read in list of parameters and extract the input parameters
    params=ascii.read('params.in', fill_values='')
    #number of columns for our table
    ncols=len(params.columns)

    lte='T'
    if(nlte):
        lte='F'
    ltg='T'
    inp_model=''

    #for each set of parameters in our input file
    for i in range(0, len(params)): 
        f=open('fort.5', 'w')
        tailf=open('tail', 'r')
        log_qgrav=params[i][0]
        log_teff=params[i][1]
        log_dmtot=params[i][2]

        if ncols>3:
            inp_model=params[i][3]
        if inp_model != '':
            ltg='F'

        print lte, ltg, inp_model

        setup(log_qgrav, log_teff, log_dmtot, lte, ltg, inp_model)
        run()


##driver; parse user input
def main():
    parser=argparse.ArgumentParser(
        description='Wrapper for running TLUSTY')
    parser.add_argument('--nlte',
        help='Switch on nlte',
        action='store_true')
    args=parser.parse_args()
    nlte=args.nlte

    tlusty_runs(nlte)


if __name__ == '__main__':
    main()