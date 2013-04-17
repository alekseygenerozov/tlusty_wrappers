import argparse
import numpy as np

import subprocess

import astropy.table as at
from astropy.io import ascii

from types import*

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
        loc='t' + str(int(log_teff*10)) + 'm' + str(int(log_dmtot*10)) + 'q' + str(int(log_qgrav*10))
        bash_command('mkdir -p ' + loc)

        #If the location of the input model is not blank then copy model atmosphere to the current directory
        if inp_model:
            bash_command('cp ' + inp_model + '/fort.7 ' './fort.8')
        #Return name of directory setup made   
        return loc



##Run TLUSTY from the bash shell
def run():
    cmd="./t202 <fort.5 >fort.6"
    bash_command(cmd)

##Gets nominal convergence from the convergence log file in tlusty
def reltot(file='fort.9'):
    dat=ascii.read(file, header_start=1)
    niter= dat['ITER'][-1]
    maxchange= dat['MAXIMUM'][-1]
    print niter, maxchange
   
    maxchange=float(maxchange)
    return maxchange

##Move all tlusty output files to the apropriate directory    
def clean(outdir,maxchange,nlte):
    maxchange= np.log10(maxchange)
    if np.isnan(maxchange):
        maxchange=1000

    #constructing destination path
    dest='./'
    if maxchange<0 and nlte:
        dest=dest + outdir + '/converged'
    elif maxchange<0:
        dest=dest + outdir + '/lteconv'
    else:
        dest=dest + outdir + '/nconv'
    print(dest)
    #delete any files in the destination folder (i.e. from a previous run of tlusty)
    bash_command('rm -f ' + dest + '/*')
    bash_command('mkdir -p ' + dest)
    #moving all tlusty i/o files to the destination
    bash_command('mv ' + 'fort* ' + dest)
    #move optional parameter file to destination
    bash_command('cp ' + 'tmp.flag ' + dest)


##Construct tlusty model 
def tlusty_runs(nlte=False):
    #read in list of parameters and extract the input parameters
    params=ascii.read('params.in')
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
        if inp_model:
            ltg='F'

        print lte, ltg, inp_model

        #set up input files, then run tlusty, finally check for nominal convergence and move all output files to 
        outdir=setup(log_qgrav, log_teff, log_dmtot, lte, ltg, inp_model)
        run()
        maxchange=reltot()
        #Move tlusty output files to the appropriate directory
        clean(outdir,maxchange,nlte)



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