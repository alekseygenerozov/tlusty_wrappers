#!/usr/bin/env python

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
def setup(log_qgrav, log_teff, log_dmtot, lte, ltg, model, copy=True):
        #files to which we will write input for tlusty
        f5=open('fort.5', 'w')
        f1=open('fort.1', 'w')
        #Tail of tlusty input file assumed to already be present in the current directory
        tailf=open('tail', 'r')

        f1.write('1\n')

        qgrav=10**log_qgrav
        teff=10**log_teff
        dmtot=10**log_dmtot
        out='{0} {1:8.7E} {2:8.7E} {3:8.7E}'.format(0, teff, qgrav, dmtot)
        print out

        #Construct fort.5 input file for tlusty
        f5.write(out + '\n')
        #Second line of input
        f5.write(lte + '  ' + ltg + '\n')
        #Tail 
        tail=tailf.read()
        f5.write(tail) 

        #Format strings with info on teff, dmtot, qgrav for creating folder
        t=str(int(np.around(log_teff*10)))
        m=str(int(np.around(log_dmtot*10)))
        q=str(int(np.around(log_qgrav*10)))
        #Create a folder in which to store the tlusty output in 
        loc='t' + t + 'm' + m + 'q' + q
        bash_command('mkdir -p ' + loc)

        #If the location of the input model is not blank then copy model atmosphere to the current directory
        if model:
            bash_command('cp ' + model + '/fort.7 ' + './fort.8')
            if copy:
                bash_command('cp ' + model + '/tmp.flag ' + './')
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

##Run tlusty based on parameters found at the location of model
def tlusty_runs_model(model, nlte=False, copy=True):
    params=ascii.read(model + '/fort.5', data_start=0, data_end=1)
    #print params
    log_teff=np.log10(params[0][1])
    log_qgrav=np.log10(params[0][2])
    log_dmtot=np.log10(params[0][3])
    print log_teff,log_qgrav,log_dmtot

    lte='T'
    if(nlte):
        lte='F'
    ltg='F'

    outdir=setup(log_qgrav, log_teff, log_dmtot, lte, ltg, model, copy)
    run()
    maxchange=reltot()
    #Move tlusty output files to the appropriate directory
    clean(outdir,maxchange,nlte)



##Construct tlusty model based on info in myfile
def tlusty_runs_file(myfile='params.in', nlte='false', copy=True):
    params=ascii.read(myfile)
    ncols=len(params.columns)

    lte='T'
    if(nlte):
        lte='F'
    ltg='T'
    model=''

    #for each set of parameters in our input file
    for i in range(0, len(params)): 
        f=open('fort.5', 'w')
        tailf=open('tail', 'r')
        log_qgrav=params[i][0]
        log_teff=params[i][1]
        log_dmtot=params[i][2]

        if ncols>3:
            model=params[i][3]
        if model:
            ltg='F'

        print lte, ltg, model

        #set up input files, then run tlusty, finally check for nominal convergence and move all output files to 
        outdir=setup(log_qgrav, log_teff, log_dmtot, lte, ltg, model, copy)
        run()
        maxchange=reltot()
        #Move tlusty output files to the appropriate directory
        clean(outdir,maxchange,nlte)



##driver; parse user input
def main():
    parser=argparse.ArgumentParser(
        description='Wrapper for running TLUSTY')
    parser.add_argument('-n','--nlte',
        help='Switch on nlte',
        action='store_true')
    parser.add_argument('-f', '--file', 
        help='File containing input parameters; default is params.in', 
        default='params.in')
    parser.add_argument('-m', '--model',
        help='Set parameters (qgrav, teff, dmtot) to those at location' +
        'and use model atmosphere at this location; higher precedence compared to file')
    parser.add_argument('-nc', '--nocopy',
        help='This flag turns off the default behavior of copying tmp.flag from model location',
        action='store_true')
    args=parser.parse_args()

    myfile=args.file
    nlte=args.nlte
    model=args.model
    copy=not args.nocopy


    if  model:
        tlusty_runs_model(model, nlte, copy)
    else:
        tlusty_runs_file(myfile, nlte, copy)


if __name__ == '__main__':
    main()