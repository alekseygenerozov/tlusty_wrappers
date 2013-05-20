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
     return process


##Setup input file for tlusty
def setup(log_qgrav, log_teff, log_dmtot, nlte, model, copy=True, tailname='tail'):
        lte='T'
        if(nlte):
            lte='F'
        ltg='T'
        if (model):
            ltg='F'

        #files to which we will write input for tlusty
        f5=open('fort.5', 'w')
        f1=open('fort.1', 'w')

        #Tail of tlusty input file assumed to already be present in the current directory
        print(tailname)
        tailf=open(tailname, 'r')

        f1.write('1\n')

        qgrav=10.**log_qgrav
        teff=10.**log_teff
        dmtot=10.**log_dmtot
   
        out='{0} {1:8.7e} {2:8.7e} {3:8.7e}'.format(0, teff, qgrav, dmtot)
        print out
        # fort.5 input file for tlusty
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
           
            process=bash_command('check ' + model + '/fort.7')

            if len(process.stdout.readlines())==0:
                print 'test'
                loc=''
                return loc
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
    dat=np.genfromtxt(file, skip_header=3)

    niter= dat[-1, 0]
    print niter
    dat_cut=dat[(dat[:,0]==niter)]


    try:
         change=dat_cut[:, 2:7]
         change=np.array(map(np.abs, change))
         maxchange=change.max()
         print maxchange
    except ValueError:
         maxchange=1000

    #Using numpy's is nan function to check if the maximum change is a nan
    if np.isnan(maxchange):
         maxchange=1000
    #Probably not good if the maximum is precisely zero
    if maxchange==0:
         return 1000

    print niter, maxchange
    return maxchange


##Move all tlusty output files to the apropriate directory    
def clean(outdir,maxchange,nlte):
    maxchange= np.log10(np.absolute(maxchange))

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

    return dest


##Run tlusty based on command line input parameters
def tlusty_runs_input(params, model, nlte=False, copy=True, tailname='tail'):
    log_qgrav=params[0]
    log_teff=params[1]
    log_dmtot=params[2]

    print log_teff,log_qgrav,log_dmtot


    outdir=setup(log_qgrav, log_teff, log_dmtot, nlte, model, copy, tailname)
    if not outdir:
        return
    run()
    maxchange=reltot()
    #Move tlusty output files to the appropriate directory
    clean(outdir,maxchange,nlte)
    #return maxchange



##Run tlusty based on parameters found at the location of model
def tlusty_runs_model(model, nlte=False, copy=True, tailname='tail'):
    params=ascii.read(model + '/fort.5', data_start=0, data_end=1)
    #print params
    log_teff=np.log10(params[0][1])
    log_qgrav=np.log10(params[0][2])
    log_dmtot=np.log10(params[0][3])
    print log_teff,log_qgrav,log_dmtot


    outdir=setup(log_qgrav, log_teff, log_dmtot, nlte, model, copy, tailname)
    if not outdir:
        return

    run()
    maxchange=reltot()
    #Move tlusty output files to the appropriate directory
    clean(outdir,maxchange,nlte)
    #return maxchange



##Construct tlusty model based on info in myfile
def tlusty_runs_file(myfile='params.in', nlte=False, copy=True, combo=False, tailname='tail'):
    params=ascii.read(myfile)
    ncols=len(params.columns)

    model=''

    #for each set of parameters in our input file
    for i in range(0, len(params)): 
        f=open('fort.5', 'w')
        tailf=open('tail', 'r')
        print(params[i])
        log_qgrav=params[i][0]
        log_teff=params[i][1]
        log_dmtot=params[i][2]

        if ncols>3:
            model=params[i][3]

        if combo:
            nlte=False    

        #set up input files, then run tlusty, finally check for nominal convergence and move all output files to 
        outdir=setup(log_qgrav, log_teff, log_dmtot, nlte, model, copy, tailname)
        if not outdir:
            continue
        run()
        maxchange=reltot()
        #Move tlusty output files to the appropriate directory
        outdir=clean(outdir,maxchange,nlte)
        if maxchange>1:
             continue

        #If we would like to calculate a combination on lte and nlte atmospheres
        if combo:
            nlte=True
            outdir=setup(log_qgrav, log_teff, log_dmtot, nlte, outdir, True, tailname)
            run()
            maxchange=reltot()
            #Move tlusty output files to the appropriate directory
            outdir=clean(outdir,maxchange,nlte)



   # return maxchange



##driver; parse user input
def main():
    parser=argparse.ArgumentParser(
        description='Wrapper for running TLUSTY')
    parser.add_argument('-n','--nlte',
        help='Switch on nlte',
        action='store_true')
    parser.add_argument('-c','--combo',
        help='Combination of lte and nlte models when reading from input file',
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
    parser.add_argument('-ft', '--tail',
        help='Stores name of file containing atomic data.',
        default='tail')
    parser.add_argument('-p', '--params',
        metavar=('qgrav', 'teff', 'dmtot'),
        help='Stores required parameters for atmosphere. Need 3 positional arguments in the order qgrav, teff, dmtot',
        nargs=3,
        type=float)
    args=parser.parse_args()

    myfile=args.file
    nlte=args.nlte
    model=args.model
    params=args.params
    copy=not args.nocopy
    combo=args.combo
    tailname=args.tail

    print(tailname)




    if params:
        tlusty_runs_input(params, model, nlte, copy, tailname)
    elif  model:
        print 'test'
        tlusty_runs_model(model, nlte, copy, tailname)
    else:
        tlusty_runs_file(myfile, nlte, copy, combo, tailname)


if __name__ == '__main__':
    main()
