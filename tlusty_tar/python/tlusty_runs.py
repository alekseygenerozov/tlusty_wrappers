#!/usr/bin/env python

import argparse
import numpy as np
import subprocess
import re
import astropy.table as at
from astropy.io import ascii
from types import*

##Run a command from the bash shell
def bash_command(cmd):
     process=subprocess.Popen(['/bin/bash', '-c', cmd],  stdout=subprocess.PIPE)
     return process.communicate()[0]
     # process.wait()
     # return process

##Routine to parse atmosphere
def parse_atm(f='fort.7'):
    end=re.compile('fort.*')
    f=re.sub(end,'',f)
    if (f):
        f=f+'/fort.7'
    else:
        f='fort.7'

    #Importing data from file
    atm=np.fromfile(f, sep=' ',dtype=float)
    nd=int(atm[0])
    blocks=int(atm[1])
    m=atm[2:2+nd]
    m=[m]
    m=np.transpose(m)
   
    #Reshaping data array
    atm=atm[2+nd:]
    atm=np.reshape(atm, (-1,blocks))
    #Concatenating mass depth grid with the actual vertical structure data
    atm=np.concatenate((m, atm), axis=1)
    return atm





##Checks atmospheric structures for density inversions greater than a certain threshold.
def inv_check(atm, thres=1.e-2):
    dens=atm[:,3]
    diff=np.diff(dens)
    dens=dens[:-1]
    #Fractional change in density
    frac=diff/dens
    frac=np.min(frac)
    print frac
    #If the 'largest' negative fractional change in density exceeds the specified threshold then return True, otherwise
    #return false
    if(frac<-thres):
        return True
    else:
        return False


##Function to parse file name in order to get list of parameters
def parse_file(f):
    refloat =r'[+-]?\d+\.?\d*'
    t=re.search("t"+refloat,f)
    t=re.search(refloat,t.group(0))
    
    m=re.search("m"+refloat,f)
    m=re.search(refloat,m.group(0))
    
    q=re.search("q"+refloat,f)
    q=re.search(refloat,q.group(0))

    params=np.array([float(t.group(0)), float(m.group(0)), float(q.group(0))])
    params=params/10.

    return params


##Go from a list of parameters to corresponding string
def reverse_parse(params):
    #Format strings with info on teff, dmtot, qgrav for creating folder
    t='{0:.8g}'.format(params[0]*10)
    m='{0:.8g}'.format(params[1]*10)
    q='{0:.8g}'.format(params[2]*10)

    f= 't' + t + 'm' + m + 'q' + q
    return f


##Setup input file for tlusty
def setup(log_qgrav, log_teff, log_dmtot, nlte, model, copy=True, tailname='tail', value=''):
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
        tailf=open(tailname, 'r')

        f1.write('1\n')

        qgrav=10.**log_qgrav
        teff=10.**log_teff
        dmtot=10.**log_dmtot

        print log_teff,log_dmtot,log_qgrav
   
        out='{0} {1:8.7e} {2:8.7e} {3:8.7e}'.format(0, teff, qgrav, dmtot)
        #print out
        # fort.5 input file for tlusty
        f5.write(out + '\n')
        #Second line of input
        f5.write(lte + '  ' + ltg + '\n')
        #Tail 
        tail=tailf.read()
        f5.write(tail) 

        #Format strings with info on teff, dmtot, qgrav for creating folder
        t='{0:.8g}'.format(log_teff*10)
        m='{0:.8g}'.format(log_dmtot*10)
        q='{0:.8g}'.format(log_qgrav*10)

        #print t
        #Create a folder in which to store the tlusty output in 
        loc='t' + t + 'm' + m + 'q' + q
        bash_command('mkdir -p ' + loc)


        #If the location of the input model is not blank then copy model atmosphere to the current directory
        if model:
            mparams=parse_file(model)

            #Check if unit 7 file is present in the specified in location.      
            process=bash_command('check ' + model + '/fort.7')
            #If it is not found then return empty string as flag
            if len(process)==0:
                loc=''
                return loc
            bash_command('cp ' + model + '/fort.7 ' + './')
            #If we are interpolating to a new mass scale then construct new mass grid for unit 8 file
            if np.abs((log_dmtot-mparams[1])/mparams[1])>1.e-2:
                bash_command('echo 70 0.001 ' + str(10**log_dmtot) + ' 0 0 > dmgrid.in')
                bash_command('rm dmgrid.out')
                bash_command('./dmgrid')
                bash_command('cat fort.7 dmgrid.out >fort.8')
            else:    
                bash_command('cat fort.7  >fort.8')

            #If the copy parameter has been set then we should copy the optional parameters file from the model location as well
            if copy:
                bash_command('cp ' + model + '/tmp.flag ' + './')

            #Amalgamating the old convergence log information into a single file
            bash_command('cat '+ model + '/fort.9_old ' + model + '/fort.9 > fort.9_old' )
        f=open('tmp.flag', 'a')
        if value:
            f.write(value+ '\n')
        #Return name of directory setup made   
        return loc



##Run TLUSTY from the bash shell
def run():
    cmd="./t202 <fort.5 >fort.6"
    bash_command(cmd)


# ##Gets nominal convergence from the convergence log file in tlusty
# def reltot(file='fort.9'):
#     dat=np.genfromtxt(file, skip_header=3)

#     niter= dat[-1, 0]
#     dat_cut=dat[(dat[:,0]==niter)]
#     #Check for non-numeric values in the convergence log
#     try:
#          change=dat_cut[:, 2:7]
#          change=np.array(map(np.abs, change))
#          maxchange=change.max()
#          #print maxchange
#     except ValueError:
#          maxchange=1000

#     #Using numpy's is nan function to check if the maximum change is a nan
#     if np.isnan(maxchange):
#          maxchange=1000
#     #Probably not good if the maximum is precisely zero
#     if maxchange==0:
#          return 1000

#     print niter, maxchange
#     return maxchange


##Function which parse tlusty convergence log
def converge_parse(f='fort.9'):
    #In case user decided to enter fort file name after directory name (i.e sanitizing user input file)
    suffix=re.compile('fort.*')
    f=re.sub(suffix,'',f)
    if (f):
        f=f+'/fort.9'
    else:
        f='fort.9'

    #Importing dat from from the specified file
    dat=np.genfromtxt(f, skip_header=3)
    #Split at specific locations (i.e. split convergence data by iteration)
    diff=np.diff(dat[:, 0])
    br=np.where((diff!=0))
    if(len(br[0])>0):
        br=br[0][0]
    else:
        br=69
    dat=np.reshape(dat, ((-1,br+1,9)))
    #Taking all of the values corresponding to the changes in state vector
    dat=dat[:,:,2:7]
    
    #Making all the log changes are positive; sometime log entries contain invalid numbers, in which case we would want to
    #raise an exception
    try:
        dat=np.abs(dat)
    except:
        raise

    #Getting the maximum change in state vector for each iteration
    chmax=np.empty([len(dat),2])
    for i in range(len(dat)):
        chmax[i]=np.array([i+1, dat[i].max()])
        if np.isnan(chmax[i, 1]):
            raise
    chmax=np.array(chmax)
    return chmax

##Check convergence log
def converge_check(f='fort.9',  thres=0.1):
    #Maximum change read from the convergence log; if an exception is raised then return false (i.e. the atmosphere
    #does not converge)
    try:
        chmax=converge_parse(f)
    except:
        return [False,False]

    end=min([4, len(chmax)])

    print chmax[-1,0],chmax[-1,1]
    conv=(chmax[-1, 1]<thres)

    #Check to see if the max relative change has been monotonically decreasing
    mono=True
    if (len(chmax)>=2):
        for i in range(1,end):
            if chmax[-i, 1]>1.01*chmax[-i-1,1]:
                mono=False

    bounded=True
    #Check to see if maximum relative change is bounded by by our threshold over the last few iterations
    for i in range(1, end+1):
        if chmax[-i, 1]>thres:
            bounded=False
    #Check if the relative change in state vector has been montonically non-increasing for the last few iterations
    #somewhat arbitrarily picked 4. If the total number of iterations is less than 4 then give the model the benefit of the doubt.
    conv2=(conv and mono) or bounded
    return [conv, conv2]



##Move all tlusty output files to the apropriate directory    
def clean(outdir, nlte, remove=False):
    #maxchange= np.log10(np.absolute(maxchange))
    [conv,conv2]=converge_check()
    atm=parse_atm('fort.7')
    nan_present=np.any(np.isnan(atm))
    if nan_present:
        print 'Warning NANs detected in the TLUSTY output atmosphere!'
        conv=conv2=False
    inv=inv_check(atm)

    prefix=''
    if not nlte:
        prefix='lte_'
    #constructing destination path
    dest='./'
    if conv2:
        if inv:
            dest=dest + outdir + '/'+prefix+'dens_inversion'
        else:   
            dest=dest + outdir + '/'+prefix+'converged'
    elif conv:
        if inv:
            dest=dest + outdir + '/'+prefix+'dens_inversion_marginal'
        else:   
            dest=dest + outdir + '/'+prefix+'converged_marginal'
    else:
        dest=dest + outdir + '/nconv'

    #check for existence of destination file
    i=2
    dest_orig=dest
    process=bash_command('check ' + dest)
    # If the remove flag has not been set and the destination already exists, we should keep trying extensions for destination 
    # until we obtain a destination name that does not already exist.
    while len(process)!=0 and not remove:
        dest=dest_orig+"_"+str(i)
        process=bash_command('check ' + dest)
        i+=1
    print dest

    #make the destination file if it doesn't already exist.
    bash_command('mkdir -p ' + dest)
    #If the remove flag has been set then remove everyting at the destination location.
    if remove:
        bash_command('rm -f ' + dest + '/*')
    #moving all tlusty i/o files to the destination
    bash_command('mv ' + 'fort* ' + dest)
    #move optional parameter file to destination
    bash_command('cp ' + 'tmp.flag ' + dest)

    return [conv2, dest]


##Run tlusty based on command line input parameters
def tlusty_runs_input(params, model, nlte=False, copy=True, combo=False, tailname='tail', remove=False, value=''):
    log_teff=params[0]
    log_dmtot=params[1]
    log_qgrav=params[2]

    print nlte

    outdir=setup(log_qgrav, log_teff, log_dmtot, nlte, model, copy, tailname,value)
    if not outdir:
        return
    run()
    #Move tlusty output files to the appropriate directory
    [conv,model]=clean(outdir, nlte, remove)       
    if not conv:
        return conv
    print model
    #If we would like to calculate a combination of lte and nlte atmospheres
    if combo:
        nlte=True
        outdir=setup(log_qgrav, log_teff, log_dmtot, nlte, model, True, tailname, value)
        if not outdir: 
            return outdir
        run()
        #Move tlusty output files to the appropriate directory
        [conv,model]=clean(outdir, nlte, remove)
    return conv



##Run tlusty based on parameters found at the location of model
def tlusty_runs_model(model, nlte=False, copy=True, tailname='tail', remove=False, value='',interp=False):
    #Extract parameters of model
    params=parse_file(model)
    log_teff=params[0]
    log_dmtot=params[1]
    log_qgrav=params[2]

    #Setup necessary input file for running and the output directory where we will store output files
    outdir=setup(log_qgrav, log_teff, log_dmtot, nlte, model, copy, tailname,value)
    if not outdir:
        return

    run()
    #Move tlusty output files to the appropriate directory
    clean(outdir, nlte, remove)
    #return [maxchange, dest]



##Construct tlusty model based on info in myfile
def tlusty_runs_file(myfile='params.in', nlte=False, copy=True, combo=False, tailname='tail', remove=False, value=''):
    params=ascii.read(myfile)
    ncols=len(params.columns)

    model=''

    #for each set of parameters in our input file
    for i in range(0, len(params)): 
        model=''
        f=open('fort.5', 'w')
        tailf=open('tail', 'r')
    
        log_teff=params[i][0]
        log_dmtot=params[i][1]
        log_qgrav=params[i][2]

        if ncols>3:
            model=params[i][3]
        if model:
            combo=False
        if combo:
            nlte=False    

        #set up input files, then run tlusty, finally check for nominal convergence and move all output files to 
        outdir=setup(log_qgrav, log_teff, log_dmtot, nlte, model, copy, tailname,value)
        if not outdir:
            continue
        run()
        #maxchange=reltot()
        #Move tlusty output files to the appropriate directory
        [conv,model]=clean(outdir, nlte, remove)       
        if not conv:
            continue

        #If we would like to calculate a combination of lte and nlte atmospheres
        if combo:
            nlte=True
            outdir=setup(log_qgrav, log_teff, log_dmtot, nlte, model, True, tailname)
            if not outdir:
                continue
            run()
            #maxchange=reltot()
            #Move tlusty output files to the appropriate directory
            clean(outdir,nlte, remove)


##Driver; parse user input
def main():
    parser=argparse.ArgumentParser(
        description='Wrapper for running TLUSTY')
    parser.add_argument('-n','--nlte',
        help='Switch on nlte',
        action='store_true')
    parser.add_argument('-c','--combo',
        help='Combination of lte and nlte models when reading from input file or using command line parameters',
        action='store_true')
    parser.add_argument('-r', '--remove',
        help='If model already exists at destination overwrite existing model',
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
        metavar=('teff', 'dmtot', 'qgrav'),
        help='Stores required parameters for atmosphere. Need 3 positional arguments in the order teff, dmtot, qgrav',
        nargs=3,
        type=float)
    parser.add_argument('-val', '--value',
        help='Sets param value to be value in optional parameters file',
        default='')
    args=parser.parse_args()

    #Extract user inputted parameters
    myfile=args.file
    nlte=args.nlte
    model=args.model
    params=args.params
    copy=not args.nocopy
    combo=args.combo
    tailname=args.tail
    remove=args.remove
    value=args.value



    if params:
        tlusty_runs_input(params, model, nlte, copy, combo, tailname, remove, value)
    elif  model:
        tlusty_runs_model(model, nlte, copy, tailname, remove, value)
    else:
        tlusty_runs_file(myfile, nlte, copy, combo, tailname, remove, value)


if __name__ == '__main__':
    main()
