#!/usr/bin/env python

from scipy.interpolate import griddata
import numpy as np
import numpy.ma as ma
import re

from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.animation as animation

import argparse
import subprocess

import graybody as gray
import tlusty_runs as tr
import warnings

import shlex
from types import*



##Defining physical constants
G=6.67*10**-8
c=3*10**10
h=6.67*10**-27
kb=1.38*10**-16
M_sun=2.*10**33



##Run a command from the bash shell
def bash_command(cmd):
    process=subprocess.Popen(['/bin/bash', '-c',cmd],  stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    return process.communicate()[0]
    # process.wait()
    # return process

##Function to extract frequency from wavelength in angstroms
def get_freq(w):
    return c*10**8/w
##Inverse of get_freq
def get_w(nu):
    return c*10**8/nu


# Function to parse file name in order to get list of parameters.
def parse_file(f, string=False):
    f2=re.sub('.14\Z', '', f)
    refloat =r'[+-]?\d+\.?\d*'
    t=re.search("t"+refloat,f2)
    t=re.search(refloat,t.group(0))
    
    m=re.search("m"+refloat,f2)
    m=re.search(refloat,m.group(0))
    
    q=re.search("q"+refloat,f2)
    q=re.search(refloat,q.group(0))
    
    if string:
        params=(t.group(0), m.group(0), q.group(0))
    else:
        params=(float(t.group(0))/10., float(m.group(0))/10., float(q.group(0))/10.)
    
    return params

##Go from a list of parameters to corresponding string
def reverse_parse(params):
    #Format strings with info on teff, dmtot, qgrav for creating folder
    t='{0:.8g}'.format(params[0]*10)
    m='{0:.8g}'.format(params[1]*10)
    q='{0:.8g}'.format(params[2]*10)

    f= 't' + t + 'm' + m + 'q' + q
    return f


# Extract spectrum from unit 14 type output file, if mu  is set to -1 then get the flux
def get_spec(file,  mu=-1, nmu=10):
    #Read file assumed to be in the format of tlusty uni 14 file. File has lines with wavelengths and h followed by blocks with 
    #intensities and polariztions
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        spec=np.genfromtxt(file, skip_header=1, invalid_raise=False)
        wlh=np.genfromtxt(file, invalid_raise=False)
    
    nfreq=len(wlh)
    spec=spec.flatten()
    wlh=wlh.flatten()
    #Skip over fluxes and polarizations
    spec=spec[::2]
    wl=wlh[::2]
    #Second moment of radiative intensity
    h=wlh[1::2]
    flux=4*np.pi*h
    flux.shape=(nfreq,1)

    #Add flux to end of intensity arrays for every frequency
    spec.shape=(nfreq, nmu)
    spec=np.append(spec, flux, axis=1)

    return (wl,spec) 



# Regrid spectrum in wavelength space. Want all spectra to be on the same grid in wavelength space when we interpolate; In frequency the default limits for the 
# the wavelength correspond to 2.4e14 and 5e19 Hz.
def regrid(spec, wlo=30, whi=3.e5, nws=300):
    wgrid=np.log(wlo)+np.log(whi/wlo)*np.arange(0, nws)/(nws-1)
    newspec=griddata(np.log(spec[0]), spec[1], wgrid, fill_value=1.e-36,method='linear')
    newspec=newspec[::-1]

    wgrid=np.exp(wgrid)
    freq=map(get_freq, wgrid)
    freq=freq[::-1]
    
    freq2=np.empty_like(newspec)
    for i in range(len(freq2)):
        freq2[i].fill(freq[i])
    spec=np.vstack([freq2, newspec])
    spec.shape=(2,nws,11)

    return spec



# Function calculates brightness temperature for a given frequncy and intensity
def Tb(freq, intens, fcol=2):
    return np.log10(h*freq/kb/fcol/np.log(1+(2*h*freq**3/fcol**4/intens/c/c)))


# Invert brightness temperature to obtain the intensity
def Tb_inv(freq, Tb, fcol=2):
    if np.any(np.isinf(np.exp(h*freq/kb/fcol/10**Tb))):
        return 0
    try:
        return (2./fcol**4)*(h*(freq**3)/c/c)/(np.exp(h*freq/kb/fcol/10**Tb)-1)
    except:
        return 0

    
# # Get brightness temperature from a spectrum, with default color correction factor of 2 as described in ApJS 164 530D
# # returns the log of the brightness temperature
# def bright_spec(spec, fcol=2):
#     freq=spec[0]
#     intens=spec[1]
#     farr=np.empty(len(freq))
#     farr.fill(fcol)
    
#     bright=map(Tb, freq, intens, farr)
#     return np.vstack((freq, bright))


# # Invert brightness temperature to obtain the spectrum
# def bright_inv(spec_inv, fcol=2):
#     freq=spec_inv[0]
#     temp=spec_inv[1]
    
#     farr=np.empty(len(freq))
#     farr.fill(fcol)
    
#     intens=map(Tb_inv, freq, temp, farr)
#     spec=np.vstack((freq, intens))
#     return spec

##Convert spectral information to photometric colors; note that this will in general be red-shift dependent
def to_colors(spec, z=0):
    bands=[[7.5*10**14,1.*10**15], [5.45455*10**14,7.5*10**14], [ 
  4.28571*10**14,5.45455*10**14], [3.52941*10**14,4.28571*10**14], [ 
  2.5*10**14,3.75*10**14]]
    colors=np.empty(5)
    for i in range(len(bands)):
        #Boundaries of color band
        nu1=bands[i][0]*(1+z)
        nu2=bands[i][1]*(1+z)
        #Filter to get band of interest
        filt=(spec[0]>nu1) & (spec[0]<nu2)

        nu=spec[0,filt]
        dnu=np.diff(nu)
        lum=spec[1,filt]
        print len(lum)
        #Take average of luminosities on either sides of frequency bin
        lum=(lum[0:-1]+lum[1:])/2
        print len(lum)
        #Integrating luminosity
        colors[i]=np.sum(lum*dnu)

    return np.log10(colors[0:-1]/colors[1:])

##Constructing table of spectra from unit 14 files
def construct_table(models, logi=False):
    models=np.genfromtxt(models, dtype='string')
    params=map(parse_file,models)
    params=np.array(params)
    params=params[:, ::2]
    
    spec=map(get_spec, models)
    spec=map(regrid, spec)
    spec=np.array(spec)

    if logi:
        spec[:, 1]=np.log10(spec[:, 1])
        spec=np.array(spec)
    else:
        spec[:, 1]=map(Tb, spec[:,0],spec[:,1])
        spec=np.array(spec)

    return (params, spec)

##Get parameters corresponding to our disk
def get_params(file):
    return np.genfromtxt(file)

##Use table construct spectra from list of parameters
def params_to_spec(params, table, method='', logi=False, mu=2):
    intens=table[1]
    nu=table[1][0,0]
    #Interpolating based on table (interpolating in teff/qg space)
    print table[0].shape,intens.shape,params.shape
    if method:
        grid=griddata(table[0], intens, params, method=method)
    else:
        grid=griddata(table[0], intens, params)

    #Angle space
    #Flag to return flux for all angles
    if mu>1:
        grid2=grid
    #Flag to return flux
    elif mu<0:
        grid2=grid[:,:,:,-1] 
        nu=nu[:,-1]    
    #Otherwise interpolate the spectrum in mu.
    else:
        muarr=np.array([0.0130467357,0.0674683167,0.160295216,\
            0.283302303,0.425562831,0.574437169, 0.716697697,\
            0.839704784,0.932531683,0.986953264])

        nu=nu[:,0]
        grid2=np.empty((len(grid),2, len(nu)))
        for i in range(len(grid)):
            grid2[i,0]=nu
            grid2[i,1]=griddata(muarr, np.transpose(grid[i,1,:,:-1]), [mu])


    #Make sure spectrum is converted back to intensity space, also store info about bad spectra (i.e. those that fall outside
    #the range of the table used.)
    vTb_inv=np.vectorize(Tb_inv)
    good=np.empty(len(grid2))
    good.fill(True)
    for i in range(len(grid2)):
        #print np.any(np.isnan(grid2[i]))
        if np.any(np.isnan(grid2[i])):
            good[i]=False
            grid2[i,0]=nu
            grid2[i,1].fill(params[i,0])
            fcol=np.ones_like(nu)
            
            grid2[i,1]=vTb_inv(nu, grid2[i,1], fcol)
        elif logi:
            grid2[i,1]=10.**grid2[i,1]
        else:
            print grid2.shape
            grid2[i,1]=vTb_inv(nu, grid2[i,1])

    #Return the interpolated spectra for parameters within our grid
    return (grid2,good)


##Compares tlusty spectrum to one that is interpolated from a table 
def test_spec(f, table=[], tablef='tmpd', method='', logi=False):
    if table==[]:
        table=construct_table(tablef, logi=logi)
    #Read tlusty spectrum from file. 
    testspec=get_spec(f)
    testspec=regrid(testspec)
    testspec=testspec[:,:,-1]
    
    #Construct interpolated spectrum for the params corresponding to file.
    params=parse_file(f)
    params_s='t'+str(10*params[0])+'m'+str(10*params[1])+'q'+str(10*params[2])
    print params
    params=params[::2]
    testspec_interp=params_to_spec([params], table, method=method, logi=logi, mu=-1)
    testspec_interp=testspec_interp[0]

    if np.any(np.isnan(testspec_interp)):
        print "Warning -- unable to interpolate spectrum for specified parameters."
        return [testspec[0], testspec[0]*testspec[1], testspec_interp[0]*testspec_interp[1]]
    return [testspec[0], testspec[0]*testspec[1], testspec_interp[0]*testspec_interp[1]]


#takes list of models creates animation of computed spectrum plus interpolated spectrum
def animate_test_spec(models, table=[], tablef='tmpd', method='', logi=False):
    def update_img(n):
        spec_plot.set_ydata(spec[n])
        spec_plot_interp.set_ydata(spec_interp[n])
        label.set_text('log(teff)= '+str(teffs[n]))
        return [spec_plot, spec_plot_interp]

    fig, ax=plt.subplots()
    plt.loglog()
    plt.xlabel(r"$\nu$ [hz]")
    plt.ylabel(r"$\nu F_{\nu}$ [ergs s$^{-1}$ cm$^{-2}$ ]")
    plt.axis([10.**14, 2*10.**18, 10.**6, 10.**16])
    label=ax.text(0.02, 0.95, '', transform=ax.transAxes)

    (spec, spec_interp, teffs)=([],[],[])

    nu=test_spec(models[0], table=table, tablef=tablef, method=method, logi=logi)[0]
    tmpspec=(test_spec(models[0], table=table, tablef=tablef, method=method, logi=logi)[1])
    tmpspec_interp=(test_spec(models[0], table=table, tablef=tablef, method=method, logi=logi)[1])
    print nu.shape,tmpspec.shape
    spec_plot,=ax.plot(nu, tmpspec)
    spec_plot_interp,=ax.plot(nu, tmpspec_interp)


    (ts,ms,qs)=parse_file(models[0], True)
    writer = animation.writers['ffmpeg'](fps=10) 
    for m in models:
        m=m.rstrip()
        (ts2,ms2,qs2)=parse_file(m, True)
        if not(ms2==ms) or not(qs2==qs):
            print len(spec)
            ani = animation.FuncAnimation(fig,update_img,len(spec),interval=50)
            ani.save('interp_m'+ms+'_q'+qs+'.mp4',writer=writer,dpi=100)
            (spec, spec_interp, teffs)=([],[],[])
            (ms,qs)=(ms2,qs2)

        #spec=test_spec(m, table=table, tablef=tablef, method=method, logi=logi) 
        teffs.append(float(ts2)/10)
        spec.append(test_spec(m, table=table, tablef=tablef, method=method, logi=logi)[1])
        spec_interp.append(test_spec(m, table=table, tablef=tablef, method=method, logi=logi)[2])


##Computes composite spectrum given array of radii, spectra. Also computes corresponding blackbody spectrum using list of Teff that are passed to the function
def sum_spec(r, specs, Teff, Qg, mu, M=10.**6): 
    r=r*(G*M*M_sun/c**2)  
    lr=np.log10(r)
    #For now assuming a logarithmically evenly spaced grid--for simplicity
    dlr=np.diff(lr)[0]
    #Calculate the dr's
    dr=np.empty_like(lr)
    for i in range(len(lr)):
        dr[i]=10.**(lr[i]+(dlr/2))-10.**(lr[i]-(dlr/2))
    
    valid=specs[1]
    print valid
    specs=specs[0]
    nu=specs[0,0]
    f=specs[:, 1]
    print Teff,Qg

    bb=np.zeros_like(f[0])
    gb=np.zeros_like(f[0])
    L=np.zeros_like(f[0])
    L2=np.zeros_like(f[0])
    for i in range(len(r)):
        #Check for any invalid entries -- these correspond to the parameters outside the edge of our table
        # valid=not np.any(np.isnan(specs[i]))
        rad=2*np.pi*r[i]*dr[i]
        #For every frequncy under consideration
        for j in range(len(nu)):
            #If we have a valid interpolated spectrum...
            #tmpbb=np.empty(len(np.atleast_1d(nu[0])))
            #tmpgb=np.empty(len(np.atleast_1d(nu[0])))
            tmpbb=(np.pi*rad*(Tb_inv(nu[j], Teff[i], fcol=1)))
            tmpgb=(rad*(gray.gb(nu[j], 10.**Teff[i], 10.**Qg[i])))
            if (0<mu<1):
                tmpbb=tmpbb/np.pi
                tmpgb=tmpgb/np.pi

            bb[j]+=tmpbb
            gb[j]+=tmpgb
            if valid[i]:
                L[j] =(L[j] +rad*f[i,j])
                L2[j]=(L2[j]+rad*f[i,j])
            #Otherwise add bb flux to L but not L2
            else:
                L[j]+=(tmpbb)

    return (nu, gb, bb, L, L2) 


##Calculates a composite disk spectrum given an file containing input radial parameters.
def disk_spec(f, table=[], tablef='tmpd', method='', logi=False, mu=0.6):
    #Construct table if necessary
    if table==[]:
        table=construct_table(tablef, logi=logi)
    #For consistency with disk_spec_gr. Not only M will do anything here
    #Get mass, spin, inclination from the header line 
    global_params=dict({'M':1.e6})
    header=bash_command('head -1 '+f)
    header=shlex.split(header)
    for p in header:
        p2=p.split('=')
        if p2[0] in global_params:
            global_params[p2[0]]=float(p2[1])
    M=global_params['M']
    a=0
    #mu=-1

    disk_params=np.genfromtxt(f, skip_header=1)
    specs=params_to_spec(disk_params[:, 1:4:2], table, method=method, logi=logi, mu=mu)

    #Disk parameters
    r=disk_params[:, 0]
    Teff=disk_params[:, 1]
    Qg=disk_params[:, 3]
    #Finding the total flux, for the given parameters
    totf=sum_spec(r, specs, Teff, Qg, mu, M=M)
    # nu=totf[0][:,mu]
    # totfg=totf[1][:,mu]
    # totfb=totf[2][:,mu]
    # totft=totf[3][:,mu]
    # totft2=totf[4][:,mu]
    inc_factor=1
    if 0<mu<1:
        inc_factor=mu*4*np.pi
    nu=totf[0]
    totfg=inc_factor*totf[1]
    totfb=inc_factor*totf[2]
    totft=inc_factor*totf[3]
    totft2=inc_factor*totf[4]

    # #r=r*(G*M*M_sun/c**2)  
    # lr=np.log10(r)
    # #For now assuming a logarithmically evenly spaced grid--for simplicity
    # dlr=np.diff(lr)[0]
    # dr=np.empty_like(r)
    # for i in range(len(lr)):
    #     dr[i]=10.**(lr[i]+(dlr/2))-10.**(lr[i]-(dlr/2))
    fig,ax=plt.subplots(nrows=1, ncols=1, figsize=(6,6), subplot_kw=dict(adjustable='datalim'))
    # # #plt.title(str(bin_params[0])+" "+str(bin_params[1])+" "+str(bin_params[2]))
    plt.xlabel(r"$\nu$ [hz]")
    #Plotting the composite disk spectrum
    ax.set_ylabel(r"$\nu L_{\nu}$ [ergs s$^{-1}$]")

    ax.set_xlim(10.**14, 10.**17)
    ax.set_ylim(10.**41, 10.**46)
    ax.set_xscale('log')
    ax.set_yscale('log')

    #ax.plot(nu, nu*totfg, label ='graybody')
    ax.plot(nu, nu*totfb,'r-', label ='blackbody')
    ax.plot(nu, nu*totft,'b-',label ='tlusty +\nblackbody')
    #ax.plot(nu, nu*totft2,label ='tlusty')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1])

    outfile='sp_M'+'{0:.3e}'.format(M)+'_a'+'{0:.3f}'.format(a)+'_mu'+'{0:.3f}'.format(mu)
    np.savetxt(outfile, np.transpose([nu, totft]))
    # #Plotting the contributions of individual annuli
    # nu=specs[0,0]
    # ax[1].set_ylim(10**-5,10)
    # ax[1].set_xlim(10**14, 10**17)
    # ax[1].set_xscale('log')
    # ax[1].set_yscale('log')
    # #ax[1].set_ylabel(r"F$_{\nu}$ [ergs s$^{-1}$ cm$^{-2}$]")
    # wl=np.array(map(get_w, nu))
    # valid=np.empty(len(specs))
    # for i in range(len(specs)):
    #     specs[i,1]=2*np.pi*r[i]*dr[i]*specs[i,1]
    #     ax[1].plot(nu, specs[i, 1])
    #     valid[i]=not np.any(np.isnan(specs[i,1]))
    
    # ax[0].plot(nu,2.6*np.mean(specs[valid==1,1], axis=0))
    # ax[1].plot(wl,10**6*np.mean(specs[valid==1,1]/wl**2, axis=0))
    plt.close()
    return fig


##Given a radial disk profile, calculates the spectra of each of the annuli, and then call the kerrtrans9 routine
def disk_spec_gr(f, table=[], tablef='tmpd', method='', logi=False, fobs=[1.e14, 1.e17], rmax=1.e4, mu=0.6):
    #Construct table if necessary
    if table==[]:
        table=construct_table(tablef, logi=logi)

    #Get mass, spin, inclination from the header line 
    global_params=dict({'M':1.e6, 'a':0., 'mu':0.6})
    header=bash_command('head -1 '+f)
    header=shlex.split(header)
    for p in header:
        p2=p.split('=')
        if p2[0] in global_params:
            global_params[p2[0]]=float(p2[1])
 
    M=global_params['M']
    #mu=global_params['mu']
    a=global_params['a']

    
    #Get the radial profile of the disk 
    disk_params=np.genfromtxt(f, skip_header=1)
    r=disk_params[:, 0]
    #Radial cut
    cut=(r<rmax)
    disk_params=disk_params[cut]
    #Extract spectra corresponding to disk parameters
    specs=params_to_spec(disk_params[:, 1:4:2], table, method=method, logi=logi, mu=2)
    specs=specs[0]


    nu=specs[0,0,:,-1]
    r=disk_params[:, 0]
    Teff=disk_params[:, 1]
    Qg=disk_params[:, 3]


    #Creating an input file for kerrtrans9
    bash_command('rm tmp.in')
    kerr_in=open('tmp.in','a')
    #Writing frequency information...
    np.savetxt(kerr_in, [-len(nu)], fmt='%i')
    np.savetxt(kerr_in, [[nu[0],nu[-1]]], fmt='%7.5e')
    np.savetxt(kerr_in, [-len(nu)], fmt='%i')
    np.savetxt(kerr_in, [fobs], fmt='%7.5e')
    #Writing disk parameter information...
    np.savetxt(kerr_in, [mu], fmt='%f')
    np.savetxt(kerr_in, [[M, 0., a]],fmt='%7.5e %f %f')
    np.savetxt(kerr_in,[[r[-1],-1.,len(r),4,1]], fmt='%f %f %i %i %i')
    kerr_in.close()

    #Creating spectral input file for kerrtrans9
    bash_command('rm emrad.in')
    emrad=open('emrad.in','a')
    for i in range(len(r)):
        np.savetxt(emrad, [[r[i], len(nu)]], fmt='%7.5f %i')
        for j in (range(len(nu)))[::-1]:
            np.savetxt(emrad, [[get_w(nu[j]), specs[i,1,j,-1]/(4*np.pi)]], fmt='%7.5e')
            intens=specs[i,1,j,:-1]
            pol=np.zeros_like(intens)

            intens=np.transpose([intens,pol])
            intens=intens.flatten()
            intens.shape=(2,10)
            np.savetxt(emrad,intens, fmt='%7.5e')
    emrad.close()

    #Running kerrtrans9 now that input files tmp.in and emrad.in have been generated. Output spectrum to file labelled by
    #disk parameters.
    outfile='sp_M'+'{0:.3e}'.format(M)+'_a'+'{0:.3f}'.format(a)+'_mu'+'{0:.3f}'.format(mu)+'_gr'
    bash_command('./kerrtrans9 <tmp.in >'+outfile)

    fig,ax=plt.subplots(nrows=1, ncols=1, figsize=(6,6), subplot_kw=dict(adjustable='datalim'))
    # # #plt.title(str(bin_params[0])+" "+str(bin_params[1])+" "+str(bin_params[2]))
    # #plt.xlabel(r"$\nu$ [hz]")
    #Plotting the composite disk spectrum
    ax.set_ylabel(r"$\nu L_{\nu}$ [ergs s$^{-1}$]")
    ax.set_xlim(10.**14, 10.**17)
    ax.set_ylim(10.**35, 10.**42)
    ax.set_xscale('log')
    ax.set_yscale('log')

    nu=np.genfromtxt(outfile, usecols=0)
    totf=np.genfromtxt(outfile, usecols=1)
    ax.plot(nu, nu*totf)

    return fig


def main():
    parser=argparse.ArgumentParser(
        description='Either takes list of disk parameters and computes composite disk spectrum from table or compares spectra interpolated'+
          ' from table to those in a test directory')
    parser.add_argument('-d', '--disk',
        help='file containing list of disk profile files ',
        default='')
    parser.add_argument('-mu', '--mu',
        help='cosine of inslincation angle for the case of a disk ',
        type=float,
        default=0.6)
    parser.add_argument('-t', '--test',
        help='directory containing the location of models to test',
        default='')
    parser.add_argument('-sf', '--skip',
        help='in case d argument is used this can be used to specify file name containing the names of models to be skipped',
        default='')
    parser.add_argument('-tf', '--tablefile',
        help='name of file containing list of table models',
        default='tmpd')
    parser.add_argument('-m', '--method',
        help='Specify order of interpolation should be linear for table',
        choices=['linear', 'cubic', 'nearest'])
    parser.add_argument('-li', '--logi',
        help='Specifies that the interpolation should be done in terms of log intensity instead of brightness temp.',
        action='store_true')
    parser.add_argument('-gr', '--gr',
        help='Specifies that we should use relativistic transfer function to calculate the composite disk spectrum.',
        action='store_true')
    parser.add_argument('-rmax', '--rmax',
        help='maximum radius of disk ',
        type=float,
        default=1.e4)
    # parser.add_argument('-a', '--animate',
    #     help='For the case of test spectra, specifies that a movie should be made rather than a static pdf.'
    #     a)

    args=parser.parse_args()
    t=args.test
    d=args.disk
    method=args.method
    tablef=args.tablefile
    logi=args.logi
    skip=args.skip
    mu=args.mu
    gr=args.gr
    rmax=args.rmax

    if d:
        table=construct_table(tablef, logi=logi)
        #pdf_pages = PdfPages('composite.pdf')
        param_files=np.genfromtxt(d, dtype=str)
        param_files=np.atleast_1d(param_files)
        for i in range(len(param_files)):
            #print param_files
            if gr:
                fig=disk_spec_gr(param_files[i], table=table, tablef=tablef, method=method,logi=logi,rmax=rmax,mu=mu)
                fig.savefig('composite_'+str(i)+'_gr'+'.png')
            else:
                fig=disk_spec(param_files[i], table=table, tablef=tablef, method=method,logi=logi,mu=mu)
                fig.savefig('composite_'+str(i)+'.png')
            #pdf_pages.savefig(fig)
        #pdf_pages.close()
    elif t:
        table=construct_table(tablef, logi=logi)
        # pdf_pages = PdfPages('interp_test.pdf')
        process=bash_command('echo '+t+'/*14')

        # bash_command('rm interp_log')
        # logfile=open('interp_log', 'a')

        # deviation_list=np.empty(0)
        #models=process.stdout.readlines()[0]
        models=shlex.split(models)
        if skip:
            skip_models=np.genfromtxt(t+'/'+skip, dtype=str)
            models=np.setdiff1d(models, skip_models)
        models=np.array(models, dtype=str)
        #Sort models used by teff, then m, and finally q
        params=map(parse_file, models)
        dtype=[('teff',float), ('m', float), ('q', float)]
        params=np.array(params, dtype=dtype)
        order=np.argsort(params, order=['q', 'm', 'teff'])
        models=models[order]

        print models

        #Create animations comparing the interpolated spectra to tlusty spectra found in the specified directory
        animate_test_spec(models, table=table, tablef=tablef, method=method, logi=logi)

        
        # for o in order:
        #     params2=params[o]
        #     print params2
        #     diff=np.diff(params2[:, 0])
        #     br=np.where((diff>10**-6)

        #     print br

            # print np.split(params2, br)
        
        #animate_test_spec(models, table=table, tablef=tablef, method=method, logi=logi)


        # params=np.array(map(tr.parse_file, models))
        # params=map(tr.parse_file, models)
        # params=np.array(params)
        # teffs=params[:,0]
        # order=np.argsort(teffs)

        # models[order]
        # for m in models:
        #     m=m.rstrip()
        #     spec=test_spec(m, table=table, tablef=tablef, method=method, logi=logi)
        #     fig=spec[0]
        #     deviation=spec[1]
        #     logfile.write(m+" "+str(deviation)+"\n")
        #     if deviation!=-1:
        #         deviation_list=np.append(deviation_list,deviation)
        #     #Save file 
        #     pdf_pages.savefig(fig)

        # print deviation_list
        # pdf_pages.close()
        # fig2=plt.figure()
        # plt.hist(deviation_list, bins=50, normed=1)
        # fig2.savefig('dev_hist.pdf')
    else:
        parser.print_help()

    #pdf_pages.close()



if __name__ == '__main__':
    main()
    
    

    

# def Knu(es, f):
#     return  np.sqrt(-4.*(-1. + 2.*es - 1.*es**2)**6 + \
#     (2. - 12.*es + 30.*es**2 - 40.*es**3 + 30.*es**4 - 12.*es**5 + \
#     2.*es**6 - 27.*es**2*f + 108.*es**3*f - 162.*es**4*f + 108.*es**5*f - \
#     27.*es**6*f)**2)
    

# <codecell>

# table=construct_table('tmpd')
# disk_params=get_params('dparams_pert_total')


# specs=params_to_spec(disk_params, table)

# loglog()
# figsize(20, 8)

# r=disk_params[:, 0:2]
# Teff=disk_params[:, 2]


# totf=sum_spec(r, specs, Teff)
# totf1=totf[1]
# totf2=totf[0]

# grayf=np.genfromtxt("graydisk")


# peak=10.**44
# #peak=np.max([max1, max2])
# #totf2=sum_spec(r[2:], specs[2:])
# #totf3=sum_spec(r[5:], specs[5:])

# loglog()
# figsize(20, 8)

# xlabel("nu [hz]")
# ylabel("nu L_nu [ergs/s]")

# axis([10.**14, 2*10.**18, 10**-6*peak, 2*peak]) 

# plot(totf1[0], totf1[0]*totf1[1])
# plot(totf2[0], totf2[0]*totf2[1])
# plot(grayf[:, 0], grayf[:, 0]*grayf[:, 1])
# #plot(totf2[0], totf2[0]*totf2[1])
# #plot(totf3[0], totf3[0]*totf3[1])

# savefig("composite_pert.pdf")




   # colors=np.empty([5,4])
    # for z in range(0,5):
    #     colors[z]=to_colors(np.array([nu, totft]), z=z)


    # fig=plt.figure()
    # ax=fig.add_subplot(221)
    # ax.set_xlim(-1,5)
    # ax.set_ylim(-1,3)
    # ax.plot(colors[0,:], colors[1,:], 'rs')

    # ax=fig.add_subplot(222)
    # ax.set_xlim(-1,3)
    # ax.set_ylim(-1,3)
    # ax.plot(colors[1,:], colors[2,:], 'rs')

    # ax=fig.add_subplot(223)
    # ax.set_xlim(-1,3)
    # ax.set_ylim(-1,2)
    # ax.plot(colors[2,:], colors[3,:], 'rs')
    # plt.show()
    # plt.close()
    # np.savetxt('spec', np.transpose(np.array([nu, totft])))