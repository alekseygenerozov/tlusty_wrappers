#!/usr/bin/env python

from scipy.interpolate import griddata
import numpy as np
import numpy.ma as ma
import re

import matplotlib as mpl
#from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
import colorsys

import argparse
import subprocess

import graybody as gray
import tlusty_runs as tr
import warnings

import shlex
from types import*

import colorsys





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


#Getting all files corresponding to a pattern pat
def get_spec_pat(pat, mu=-1):
    files=bash_command('echo '+pat)
    files=shlex.split(files) 
    spec=[]
    params=[]
    for f in files:
        params.append(parse_file(f))
        spec.append(regrid(get_spec(f, mu=-1)))
    spec=np.array(spec)
    params=np.array(params)
    spec=spec[:,:,:,-1]
    return [spec, params]

#Generating rainbow of colors of size n
def rainbow(n=10):
    sat=1.
    value=1.
    cols=np.empty([n, 3])
    for i in range(len(cols)):
        hue=1./(i+1)
        cols[i]=np.array(colorsys.hsv_to_rgb(hue, sat, value))

    return cols

#Plot flux of all files matching a pattern
def plotf(pat):
    [spec, params]=get_spec_pat(pat)
    fig, ax=plt.subplots(figsize=(12,8))
    ax.set_xscale('log')
    ax.set_yscale('log')
    cols=rainbow(n=len(spec))

    peak=0.
    for i in range(len(spec)):
        peak=np.max([peak, np.max(spec[i,0]*spec[i,1])]) 
        ax.plot(spec[i,0], spec[i,0]*spec[i,1], color=cols[i])
    plt.axis([1.E14, 3.E16, 10.**-4*peak, 4*peak])
    # plt.close()
    return fig


# Regrid spectrum in wavelength space. Want all spectra to be on the same grid in wavelength space when we interpolate; In frequency the default limits for the 
# the wavelength correspond to...If the keep keyword argument is set, then we simply keep the old frequency grid.
def regrid(spec, wlo=30, whi=3.e5, nws=300, keep=False):
    #Wavelength grid onto which we would like to interpolate. 
    wgrid=[1.47712, 1.49495, 1.51278, 1.53061, 1.54844, 1.56627, 1.5841, \
1.60193, 1.61976, 1.63759, 1.65542, 1.67325, 1.69108, 1.70891, \
1.72674, 1.74457, 1.7624, 1.78023, 1.79806, 1.81589, 1.83372, \
1.85155, 1.86938, 1.88721, 1.90504, 1.92288, 1.94071, 1.95854, \
1.97637, 1.9942, 2.01203, 2.02986, 2.04769, 2.06552, 2.08335, \
2.10118, 2.11901, 2.13684, 2.15467, 2.1725, 2.19033, 2.20816, \
2.22599, 2.24382, 2.26165, 2.27948, 2.29731, 2.31514, 2.33297, \
2.3508, 2.36231, 2.36961, 2.37691, 2.38421, 2.39151, 2.39881, 2.4061, \
2.4134, 2.4207, 2.428, 2.4353, 2.4426, 2.4499, 2.4572, 2.4645, \
2.4718, 2.4791, 2.4864, 2.4937, 2.50099, 2.50829, 2.51559, 2.52289, \
2.53019, 2.53749, 2.54479, 2.55209, 2.55939, 2.56669, 2.57399, \
2.58129, 2.58859, 2.59589, 2.60318, 2.61048, 2.61778, 2.62508, \
2.63238, 2.63968, 2.64698, 2.65428, 2.66158, 2.66888, 2.67618, \
2.68348, 2.69078, 2.69807, 2.69865, 2.70733, 2.71463, 2.72193, \
2.72923, 2.73653, 2.74383, 2.75113, 2.75843, 2.76573, 2.77302, \
2.78032, 2.78762, 2.79492, 2.80222, 2.80952, 2.81682, 2.82412, \
2.83142, 2.83872, 2.84602, 2.85332, 2.86062, 2.86791, 2.87521, \
2.88251, 2.88981, 2.89711, 2.90441, 2.91171, 2.91901, 2.92631, \
2.93361, 2.94091, 2.94821, 2.95551, 2.95568, 2.96437, 2.97167, \
2.97897, 2.98627, 2.99357, 3.00087, 3.00816, 3.01546, 3.02276, \
3.03006, 3.03736, 3.04466, 3.05196, 3.05926, 3.06656, 3.07386, \
3.08116, 3.08846, 3.09576, 3.10305, 3.11035, 3.11765, 3.12495, \
3.13225, 3.13955, 3.14685, 3.15415, 3.16145, 3.16875, 3.17605, \
3.18335, 3.19065, 3.19795, 3.20524, 3.21254, 3.21984, 3.22714, \
3.23444, 3.24174, 3.24904, 3.25634, 3.26364, 3.27094, 3.27824, \
3.28554, 3.29284, 3.30013, 3.30743, 3.30787, 3.31655, 3.32385, \
3.33115, 3.33845, 3.34575, 3.35305, 3.36035, 3.36765, 3.37495, \
3.38225, 3.38954, 3.39684, 3.40414, 3.41144, 3.41874, 3.42604, \
3.43334, 3.44064, 3.44794, 3.45524, 3.46254, 3.46984, 3.47714, \
3.48443, 3.49173, 3.49903, 3.50633, 3.51363, 3.52093, 3.52823, \
3.53553, 3.54283, 3.55013, 3.55743, 3.55774, 3.56643, 3.57373, \
3.58103, 3.58833, 3.59563, 3.60293, 3.61022, 3.61752, 3.62482, \
3.63212, 3.63942, 3.64672, 3.65402, 3.66132, 3.66862, 3.67592, \
3.68322, 3.69052, 3.69782, 3.70511, 3.71241, 3.71971, 3.72701, \
3.73431, 3.74161, 3.74891, 3.75621, 3.76351, 3.77081, 3.77811, \
3.78541, 3.79271, 3.80001, 3.8073, 3.8146, 3.8219, 3.8292, 3.8365, \
3.8438, 3.8511, 3.8584, 3.8657, 3.873, 3.8803, 3.8876, 3.8949, \
3.90219, 3.9101, 3.91879, 3.92609, 3.93339, 3.94069, 3.94799, \
3.95529, 3.96259, 3.96988, 3.97718, 3.98448, 3.99178, 3.99908, \
4.00638, 4.01368, 4.02098, 4.02828, 4.03558, 4.04288, 4.05018, \
4.05748, 4.06477, 4.07207, 4.07937, 4.08667, 4.09397, 4.10127, \
4.10857, 4.11587, 4.12317, 4.13047, 4.13777, 4.14507, 4.15237, \
4.15998, 4.16867, 4.17597, 4.18327, 4.19057, 4.19786, 4.20516, \
4.21246, 4.21976, 4.22706, 4.23436, 4.24166, 4.24896, 4.25626, \
4.26356, 4.27086, 4.27816, 4.28546, 4.29275, 4.30005, 4.30735, \
4.31465, 4.32195, 4.32925, 4.33655, 4.34385, 4.35115, 4.3538, \
4.36249, 4.36979, 4.37709, 4.38439, 4.39168, 4.39898, 4.40628, \
4.41358, 4.42088, 4.42818, 4.43548, 4.44278, 4.45008, 4.45738, \
4.46468, 4.47198, 4.47928, 4.48657, 4.49387, 4.50117, 4.50847, \
4.51216, 4.52085, 4.52815, 4.53545, 4.54275, 4.55005, 4.55735, \
4.56465, 4.57194, 4.57924, 4.58654, 4.59384, 4.60114, 4.60844, \
4.61574, 4.62304, 4.63034, 4.63764, 4.64606, 4.65474, 4.66204, \
4.66934, 4.67664, 4.68394, 4.69124, 4.69854, 4.70584, 4.71314, \
4.72044, 4.72774, 4.73504, 4.74233, 4.74963, 4.75693, 4.76204, \
4.76423, 4.77073, 4.77791, 4.78521, 4.79263, 4.80019, 4.80787, \
4.8157, 4.82367, 4.83179, 4.84006, 4.8485, 4.8571, 4.86588, 4.87483, \
4.88398, 4.89332, 4.90287, 4.91263, 4.92261, 4.93284, 4.9433, \
4.95403, 4.96503, 4.97631, 4.98789, 4.9998, 5.01203, 5.02463, \
5.03759, 5.05096, 5.06475, 5.079, 5.09373, 5.10897, 5.12477, 5.14116, \
5.1582, 5.17594, 5.19443, 5.21374, 5.23395, 5.25515, 5.27744, 5.30093,
 5.32576, 5.35211, 5.38015, 5.41013, 5.44234, 5.47712]
    nws=len(wgrid)
    # if keep:
    #     wgrid=np.log(spec[0])
    #     nws=len(wgrid)
    # else:
    #     wgrid=np.log(wlo)+np.log(whi/wlo)*np.arange(0, nws)/(nws-1)
    newspec=griddata(np.log10(spec[0]), spec[1], wgrid, fill_value=1.e-36,method='linear')
    newspec=newspec[::-1]


    wgrid=10.**(np.array(wgrid))
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
    nu=np.copy(table[1][0,0])
    #Interpolating based on table (interpolating in teff/qg space)
    print table[0].shape,intens.shape,params.shape
    if method:
        grid=griddata(table[0], intens, params, method=method)
    else:
        grid=griddata(table[0], intens, params)

    #Angle space
    #Flag to return flux for all angles
    if mu>1:
        grid2=np.copy(grid)
    #Flag to return flux
    elif mu<0:
        grid2=np.copy(grid[:,:,:,-1])
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
    valid=np.copy(specs[1])

    specs=specs[0]
    specs=np.array(specs)
    nu=specs[0,0]
    f=specs[:, 1]


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


            if valid[i]:
                bb[j]+=tmpbb
                gb[j]+=tmpgb
                L[j] =(L[j] +rad*f[i,j])
                L2[j]=(L2[j]+rad*f[i,j])
            #Otherwise add bb flux to L but not L2
            # else:
            #     L[j]+=(tmpbb)

    return (nu, gb, bb, L, L2) 


##Calculates a composite disk spectrum given an file containing input radial parameters.
def disk_spec(f, table=[], tablef='tmpd', method='', logi=False, mu=-1, ymax=10.**46, ind=False):
    #ind=False
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
    print f
    disk_params=np.genfromtxt(f, skip_header=1)
    print disk_params
    specs=params_to_spec(disk_params[:, 1:4:2], table, method=method, logi=logi, mu=mu)

    #Disk parameters
    r=disk_params[:, 0]
    Teff=disk_params[:, 1]
    Qg=disk_params[:, 3]


    #Finding the total flux, for the given parameters
    totf=sum_spec(r, specs, Teff, Qg, mu, M=M)
    specs=specs[0]

    inc_factor=1
    if 0<mu<1:
        inc_factor=mu*4*np.pi
    nu=totf[0]
    totfg=inc_factor*totf[1]
    totfb=inc_factor*totf[2]
    totft=inc_factor*totf[3]
    totft2=inc_factor*totf[4]

    #r=r*(G*M*M_sun/c**2)  
    lr=np.log10(r)
    #For now assuming a logarithmically evenly spaced grid--for simplicity
    dlr=np.diff(lr)[0]
    dr=np.empty_like(r)
    for i in range(len(lr)):
        dr[i]=10.**(lr[i]+(dlr/2))-10.**(lr[i]-(dlr/2))
    outfile='sp_M'+'{0:.3e}'.format(M)+'_a'+'{0:.3f}'.format(a)+'_mu'+'{0:.3f}'.format(mu)
    np.savetxt(outfile, np.transpose([nu, totft]))

    #Plotting individual spectra
    if ind: 
        fig,ax=plt.subplots(nrows=2, ncols=1, figsize=(16,32),sharex=True, subplot_kw=dict(adjustable='datalim'))
        #plt.title(str(bin_params[0])+" "+str(bin_params[1])+" "+str(bin_params[2]))
        plt.xlabel(r"$\nu$ [hz]")

        #Plotting the composite disk spectrum
        ax[0].set_ylabel(r"$\nu L_{\nu}$ [ergs s$^{-1}$]")
        ax[0].set_xlim(10.**14, 3*10.**16)
        ax[0].set_ylim(10.**38, 10.**45)
        ax[0].set_xscale('log')
        ax[0].set_yscale('log')
        
        ax[0].plot(nu, nu*totfg)
        ax[0].plot(nu, nu*totfb)
        ax[0].plot(nu, nu*totft)
        ax[0].plot(nu, nu*totft2)

        #Plotting the contributions of individual annuli
        ax[1].set_ylim(10.**6, 10.**16)
        ax[1].set_xlim(10.**14, 3*10.**16)
        ax[1].set_xscale('log')
        ax[1].set_yscale('log')
        for i in range(len(specs)):
            ax[1].plot(nu, specs[i, 1]*nu)
        plt.savefig('ind.png')
        plt.close()
    return [[nu,totfb],[nu,totfg],[nu,totft], specs]


##Given a radial disk profile, calculates the spectra of each of the annuli, and then call the kerrtrans9 routine
def disk_spec_gr(f, table=[], tablef='tmpd', method='', logi=False, fobs=[1.e14, 1.e17], rmax=1.e4, mu=0.6, ymax=1e46):
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
    ax.set_ylim(10.**-5*ymax, ymax)
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
        help='file(s) containing list of disk profile files ',
        nargs='*')
    parser.add_argument('-c', '--colors',
        help='color to use for plotting ',
        nargs='*',
        default=['b','r','k'])
    parser.add_argument('-s', '--symbols',
        help='symbols used for plotting ',
        nargs='*',
        default=['-', '--', '--', '--'])
    parser.add_argument('-mu', '--mu',
        help='cosine of inslincation angle for the case of a disk ',
        type=float,
        default=-1)
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
    parser.add_argument('-ymax', '--ymax',
        help='maximum y for plotting ',
        type=float,
        default=10.**45)
    parser.add_argument('-ind', '--individual',
        help='specifies whether to plot individual annuli',
        action='store_true')
    parser.add_argument('-bb', '--bb',
        help='Whether we should overplot blackbody and graybody spectra',
        action='store_true')


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
    ymax=args.ymax
    ind=args.individual
    bb=args.bb
    cols=args.colors
    syms=args.symbols
    #cols=['k','0.5','r']
    #if args.colors:
    labels=['BB', 'GB', 'RT']

    if d:
        #Set-up for plotting spectra
        plt.close('all')
        #fig,ax=plt.subplots(len(d), figsize=(10,8*len(d)), sharex=True)
        gs1 = gridspec.GridSpec(1, len(d))
        fig=plt.figure(figsize=(12*len(d), 9))
        #fig.subplots_adjust(bottom=0.2)
        #fig.subplots_adjust(left=0.2)

        size=40
        mpl.rcParams['axes.labelsize']=size
        mpl.rcParams['xtick.labelsize']=size
        mpl.rcParams['ytick.labelsize']=size
        mpl.rcParams['legend.fontsize']=size
        mpl.rcParams['xtick.major.pad']=12
        mpl.rcParams['lines.linewidth']=2
        

        for k in range(len(d)):
            ax = fig.add_subplot(gs1[k])
            if k==0:
                ax.set_ylabel(r"$\nu \rm L_{\nu}$ [ergs s$^{-1}$]")
            else:
                plt.setp(ax.get_yticklabels(), visible=False)

            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlabel(r"$\nu$ [Hz]")
            ax.set_ylim(10.**-5*ymax, ymax)
            ax.set_xlim(10.**14, 3.*10**16)

            table=construct_table(tablef, logi=logi)
            #pdf_pages = PdfPages('composite.pdf')
            param_files=np.genfromtxt(d[k], dtype=str)
            param_files=np.atleast_1d(param_files)

            for i in range(len(param_files)):
                #print param_files
                if gr:
                    fig=disk_spec_gr(param_files[i], table=table, tablef=tablef, method=method,logi=logi,rmax=rmax,mu=mu,ymax=ymax)
                    fig.savefig('composite_'+str(i)+'_gr'+'.png')
                else:

                    spec=disk_spec(param_files[i], table=table, tablef=tablef, method=method,logi=logi,mu=mu,ymax=ymax,ind=ind)
                    #print spec[1][1],spec[0][1]
                    if bb:
                        for j in range(3):
                            ax.plot(spec[j][0],spec[j][0]*spec[j][1], syms[i%len(syms)], color=cols[j], label=labels[j])
                            if k==0:
                                ax.annotate('RT', xy=(3.28E15, 1.E43), xytext=(3.28E15, 4.E43), fontsize=size, arrowprops=dict(arrowstyle="->"), color='black')
                                ax.annotate('BB', xy=(2.9E15, 1.E42), xytext=(1.E15, 2.E42),fontsize=size, arrowprops=dict(arrowstyle="->"), color='blue')
                                ax.annotate('GB', xy=(5.E15, 1E41), xytext=(7.E15, 1.E41), fontsize=size, arrowprops=dict(arrowstyle="->"), color='red')  
                        # if i==0:
                        #     ax.legend(bbox_to_anchor=(0.3, 1))
                    else:
                        ax.plot(spec[2][0],spec[2][0]*spec[2][1], syms[i%len(syms)], color=cols[i%len(cols)])
                        np.savetxt('spec.txt', np.transpose([spec[2][0], spec[2][0]*spec[2][1]]))
                    # for j in range(3):
                    #     if bb or j==2:
                    #         ax.plot(spec[j][0],spec[j][0]*spec[j][1], syms[i%len(syms)], color=cols[j%len(cols)])

        fig.tight_layout()
        fig.savefig('composite.eps')
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
 
    else:
        parser.print_help()

    #pdf_pages.close()



if __name__ == '__main__':
    main()
    
    

    

