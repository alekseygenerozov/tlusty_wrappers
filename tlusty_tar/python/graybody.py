import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


kes=0.4
h = 6.63*10.**-27;
kb=1.38*10**-16
c=3.*10**10
gamma=4./3.
sigma=5.67*10**-5
kr=1.6*10**24

def xi(t, nu):
    return h*nu/(kb*t)

def fnu(t, nu):
    xi1=xi(t, nu)
    return xi1**-3/(1.-np.exp(-xi1))

def knu(es, t, nu):
    f=fnu(t, nu)
    a=(1-es**-1)**-2

    #Finding knu by solving cubic equation as described in Tanaka & Menou 2010
    roots=np.roots([1,1,0,-a*f])

    #Select only positive real roots from cubic 
    good_roots=np.empty(0)
    for r in roots:
        if r > 0 and np.allclose(0.,np.imag(r)):
            good_roots=np.append(good_roots, r)


    #Check if there is more than on positive root. By Descartes' rule of signs this should not occur.
    if len(good_roots)>1:
        raise Exception('More than one positive, real root in cubic!')
    #Return the positive, real root
    good_roots=np.real_if_close(good_roots)
    return good_roots[0]


def eps(es, t, nu):
    try:
        return (1+knu(es, t, nu)**-1)**-1
    except Exception as e:
        print e.args[0]
        return 0
   
def ks(tp, qg):
    np.seterr(all='raise')
    norm=np.sqrt(140.*c*kr/(gamma*sigma*kes))
    return norm*np.sqrt(qg)*tp**(-15./4.)

def es(tp, qg):
    # kes=0.4
    return ks(tp, qg)/(ks(tp, qg)+kes)

def Chi(es):
    return 0.873*es**(-1./6.)/(1-0.127*es**(5./6.))*1/(1+(es**-1-1)**(2./3.))

def F(tp, teff, qg):
    return teff-tp*(Chi(es(tp, qg)))**0.25

  
#Standard blackbody
def bb(nu, t):
    try:
        return 2*h*nu**3/c**2/(np.exp(h*nu/(kb*t))-1)
    except:
        return 0

#Calculates the graybody flux for a given effective temperature and gravity parameter
def gb(nu, teff, qg):
    tp=fsolve(F, teff, args=(teff, qg))[0]
    es1=es(tp, qg)
    eps1=eps(es1, tp, nu)
    try:
        return 2*np.pi*eps1**0.5/(1+eps1**0.5)*bb( nu, tp)
    except:
        return 0




# nus=np.logspace(14, 16, num=300)
# #gb1=np.zeros(300)
# bb1=np.zeros(300)
# gb1=np.zeros(300)
# for i in range(len(nus)):
#      gb1[i]=gb(nus[i], 10.**4.5, 10.**-12)
#      bb1[i]=np.pi*bb(nus[i],10.**4.5)



# plt.plot(nus, bb1)
# plt.plot(nus, gb1)
# plt.loglog()
# plt.show()





