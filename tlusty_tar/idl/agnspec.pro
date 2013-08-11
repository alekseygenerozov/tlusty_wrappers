;
; PURPOSE:  contains IDL procedure which generates a full relativistic
;           disk spectrum from TLUSTY annuli spectra.   
;
; ****************************************************************
; ****************************************************************
; PRO: agnspec -- procedure which constructs full relativistic
;                 accretion disk spectrum from input annuli spectra.
; used with:
; readgrid.pro -- must be called first to read in data model common block
; ringspec.pro -- procedure which implements interpolation method
; gridgen.pro -- wrapper for agnspec
;
pro agnspec,mass=mass,mdot=mdot,angm=angm,alpha=alpha,mu=mu,    $
            rout=rout,tmin=tmin,deltar=deltar,nre=nre,          $
            rcut=rcut,tlim=tlim,dirprg=dirprg,                  $
            restore=restore,file=file,nocomp=nocomp,            $
            frmax=frmax,frmin=frmin,nfobs=nfobs,                $
            plmod=plmod,plspec=plspec,opl=opl,_extra=e,         $
            xh=xh,fc=fc,lumin=lumin,kerrout=kerrout,inc=inc
;
common model,tl,ml,ql,wls,int,pol,iconv
;
; initialization; default values
;
if n_elements(file) eq 0 then file='agnspec.save'
if n_elements(restore) gt 0 then restore,file
if n_elements(dirprg) eq 0 then dirprg='./'
;
if n_elements(mass) eq 0 then mass=10.0
if n_elements(mdot) eq 0 then mdot=3.9e-8
if n_elements(angm) eq 0 then angm=0.0
if n_elements(alpha) eq 0 then alpha=0.1
if n_elements(mu) eq 0 then mu=0.34
if n_elements(rout) eq 0 then rout=0.
if n_elements(tmin) eq 0 then tmin=1250000.
if n_elements(deltar) eq 0 then deltar=0.1
if n_elements(nre) eq 0 then nre=0
if n_elements(tlim) eq 0 then tlim=1000.
if n_elements(frmax) eq 0 then frmax=1.e19
if n_elements(plspec) eq 0 then plspec=1
if n_elements(xh) eq 0 then xh=0.7
if n_elements(fc) eq 0 then fc=1.0
if n_elements(kerrout) eq 0 then kerrout='sp.out'
if n_elements(inc) eq 0 then inc=[mu]

; physical constants
c=2.9979d10
G=6.67259d-8
h=6.6262d-27
kb=1.3807d-16
me=9.1095d-28
msun=1.9891d33
mp=1.6726d-24
sigmat=6.65248d-25

if(angm ge 0.) then begin
    signa=1.
endif else begin
    signa=-1.
endelse

ledd=4.0*!pi*G*msun*mass*c*mp/sigmat
a2=angm*angm
z1=1.0+(1.0-a2)^(1.0/3.0)*((1.0+angm)^(1.0/3.0)+(1.0-angm)^(1.0/3.0))
z2=sqrt(3.0*a2+z1*z1)
rms=3.0+z2-signa*sqrt((3.0-z1)*(3.0+z1+2.0*z2))
eff=1.0-(1.0-2.0/rms+angm/rms^1.5)/sqrt(1.0-3.0/rms+2.0*angm/rms^1.5)
;print,eff,ledd*lumin
; If desired, you Eddington luminosity for input rather than accretion
; rate
if n_elements(lumin) ne 0 then begin
    mdot=lumin*ledd/c/c/eff/6.3029D25
;    print,mdot
endif else begin
    lumin=mdot*c*c*eff*6.3029D25
endelse

nfs=n_elements(wls)
freqs=dblarr(nfs)
freqs=c*1.0d8/exp(wls)
wfs=dblarr(nfs)
wfs[0]=-freqs[1]+freqs[0]
for i=1,nfs-2 do begin
    wfs[i]=0.5*(-freqs[i+1]+freqs[i-1])
endfor
wfs[nfs-1]=-freqs[nfs-1]+freqs[nfs-2]
cosarr=[0.0130467357,0.0674683167,0.160295216,0.283302303,0.425562831, $
        0.574437169,0.716697697,0.839704784,0.932531683,0.986953264]
warr=[0.0333356722,0.0747256746,0.109543181,0.13463336,0.147762112, $
      0.147762112,0.13463336,0.109543181,0.0747256746,0.0333356722]

nfgrid=300
frgrid0=1.e19
frgrid1=1.e14
nretype=4
gr=1
;
; -----------------------------------------------
; 1. parameters and spectra of individual annuli:
; -----------------------------------------------
;
; -----------------------------------------------
; 1a. set up input data for dispar3.f
; -----------------------------------------------
;

get_lun,lun1
openw,lun1,'tmp.5'
printf,lun1,mass,mdot,angm,alpha
printf,lun1,rout,tmin,deltar,nre
printf,lun1,xh
close,lun1
free_lun,lun1
;
; -----------------------------------------------
; 1b. run dispar3 and read results
; -----------------------------------------------
;
a=dirprg+'dispar3 <tmp.5 > tmp.6'
spawn,a
;
get_lun,lun1
openr,lun1,'tmp.6'
ii=0
r=fltarr(300) & tr=r & mr=r & qr=r
while not eof(lun1) do begin
   readf,lun1,r0,t0,m0,q0,teff,dm,qgrav
   r(ii)=r0 & tr(ii)=t0 & mr(ii)=m0 & qr(ii)=q0
;   print,ii,r(ii),tr(ii),mr(ii),qr(ii),teff,dm,qgrav,$
;         format='(i3,4f10.2,f10.1,2e10.2)'
   ii=ii+1
endwhile
close,lun1
free_lun,lun1
r=r(0:ii-1)
tr=tr(0:ii-1)
mr=mr(0:ii-1)
qr=qr(0:ii-1)
nr=n_elements(r)
;if (nr gt 50) then print,nr,angm,mass,lumin
if n_elements(rcut) eq 0 then rcut=r(nr-1)*(tlim/10.^tr(nr-1))^(-1.333)
;
; -----------------------------------------------
; 1c spectra for individual annuli
; -----------------------------------------------
;

tmpin=fltarr(20,nfs)
totr=fltarr(n_elements(wls))
get_lun,lun1
get_lun,l2
get_lun,l3
openw,lun1,'emrad.in'
openu,l2,'nocomp1.out',/append
openu,l3,'nocomp2.out',/append
for i=0,nr-1 do begin
;   print,r[i],tr[i],mr[i],qr[i]
   lnn=i-(i/5)*5
   nq=n_elements(ql)
   nm=n_elements(ml)
   nt=n_elements(tl)
   flag=0
   if(qr[i] lt ql[0]) then flag=1
   if(qr[i] gt ql[nq-1]) then flag=1
   if(mr[i] lt ml[0]) then flag=1
   if(mr[i] gt ml[nm-1]) then flag=1
   if(tr[i] lt tl[0]) then flag=1
   if(tr[i] gt tl[nt-1]) then flag=1
   if(tr[i] lt alog10(tmin)) then begin
       print,tr[i],' less than tmin: ',tmin
       print,r[i],tr[i],mr[i],qr[i]
   endif
   if(flag eq 1) then begin
       printf,l2,'parameters not within grid:'
       printf,l2,r[i],tr[i],mr[i],qr[i]
       printf,l2,alpha,mass,angm,lumin
       close,lun1
       close,l2
       close,l3
       free_lun,lun1
       free_lun,l2
       free_lun,l3
       return
   endif
   ringspec,i,tr(i),mr(i),qr(i),int0,out
   if(out eq 0) then begin
   ;    print,r[i],tr[i],mr[i],qr[i]
       printf,l3,'parameters not within converged grid:'
       printf,l3,r[i],tr[i],mr[i],qr[i]
       printf,l3,alpha,mass,angm,lumin
       close,lun1
       close,l2
       close,l3
       free_lun,lun1
       free_lun,l2
       free_lun,l3
       return
   endif
   if n_elements(nocomp) eq 0 then begin
       printf,lun1,r(i),n_elements(wls)
       fluxtot=0.0
       for j=0,n_elements(wls)-1 do begin
           for k=0,9 do begin
               tbright=exp(int0(0,0,0,j,k))
               tmpin(2*k,j)=2.0*h*(freqs(j))^3/fc^4/c^2/(exp(h*freqs(j)/(kb*tbright*fc))-1.0)
               fluxtot=fluxtot+tmpin(2*k,j)*wfs[j]*warr[k]*2.0*!pi*cosarr[k]
           endfor
           for k=1,9 do tmpin(2*k+1,j)=0.
       endfor
       norm=fluxtot/5.67051d-5/10^(4.0*tr[i])
       for j=0,n_elements(wls)-1 do begin
           if(norm gt 0) then begin
               for k=0,9 do begin
                   tmpin(2*k,j)=tmpin(2*k,j)/norm
               endfor
           endif
           printf,lun1,exp(wls(j)),tmpin(0,j)
           printf,lun1,tmpin[*,j]
       endfor
   endif
endfor
;if n_elements(plmod) gt 0 then oplot,3.e18/exp(wls),totr,thick=2,$
;                                     _extra=e
close,lun1
close,l2
close,l3
free_lun,lun1
free_lun,l2
free_lun,l3
;
; -----------------------------------------------
; 2. integration over annuli to get the total spectrum
; -----------------------------------------------
;
; -----------------------------------------------
; 2a. set up input data for kerrtrans9: std input
; -----------------------------------------------
;

if n_elements(nocomp) gt 0 then return

; Loop over inclination array running kerrtrans for each inclination
;
kerrout1=kerrout+'.out'
get_lun,l2
openw,l2,kerrout1
close,l2
free_lun,l2

for iinc=0,n_elements(inc)-1 do begin

get_lun,lun1
openw,lun1,'tmp.in'
printf,lun1,-fix(nfgrid)
printf,lun1,frgrid1,frgrid0
;
if frmax gt 0 then  begin
  if n_elements(frmin) eq 0 then frmin=2.418e16  
  if n_elements(nfobs) eq 0 then nfobs=350
  printf,lun1,-fix(nfobs)
  printf,lun1,frmin,frmax
endif else begin
  get_lun,lun2
  a='wc -l freq.in > tmpf'
  spawn,a
  get_lun,lun2
  openr,lun2,'tmpf'
  readf,lun2,nko
  close,lun2
  free_lun,lun2
  printf,lun1,fix(nko)
endelse
;
printf,lun1,inc[iinc]
printf,lun1,mass,mdot,angm
printf,lun1,rcut,tmin,nr,nretype,gr
close,lun1
free_lun,lun1
;
; -----------------------------------------------
; 2b. run kerrtrans9
; -----------------------------------------------
;
li=strtrim(string(inc[iinc],FORMAT='(f4.2)'),2)
kerrout2=kerrout+'i'+li+'.out'
a='nice -n 19 '+dirprg+'kerrtrans9  <tmp.in > '+kerrout2
spawn,a
com='cat '+kerrout1+' '+kerrout2+' > kerrtemp'
spawn,com
com='mv kerrtemp '+kerrout1
spawn,com
com='rm '+kerrout2
spawn,com
;if n_elements(plspec) gt 0 then plt,'sp.out',nc=5,$
;                           /xlog,opl=opl,_extra=e
;
endfor

return
end

; **********************************************************************
; **********************************************************************

pro specgen,inc=inc,fc=fc,spin=spin,l=l,m=m,alpha=alpha,infile=infile

common model,tl,ml,ql,wls,int,pol,iconv

if n_elements(fc) eq 0 then fc=2.0
if n_elements(inc) eq 0 then inc=[0.5]
if n_elements(l) eq 0 then l=-0.5
if n_elements(m) eq 0 then m=1.0
if n_elements(alpha) eq 0 then alpha=-2
if n_elements(spin) eq 0 then a=0. else a=spin

if ((n_elements(wls) eq 0) or (n_elements(infile) ne 0)) then begin
    if n_elements(infile) eq 0 then infile='agngrid.save'
    readgrid,name=infile
endif
tmin=100000.0
rcut=10000.0

lp=strtrim(string(-alpha,FORMAT='(f4.2)'),2)
lm=strtrim(string(m,FORMAT='(f5.3)'),2)
la=strtrim(string(a,FORMAT='(f7.5)'),2)
if(l gt 0.0) then begin
    ll=strtrim(string(l,FORMAT='(f4.2)'),2)
endif else begin
    ll=strtrim(string(l,FORMAT='(f5.2)'),2)    
endelse
outfile='p'+lp+'m'+lm+'a'+la+'l'+ll

agnspec,rcut=rcut,inc=inc,mass=10^(m),angm=a,fc=fc, $
        lumin=10^(l),alpha=10^(alpha),kerrout=outfile, $
        tmin=tmin,deltar=0.1

return

end


