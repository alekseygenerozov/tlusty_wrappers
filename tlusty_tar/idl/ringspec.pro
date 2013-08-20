; ****************************************************************
; ****************************************************************
; PRO ringspec -- procedure for interpolating spectra for a specified
;     set of TEFF, QGRAV, AND M0 using a precomputed table of annuli 
;     spectra.  Called by agnspec.pro
;
pro ringspec,irad,teff,dm0,qgrav,int0,out

common model,tl,ml,ql,wls,int,pol,iconv
;
out=1
coeff=dblarr(8)
nt=n_elements(tl)
;maxt=tl[nt-1]
;mint=tl[0]
it=findgen(nt)
ii=interpol(it,tl,teff) 
i=fix(ii)+1
i=i>1<(nt-1)
i1=i-1
a=(teff-tl(i1))/(tl(i)-tl(i1))
a1=1.-a
;
nm=n_elements(ml)
;maxm=ml[nm-1]
;minm=ml[0]
im=findgen(nm)
ii=interpol(im,ml,dm0) 
j=fix(ii)+1
j=j>1<(nm-1)
j1=j-1
b=(dm0-ml(j1))/(ml(j)-ml(j1))
b1=1.-b
;
nq=n_elements(ql)
;maxq=ql[nq-1]
;minq=ql[0]
iq=findgen(nq)
ii=interpol(iq,ql,qgrav) 
k=fix(ii)+1
k=k>1<(nq-1)
k1=k-1
c=(qgrav-ql(k1))/(ql(k)-ql(k1))
c1=1.-c

if(iconv[i,j,k] eq 0) then begin
    print,tl[i],ml[j],ql[k]
    out=0
;    iconv[i,j,k]=1
endif
if(iconv[i1,j,k] eq 0) then begin
    print,tl[i1],ml[j],ql[k]
    out=0
;    iconv[i1,j,k]=1
endif
if(iconv[i,j1,k] eq 0) then begin
    print,tl[i],ml[j1],ql[k]
    out=0
;    iconv[i,j1,k]=1
endif
if(iconv[i,j,k1] eq 0) then begin
    print,tl[i],ml[j],ql[k1]
    out=0
;    iconv[i,j,k1]=1
endif
if(iconv[i,j1,k1] eq 0) then begin
    print,tl[i],ml[j1],ql[k1]
    out=0
;    iconv[i,j1,k1]=1
endif
if(iconv[i1,j,k1] eq 0) then begin
    print,tl[i1],ml[j],ql[k1]
    out=0
;    iconv[i1,j,k1]=1
endif
if(iconv[i1,j1,k] eq 0) then begin
    print,tl[i1],ml[j1],ql[k]
    out=0
;    iconv[i1,j1,k]=1
endif
if(iconv[i1,j1,k1] eq 0) then begin
    print,tl[i1],ml[j1],ql[k1]
    out=0
;    iconv[i1,j1,k1]=1
endif

in0=b*(a*int(i,j,k,*,*)+a1*int(i1,j,k,*,*)) + $
     b1*(a*int(i,j1,k,*,*)+a1*int(i1,j1,k,*,*)) 
in1=b*(a*int(i,j,k1,*,*)+a1*int(i1,j,k1,*,*)) + $
     b1*(a*int(i,j1,k1,*,*)+a1*int(i1,j1,k1,*,*)) 
int0=c*in0+c1*in1

return
end
