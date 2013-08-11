; ****************************************************************
; ****************************************************************
; PRO: readgrid -- procedure which generates/reads an idl .save file
;      containing at set of annuli spectra tabulated by TEFF, QGRAV,
;      and M0.  Used with agnspec.pro to initiall model common block 
;
pro readgrid,save=save,name=name,$
             wlo=wlo,whi=whi,nws=nws,fc=fc
;
common model,tl,ml,ql,wls,int,pol,iconv
;
if n_elements(name) eq 0 then name='agngrid.save'
;
if n_elements(save) gt 0 then begin
;
if n_elements(wlo) eq 0 then wlo=0.06
if n_elements(whi) eq 0 then whi=12400
if n_elements(nws) eq 0 then nws=300
if n_elements(fc) eq 0 then fc=2.0

c=2.9979d10
h=6.6262d-27
kb=1.3807d-16

wls=alog(wlo)+alog(whi/wlo)/(nws-1)*findgen(nws)
freqs=dblarr(nws)
freqs=c*1.0d8/exp(wls)
na=10
x0=fltarr(2*na+2)
na1=na-1
na2=na+2
nw0=500
wl=fltarr(nw0)
int0=fltarr(nw0,na)
pol0=int0


tl=5.0+findgen(26)*0.1
ml=[2.0,2.25,2.5,2.75,3.0,4.0,5.0,6.0,7.0]
ql=-5.0+findgen(17)
tlab=['50','51','52','53','54','55','56','57','58','59','60','61','62', $
      '63','64','65','66','67','68','69','70','71','72','73','74','75']
;tlab=['60','61','62', $
;      '63','64','65','66','67','68','69','70','71','72','73','74']
mlab=['20','22','25','27','30','40','50','60','70']
qlab=['-50','-40','-30','-20','-10','00','10','20','30','40','50','60','70','80','90','100','110']


int=fltarr(n_elements(tl),n_elements(ml),n_elements(ql),nws,na)
pol=int
iconv=fltarr(n_elements(tl),n_elements(ml),n_elements(ql))
;

spawn,'ls *.14 > tmpd'
get_lun,lun1
openr,lun1,'tmpd'
iex=0
ii=0
a=''
modl=strarr(1000)
while not eof(lun1) do begin
  readf,lun1,a
  modl(ii)=strtrim(a,2)
  ii=ii+1
endwhile
modl=modl(0:ii-1)
;

;  This is the part that produces the save file.  After the models are
;  read in from the fort.10 files they are interpolated onto the wls
;  grid used by agnspec.

for i=0,n_elements(tl)-1 do begin
;  td=tl(i)*10.
;  if td eq fix(td)/2*2. then dire='' else dire='add/'
;  if td gt 51. then dire='add/'
    dire='./'
    for j=0,n_elements(ml)-1 do begin
        for k=0,n_elements(ql)-1 do begin
            a='t'+tlab(i)+'m'+mlab(j)+'q'+qlab(k)
            iex=max(a+'.14' eq modl)
;      print,a,iex
            if iex eq 1 then begin
                iconv(i,j,k)=1
            endif
            print,a,iconv(i,j,k)
            a=dire+a
            if iconv(i,j,k) eq 1 then begin
                get_lun,lun1
                openr,lun1,a+'.14'
                n=0L
                w0=0.
                wl=fltarr(nw0)
                int0=fltarr(nw0,na)
                pol0=int0
                while not eof(lun1) do begin
                    readf,lun1,x0
                    if x0(0) lt w0 then goto,endread 
                    wl(n)=x0(0)
                    for m=0,na1 do begin
                        int0(n,m)=x0(2*m+2)
                        pol0(n,m)=x0(2*m+3)
                    endfor
                    n=n+1
                    w0=x0(0)
                endwhile
                endread: nw=n
                close,lun1
                free_lun,lun1
                wl=alog(wl(0:nw-1))
                int0=(int0(0:nw-1,0:na1)>1.e-36)
                pol0=pol0(0:nw-1,0:na1)
                for m=0,na1 do begin
                    in0=int0(*,m)
                    ins=interpol(in0,wl,wls)
                    for n=0,nws-1 do begin
                        if(wls(n) lt wl(0)) then ins(n)=(1.e-36)
                        if(wls(n) gt wl(nw-1)) then ins(n)=(1.e-36)
                        if(ins(n) lt 0.0) then ins(n)=(1.e-36)
                    endfor
                                ;     print,'in0:',in0
                                ;     print,'ins:',ins
                    pl0=pol0(*,m)
                    pls=interpol(pl0,wl,wls)
                    for n=0,nws-1 do begin
                        int(i,j,k,n,m)=alog(h/fc/kb*freqs(n)/alog(1.0+2.0*h*(freqs(n))^3/fc^4/c^2/ins(n)) )
                        if(n eq 100 and m eq 4) then print,int(i,j,k,n,m),ins(n)
                        pol(i,j,k,n,m)=pls(n)
                    endfor
                endfor
            endif else begin
                for m=0,na1 do begin 
                    for n=0,nws-1 do begin
                        int(i,j,k,n,m)=alog(1.e-36)
                        pol(i,j,k,n,m)=alog(1.e-36)
                    endfor
                endfor
            endelse
        endfor
    endfor
endfor
;
;help
save,tl,ml,ql,tlab,mlab,qlab,wls,int,pol,iconv,file=name
;
endif else restore,file=name
;
return
end


