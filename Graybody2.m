(* ::Package:: *)

(************************************************************************)
(* This file was generated automatically by the Mathematica front end.  *)
(* It contains Initialization cells from a Notebook file, which         *)
(* typically will have the same name as this file except ending in      *)
(* ".nb" instead of ".m".                                               *)
(*                                                                      *)
(* This file is intended to be loaded into the Mathematica kernel using *)
(* the package loading commands Get or Needs.  Doing so is equivalent   *)
(* to using the Evaluate Initialization Cells menu command in the front *)
(* end.                                                                 *)
(*                                                                      *)
(* DO NOT EDIT THIS FILE.  This entire file is regenerated              *)
(* automatically each time the parent Notebook file is saved in the     *)
(* Mathematica front end.  Any changes you make to this file will be    *)
(* overwritten.                                                         *)
(************************************************************************)



Clear["Global`*"]
Needs["SigFig`"]
Needs["PlotLegends`"]
Needs["CustomTicks`"]


K\[Nu]=Solve[y^2 (1+y)==(x)^2 f, y]//StandardForm
(*Physical constants in cgs units*)
G=6.67 10^-8;  (*Newton's constant in cgs*)
c=3 10^10 ; (*Speed of light in cgs*)
\[Kappa]r=1.6 10^24;(*Rosseland mean opacity at solar metallicity*)
kb=1.38 10^-16;
\[Mu]0=0.615;
mp=1.67 10^-24;
me=9.11 10^-28;
(*Stefan-Boltzmann constant in cgs*)
\[Sigma]=5.67 10^-5;
h=6.63 10^-27;
Msun=2 10^33;
(*keV in cgs units*)
eV=1.6 10^-12;
keV=10^3 eV;
R=13.6 eV;
K=2.815 10^29;
nmax=8;
(*Table of ionization energies for various species*)
ions={{"H11",2.1785304`*^-11},{"H12",5.446326`*^-12},{"H13",2.4205893`*^-12},{"H14",1.3615815`*^-12},{"H15",1,8.7141216`*^-13},{"H16",6.0514734`*^-13},{"H17",4.4459804`*^-13},{"H18",3.4039538`*^-13},{"He1",3.9389425`*^-11},{"He21",8.7176979`*^-11},{"He22",2.1794245`*^-11},{"He23",2,9.686331`*^-12},{"He24",2,5.4485612`*^-12},{"c1",1.8042644`*^-11},{"c2",3.9063653`*^-11},{"c3",7.6719252`*^-11},{"c4",1.0331629`*^-10},{"c5",6.2816386`*^-10},{"c6",7.8427095`*^-10},{"n1",2.3296135`*^-11},{"n2",4.7423235`*^-11},{"n3",7.6016697`*^-11},{"n4",1.2411747`*^-10},{"n5",1.5682464`*^-10},{"n6",8.844911`*^-10},{"n7",1.0674799`*^-9},{"o1",2.1812346`*^-11},{"o2",5.6260111`*^-11},{"o3",8.8010408`*^-11},{"o4",1.2402172`*^-10},{"o5",1.8247366`*^-10},{"o6",2.2124336`*^-10},{"o7",1.1843989`*^-9},{"o8",1.3942595`*^-9},{"ne1",3.6128602`*^-11},{"ne2",6.8628627`*^-11},{"ne3",1.0631106`*^-10},{"ne4",1.6270685`*^-10},{"ne5",2.1145801`*^-10},{"ne6",2.6459267`*^-10},{"ne7",3.472636`*^-10},{"ne8",4.0057914`*^-10},{"ne9",2.0034556`*^-9},{"ne10",2.1785304`*^-9},{"mg1",1.2810266`*^-11},{"mg2",2.5189736`*^-11},{"mg3",1.3427176`*^-10},{"mg4",1.830216`*^-10},{"mg5",2.3667571`*^-10},{"mg6",3.1247191`*^-10},{"mg7",3.7686545`*^-10},{"mg8",4.4549541`*^-10},{"mg9",5.4945902`*^-10},{"mg10",6.1576339`*^-10},{"mg11",2.9517446`*^-9},{"mg12",3.1370838`*^-9},{"si1",1.3059588`*^-11},{"si2",2.6187021`*^-11},{"si3",5.3657819`*^-11},{"si4",7.2319749`*^-11},{"si5",2.6717798`*^-10},{"si6",3.2851514`*^-10},{"si7",3.9495665`*^-10},{"si8",4.8571081`*^-10},{"si9",5.6250151`*^-10},{"si10",6.4312579`*^-10},{"si11",7.6269932`*^-10},{"si12",8.3869495`*^-10},{"si13",3.9054135`*^-9},{"si14",4.2699196`*^-9},{"s1",1.7356621`*^-11},{"s2",3.9092763`*^-11},{"s3",5.8348346`*^-11},{"s4",7.925344`*^-11},{"s5",1.2176528`*^-10},{"s6",1.4751898`*^-10},{"s7",4.7067136`*^-10},{"s8",5.4991619`*^-10},{"s9",6.3514327`*^-10},{"s10",7.4905721`*^-10},{"s11",8.4571146`*^-10},{"s12",9.4601642`*^-10},{"s13",1.0917531`*^-9},{"s14",1.18475`*^-9},{"s15",5.4012554`*^-9},{"s16",5.5770378`*^-9},{"ar1",2.6403215`*^-11},{"ar2",4.6289821`*^-11},{"ar3",6.8256269`*^-11},{"ar4",1.0020359`*^-10},{"ar5",1.2569094`*^-10},{"ar6",1.5247427`*^-10},{"ar7",2.0828634`*^-10},{"ar8",2.4034762`*^-10},{"ar9",7.0775984`*^-10},{"ar10",8.019825`*^-10},{"ar11",9.0295664`*^-10},{"ar12",1.0358132`*^-9},{"ar13",1.1494753`*^-9},{"ar14",1.2661522`*^-9},{"ar15",1.4320572`*^-9},{"ar16",1.5380403`*^-9},{"ar17",6.9040077`*^-9},{"ar18",7.0584385`*^-9},{"ca1",1.0241786`*^-11},{"ca2",1.9889521`*^-11},{"ca3",8.5297974`*^-11},{"ca4",1.1270874`*^-10},{"ca5",1.4158178`*^-10},{"ca6",1.8225303`*^-10},{"ca7",2.1312037`*^-10},{"ca8",2.4668765`*^-10},{"ca9",3.158662`*^-10},{"ca10",3.539634`*^-10},{"ca11",9.9165355`*^-10},{"ca12",1.1011217`*^-9},{"ca13",1.2174407`*^-9},{"ca14",1.3699091`*^-9},{"ca15",1.4986975`*^-9},{"ca16",1.6326737`*^-9},{"ca17",1.8216955`*^-9},{"ca18",1.9396842`*^-9},{"ca19",8.5925419`*^-9},{"ca20",8.7141216`*^-9},{"fe1",1.3239406`*^-11},{"fe2",2.7120436`*^-11},{"fe3",5.1352375`*^-11},{"fe4",9.1812264`*^-11},{"fe5",1.2567041`*^-10},{"fe6",1.6596797`*^-10},{"fe7",2.0938089`*^-10},{"fe8",2.5308136`*^-10},{"fe9",3.9134305`*^-10},{"fe10",4.3911893`*^-10},{"fe11",4.8627134`*^-10},{"fe12",5.54195`*^-10},{"fe13",6.0487885`*^-10},{"fe14",6.5701636`*^-10},{"fe15",7.6565434`*^-10},{"fe16",8.1983187`*^-10},{"fe17",2.1145801`*^-9},{"fe18",2.281797`*^-9},{"fe19",2.4614767`*^-9},{"fe20",2.6396987`*^-9},{"fe21",2.8291378`*^-9},{"fe22",3.0140053`*^-9},{"fe23",3.2813417`*^-9},{"fe24",3.4273699`*^-9},{"fe25",1.4790261`*^-8},{"fe26",1.4726866`*^-8}};

ionshhe={{"H11",2.1785304`*^-11},{"H12",5.446326`*^-12},{"H13",2.4205893`*^-12},{"H14",1.3615815`*^-12},{"H15",1,8.7141216`*^-13},{"H16",6.0514734`*^-13},{"H17",4.4459804`*^-13},{"H18",3.4039538`*^-13},{"He1",3.9389425`*^-11},{"He21",8.7176979`*^-11},{"He22",2.1794245`*^-11},{"He23",2,9.686331`*^-12},{"He24",2,5.4485612`*^-12}}

mydir="/home/aleksey/First_Year_Project/tlusty_tar/examples/aleksey_tlusty_runs"

Options[Edge]={metals->False}
Edge[\[Nu]_, OptionsPattern[]]:=Module[{en, pos,close, ion,m, ions2},
m=OptionValue[metals];
If[m, ions2=ions, ions2=ionshhe];
(*Calculate energy in cgs units*)
en=h \[Nu];
close=Nearest[ions2[[All,2]], en][[1]];
pos=Position[ions2[[All, 2]],x_/; Abs[x-close]/close<10^-6];
ion=Extract[ions2, pos]//Flatten;
{ion[[1]], ion[[2]]/h}

]

ColGrad[n_, grad_:"Rainbow"]:=Module[{colorfunc, colors, colorscale},
colorscale=Range[0, 1, 1/(n-1)];
colorfunc=ColorData[grad][#]&;
colors=colorfunc/@colorscale
]

B[Teff_, fcol_:1]:=1/fcol^4 (2 h \[Nu]1^3/c^2)/(E^((h \[Nu]1)/( kb  fcol Teff))-1)
B2[Teff_, \[Xi]_]:=(kb Teff)^3/(h c)^2 (2\[Xi]^3)/(E^\[Xi]-1)
\[Nu]peak[Teff_]:=2.82 kb Teff/h


u[n_, T_]:=(R)(1/n^2-1)/(kb T)
g[n_]:=2 n^2
nstar[\[Nu]_]:=(R/(h \[Nu]))^(1/2)//Ceiling
\[Alpha][\[Nu]_, T_]:=Sum[K/(n^5 \[Nu]^3) g[n] E^(u[n, T]-u[1, T]), {n,nstar[\[Nu]], nmax}]
(*Saha equation*)
fh[n_, T_]:=(1/n ((2\[Pi] me kb T)/h^2)^(3/2) E^(-R/(kb T))+1)^-1

Graybody3[T_, Q_, Dmtot_]:=Module[{\[Rho], H, \[Kappa]es, n, n2, \[Kappa], \[Kappa]r, \[Kappa]\[Nu], \[Nu]}, 
\[Kappa]es=0.4;
H=0.4/c \[Sigma] T^4/Q;
\[Rho]=Dmtot/H;
n=\[Rho] /mp;
n2=fh[ n, T]/mp;
\[Nu]=Range[Log[10,0.1\[Nu]peak[T]],Log[10, 10 \[Nu]peak[T]], 0.02];
\[Nu]=10^\[Nu];
\[Kappa]=n2 \[Alpha][#, T]&/@\[Nu];
\[Kappa]\[Nu]=Transpose[{\[Nu], \[Kappa]}];
\[Kappa]\[Nu]=DeleteCases[\[Kappa]\[Nu], x_/;x[[2]]==0];


\[Kappa]r= Rosseland[\[Kappa]\[Nu][[All, 1]], \[Kappa]\[Nu][[All, 2]], T];

Transpose[{\[Nu], \[Nu] Sqrt[\[Kappa] (\[Kappa]+\[Kappa]es)]/(\[Kappa]r+\[Kappa]es) (B[T]/.\[Nu]1->\[Nu])}]
(*\[Kappa]=With[{nc=n2},Function[\[Nu], nc (\[Alpha][\[Nu], T])]
];
\[Kappa]r=\[Kappa]ross[\[Kappa], T];
With[{\[Kappa]c=\[Kappa], \[Kappa]rc=\[Kappa]r},Function[\[Nu], Sqrt[\[Kappa]c[\[Nu]](\[Kappa]c[\[Nu]]+0.4)]/(\[Kappa]rc+0.4)]]*)

]
Rosseland[\[Nu]_, \[Kappa]_, T_]:=Module[{\[Delta]\[Nu], dB\[Nu]dT, int, \[Kappa]2},
dB\[Nu]dT=(B'[T]/.\[Nu]1->\[Nu])[[2;;]];

\[Kappa]2=\[Kappa][[2;;]];
\[Delta]\[Nu]=\[Nu]//Differences;
(*Calculating the integral in Rosseland mean opacity by right hand summation*)

int=(1/\[Kappa]2)dB\[Nu]dT \[Delta]\[Nu]//Total;
int= int/(4 \[Sigma] T^3/\[Pi]);
int=1/int
]

(*Calculates graybody atmosphere using Taka's formalism for a particular Teff and gravity parameter. Free is a flag to switch between the free-free and bound-free opacities*)
Graybody[Teff_, Qg_, Free_:False, prec_:10^-8]:=Module[{\[Kappa]es, \[Kappa]s, \[Kappa]sb, \[Kappa]sf, \[Epsilon]s,  \[Xi], f, x  , K, K\[Nu], \[Epsilon], Tp, \[CapitalXi], \[Rho]p, pr, pg, r, error, flux1, flux2},
\[Kappa]es=0.4;
(*Note that the below expression for \[Kappa]s assumes radiation pressure dominance. We also assume solar metallicities*)
\[Kappa]sb=4.7 10^20 Sqrt[Qg] Tp^(-15/4);
\[Kappa]sf=0.04 \[Kappa]sb;
(*If the free-free opacity dominates; Note we assume a H-He atmosphere with solar abundances*)
If[Free, \[Kappa]s= \[Kappa]sf, \[Kappa]s=\[Kappa]sb];

(*Finding \[Epsilon] as a function of \[Epsilon]s*)

K\[Nu]=(y/.(Solve[y^2 (1+y)==(1-1/\[Epsilon]s)^-2 f, y][[1]]));
\[Epsilon]=(1+K\[Nu]^-1)^-1;

\[Epsilon]s=\[Kappa]s/(\[Kappa]s+\[Kappa]es);

\[CapitalXi]=(0.873 \[Epsilon]s^(-1/6))/(1-0.127 \[Epsilon]s^(5/6)) 1/(1+(\[Epsilon]s^-1-1)^(2/3));
r=FindRoot[Teff^4-\[CapitalXi]  Tp^4==0, {Tp, Teff}(*, EvaluationMonitor:> Print[ Tp, "  ", Tp//Im, "  ",-\[CapitalXi]  Tp^4//Im] *)];

(*Plot[Teff^4- \[CapitalXi] Tp^4, {Tp, 0.99 Tp/.r, 1.01 Tp/.r}]*)
(*Correction to total flux from scattering*)

Tp=Tp/.r;

(*\[Kappa]s=(1+1/\[Epsilon]s)^-1;*)

(*x=\[Kappa]s/\[Kappa]es;*)
(*Solve cubic equation for the ratio of absorption opacity to scattering opacity*)

(*K=-(1/3)+2^(1/3)/(3 (-2+27 f x^2+Sqrt[-4+(-2+27 f x^2)^2])^(1/3))+(-2+27 f x^2+Sqrt[-4+(-2+27 f x^2)^2])^(1/3)/(3 2^(1/3));*)

(*flux1=15/\[Pi]^4  NIntegrate[(2 \[Epsilon]^(1/2))/(1+\[Epsilon]^(1/2)) E^-\[Xi]/f, {\[Xi],0, \[Infinity]}];
flux2= \[CapitalXi] ;
Print[flux1, flux2];*)


\[Xi]=h \[Nu]1/(kb Tp);
f=(\[Xi])^-3 (1-E^-\[Xi]);
(*K=-(1/3)+(2^(1/3) \[Kappa]es^2)/(3 (-2 \[Kappa]es^6+27 f \[Kappa]es^4 \[Kappa]s^2+3 Sqrt[3] Sqrt[-4 f \[Kappa]es^10 \[Kappa]s^2+27 f^2 \[Kappa]es^8 \[Kappa]s^4])^(1/3))+(-2 \[Kappa]es^6+27 f \[Kappa]es^4 \[Kappa]s^2+3 Sqrt[3] Sqrt[-4 f \[Kappa]es^10 \[Kappa]s^2+27 f^2 \[Kappa]es^8 \[Kappa]s^4])^(1/3)/(3 2^(1/3) \[Kappa]es^2);*)
(*Estimate density at the base of the photosphere for the peak frequency*)
\[Rho]p= (3 c Qg)/(4 (\[Gamma]) \[Sigma] Tp^4    \[Kappa]es^2 (1+K\[Nu])K\[Nu])/.{\[Nu]1->\[Nu]peak[Teff], \[Gamma]->4/3};
(*Estimate ratio of gas pressure to radiation pressure*)
pr=(4\[Sigma])/(3c) Tp^4;
pg= (\[Rho]p kb Tp)/(\[Mu]0 mp);

error=Abs[((\[CapitalXi] Tp^4/.r)-Teff^4)/Teff^4];

(*Print[NIntegrate[\[Pi] B[Teff], {\[Nu]1, 0, 10^20}], NIntegrate[2\[Pi] \[Epsilon]^(1/2)/(1+\[Epsilon]^(1/2)) B[Tp], {\[Nu]1, 0, 10^20}]//Chop," ",  \[Sigma] Teff^4];*)
If[error>prec, Print["Warning! Error in flux has exceeded the specified threshold. Error is ", error]];


{pr/pg, Tp, \[Epsilon], ((\[CapitalXi] Tp^4/.r)-Teff^4)/Teff^4}
]

(*Function takes list of annuli parameters and outputs a list of the total disk flux*)
GraybodyDiskFlux[params_]:=Module[{\[Epsilon],Tp,\[Nu]lo, \[Nu]hi, \[Nu]pts, tq, gb}, tq=dat[[All, {3,5}]];
tq=params[[All, {3,5}]];
\[Epsilon]=Map[Graybody[#[[1]],#[[2]]]&,10^tq];
Tp=\[Epsilon][[All,2]];
\[Epsilon]=\[Epsilon][[All,3]];

\[Nu]lo=Log10[2.4 10^14];
\[Nu]hi=Log10[5 10^19];
\[Nu]pts=Range[\[Nu]lo, \[Nu]hi, 0.02];

gb=2\[Pi] \[Epsilon]^(1/2)/(\[Epsilon]^(1/2)+1)*(B/@(Tp));
gb=gb/.\[Nu]1->10^\[Nu]pts;
gb=gb//Chop;

Transpose[{10^\[Nu]pts, 2*\[Pi]*params[[All,1]]*params[[All,2]]*gb//Total}]
]

(*Plots spectrum from unit 14 output file; this files is only produced when comptonization is turned on in tlusty*)
PlotSpec[ dir_, color_:Black,xrange_:{}, nmu_:10, mu_:1]:=Module[{dir2, dat, int, \[Lambda], \[Nu], \[Lambda]h, \[Nu]h, \[Nu]peak, imu, peak, prange},
dir2=StringReplace[dir, "/fort"~~__->""];
dat=ReadList[dir2<>"/fort.14", Number];
(*wavelengths in angstroms + fluxes*)
\[Lambda]h=Extract[dat, Position[dat, x_/;Length[x]==2]];
\[Lambda]h[[All, 2]]=\[Lambda]h[[All,2]] 4\[Pi];
\[Nu]h=\[Lambda]h;
(*Getting frequency information*)
\[Nu]h[[All,1]]=c/\[Lambda]h[[All, 1]] 10^8;
\[Nu]=\[Nu]h[[All,1]];
peak=\[Nu] \[Nu]h[[All,2]]//Max;
\[Nu]peak=\[Nu][[(\[Nu]h[[All,2]]//Ordering//#[[-1]]&)]];

(*Extract the specific intensities*)
imu=DeleteCases[dat, x_/;Length[x]==2];

imu=imu//Flatten//Partition[#, 2*nmu]&//#[[All, mu]]&;

If[xrange=={}, prange=Automatic, prange=xrange];
ListLogLogPlot[Transpose[{\[Nu], \[Nu]h[[All,2]] \[Nu]}], AxesLabel->{"\[Nu] [\!\(\*SuperscriptBox[\"s\", 
RowBox[{\"-\", \"1\"}]]\)]", "\[Nu] \!\(\*SubscriptBox[\"F\", \"\[Nu]\"]\) [erg \!\(\*SuperscriptBox[\"cm\", 
RowBox[{\"-\", \"2\"}]]\) \!\(\*SuperscriptBox[\"s\", 
RowBox[{\"-\", \"1\"}]]\)]"},  Joined->True, PlotStyle->Directive[color], PlotRange->{prange, {10^-6 peak , peak}}, PlotRangeClipping->False, ImageSize->Large]
]

(*Parses the opacity output from tlusty*)
Opac[ dir_, nd_:70]:=Module[{dir2,\[Kappa], \[Nu], dat},
dir2=StringReplace[dir, "/fort"~~__->""];
dat=ReadList[dir2<>"/fort.85", Number, RecordLists->True];
(*frequencies*)
\[Nu]=Extract[dat, Position[dat, x_/;Length[x]==2]];
\[Nu]=\[Nu][[All,2]];
\[Kappa]=DeleteCases[dat,x_/;Length[x]==2]//Flatten//Partition[#, nd]&;
Table[{\[Nu][[i]], \[Kappa][[i]]}, {i, 1, Length[\[Nu]]}]

]

(*Parses a model atmosphere from a unit 7 input file*)
ParseAtm[dir_]:=Module[{dir2, atm, m,nd, blocks},
dir2=StringReplace[dir, "/fort"~~__->""];
atm=Import[dir2<>"/fort.7", "Table"]//Flatten;
(*atm=ReadList[dir<>"/fort.7", Number];*)
nd=atm[[1]];
blocks=atm[[2]];
atm=atm[[3;;]];
m=atm[[1;;nd]];
atm=atm[[nd+1;;]];
atm=Partition[atm, blocks][[All, 1;;4]];
atm=Flatten/@Transpose[{m, atm}]

]

(*Plots convergence files outputted by tlusty*)
PConv[dir_]:=Module[{dir2, ConvData, maxchange, maxchange2, iter, colorscale, color, colors, params, t, q, Teff, Qg},
dir2=StringReplace[dir, "/fort"~~__->""];
params=StringCases[dir2, RegularExpression["t[\-0-9]*m[\-0-9]*q[\-0-9]*"]];
ConvData=Import[dir2<>"/fort.9", "Table", "HeaderLines"->3]//SplitBy[#, First]&;
ConvData=Reverse/@ConvData;

maxchange=Abs[ConvData[[All, All, -3]]];
maxchange2= Max/@maxchange;

iter= ConvData//Length;
If[iter>1,colorscale=Range[0, 1, 1/(iter-1)], colorscale={1}];
color=ColorData["Rainbow"][#]&;
colors=color/@colorscale;
GraphicsGrid[{{ListLogPlot[maxchange, Joined->True, Frame->True,FrameLabel->params, PlotStyle->colors, ImageSize->Medium],
ListLogPlot[maxchange2, Joined->True, Frame->True]}}]
]

(*Plot tlusty output for flux from TLUSTY unit 13 file*)
Options[PlotF]={optcol->Green, optrange->{10^14, 10^17}, optsize->Large}
PlotF[dir_?StringQ, OptionsPattern[]]:=Module[
{dir2, fluxdat, peak, sed, peakf, peakloc,bad, prange, color, xrange, size, myticks, myticks2, params, t, q,Qg, Teff,x},
color=OptionValue[optcol];
xrange=OptionValue[optrange];
size=OptionValue[optsize];

(*import data file*)
dir2=StringReplace[dir, "/fort"~~__->""];
fluxdat= Check[Import[dir2<>"/fort.13", "Table"][[All, 1;;2]],err];
If[fluxdat==err, Return[err]];

params=ParseFile[dir2];
t=params[[1]];
q=params[[3]];
Teff=10^t;
Qg=10^q;

fluxdat[[All,2]]=4\[Pi] fluxdat[[All,2]];
sed=Transpose[{fluxdat[[All,1]],fluxdat[[All, 1]]*fluxdat[[All,2]]}];
(*Output is sometimes formatted in such a way that mathematica does not recognize it as a real number. Tends to be at edges of frequency space*)
bad=DeleteCases[sed, x_/;x[[2]]\[Element]Reals];
sed=Complement[sed, bad];
(*Get location of peak for setting the range of the plot*)
peakloc=sed[[All,2]]//Ordering//#[[-1]]&;
peak=sed[[peakloc,2]];

(*Manually entered tick marks for plotting the spectrum*)
myticks={{1.`*^14, Style["\!\(\*SuperscriptBox[\"10\", \"14\"]\)", FontFamily->Times, FontSize->14]},{2.`*^14, ""},{3.`*^14, ""},{4.`*^14, ""},{5.`*^14,""},{6.`*^14, ""},{7.`*^14, "" },{8.`*^14, ""},{9.`*^14, "" },{1.`*^15, Style["\!\(\*SuperscriptBox[\"10\", \"15\"]\)", FontFamily->Times, FontSize->14]},{2.`*^15, ""},{3.`*^15, ""},{4.`*^15,""},{5.`*^15, ""},{6.`*^15, ""},{7.`*^15, "" },{8.`*^15, ""},{9.`*^15, ""},{1.`*^16,Style["\!\(\*SuperscriptBox[\"10\", \"16\"]\)", FontFamily->Times, FontSize->14]},{2.`*^16, ""},{3.`*^16, ""},{4.`*^16, ""},{5.`*^16, ""},{6.`*^16, ""},{7.`*^16, ""},{8.`*^16, ""},{9.`*^16, ""}, {1.`*^17,Style["\!\(\*SuperscriptBox[\"10\", \"17\"]\)", FontFamily->Times, FontSize->14] }};
myticks2={{1.`*^14, ""},{2.`*^14, ""},{3.`*^14, ""},{4.`*^14, ""},{5.`*^14,""},{6.`*^14, ""},{7.`*^14, "" },{8.`*^14, ""},{9.`*^14, "" },{1.`*^15, ""},{2.`*^15, ""},{3.`*^15, ""},{4.`*^15,""},{5.`*^15, ""},{6.`*^15, ""},{7.`*^15, "" },{8.`*^15, ""},{9.`*^15, ""},{1.`*^16,""},{2.`*^16, ""},{3.`*^16, ""},{4.`*^16, ""},{5.`*^16, ""},{6.`*^16, ""},{7.`*^16, ""},{8.`*^16, ""},{9.`*^16, ""}, {1.`*^17,"" }};

(*Actually plotting the spectrum*)
ListLogLogPlot[sed, Joined->True, PlotRange->{xrange, {10^-4  peak, 4 peak}},FrameLabel->{{ "\[Nu] \!\(\*SubscriptBox[\"F\", \"\[Nu]\"]\) [erg \!\(\*SuperscriptBox[\"cm\", 
RowBox[{\"-\", \"2\"}]]\) \!\(\*SuperscriptBox[\"s\", 
RowBox[{\"-\", \"1\"}]]\)]", None},{ "\[Nu] [\!\(\*SuperscriptBox[\"s\", 
RowBox[{\"-\", \"1\"}]]\)]", None}}, Frame->True,FrameTicks->{{Automatic, Automatic},{myticks, myticks2}},  PlotStyle->Directive[color], ImageSize->size]

]

(*User can also pass list of directories to PlotF ot easily plot multiple spectra. This function then call the PlotF function above for each of the directories*)
PlotF[dir_?ListQ, OptionsPattern[]]:=Module[{color, xrange, lc, lr, params,dir2, spectra, spectra2, reorder,size, speclegend, line}, 
color=OptionValue[optcol];

If[Not[ListQ[color]] , color={color}];
If[color=={"Grad"}, color=ColGrad[dir//Length]];
lc=color//Length;

xrange=OptionValue[optrange];
(*reorder=OptionValue[optorder];*)
size=OptionValue[optsize];
(*If[reorder,params=ParseFile/@dir;
dir2=dir[[params[[All, 1]]//Ordering]]//Reverse, dir2=dir];*)
(*lr=xrange//Length;*)
color=Table[color[[Mod[i, lc, 1]]],{i, 1, dir//Length}];
spectra=Table[PlotF[dir[[i]], optcol->color[[i]], optrange->xrange, optsize->size], {i, 1, dir//Length}];
spectra=DeleteCases[spectra, err];

params=ParseFile/@dir;

line=Line[{{0,0}, {1,0}}];
speclegend=Table[{Graphics[{color[[i]], line} ], "\!\(\*SubscriptBox[\"T\", \"eff\"]\)="<>ToString[params[[i, 1]]//N]<>" "<> "m="<>ToString[params[[i, 2]]//N]<>" "<> "\!\(\*SubscriptBox[\"Q\", \"g\"]\)="<>ToString[params[[i, 3]]//N]   } , {i, 1, dir//Length} ];
speclegend=Graphics[Legend[speclegend, LegendShadow->False]];

spectra2=spectra[[params[[All, 1]]//Ordering//Reverse]];
GraphicsRow[{Show[spectra2], speclegend}]
]

(*Parses file name to infer parameters associated with the model, differs from the above function only in t*)
ParseFile[file_]:=Module[{s, m,t, q, sm, st, sq},
s=StringCases[file, RegularExpression["t[\-0-9]*m[\-0-9]*q[\-0-9]*"]];
st=StringCases[s, RegularExpression["t[\-0-9]*"]]//Flatten;
t=StringCases[st, "t"~~x__->x]//Flatten//#[[1]]&//ToExpression//(# (1/10))&;
sm=StringCases[s, RegularExpression["m[\-0-9]*"]]//Flatten;
m=StringCases[sm, "m"~~x__->x]//Flatten//#[[1]]&//ToExpression//(# (1/10))&;
sq=StringCases[s, RegularExpression["q[\-0-9]*"]]//Flatten;
q=StringCases[sq, "q"~~x__->x]//Flatten//#[[1]]&//ToExpression//(# (1/10))&;
{t, m, q}
]

(*Function which gives model for a given set of parameters, with the assumed format corresponding to the output form from ParseFile*)
ReverseParse[params_?ListQ]:=Module[{},
"t"<>(params[[1]]*10//ToString)<>"m"<>(params[[2]]*10//ToString)<>"q"<>(params[[3]]*10//ToString)<>"/converged"
]

(*Compare tlusty unit 13 file to graybody output*)
Options[GraybodyCompare]={ optsize->Large}
GraybodyCompare[dir_?StringQ, OptionsPattern[]]:=
Module[{dir2, Tp, \[Epsilon], bb, gb, tlustyspec, gp, t, q, params, Teff, Qg, size}, 
size=OptionValue[optsize];

dir2=StringReplace[dir, "/fort"~~__->""];
params=ParseFile[dir2];
t=params[[1]];
q=params[[3]];
Teff=10^t;
Qg=10^q;

Tp=Graybody[Teff, Qg][[2]];
\[Epsilon]=Graybody[Teff, Qg][[3]];
bb=\[Pi] B[Teff] ;
gb=2\[Pi] \[Epsilon]^(1/2)/(1+\[Epsilon]^(1/2)) B[Tp] //Re;
(*Print[NIntegrate[{bb, gb}, {\[Nu]1, 0, 10^19}]];*)
tlustyspec=PlotF[dir2 <> "/fort.13",optrange->{0.1 \[Nu]peak[Teff], 10\[Nu]peak[Teff]}, optsize->size, optcol->Red];
gp=LogLogPlot[{bb \[Nu]1, gb \[Nu]1}, {\[Nu]1, 0.1 \[Nu]peak[Teff], 10\[Nu]peak[Teff]},  PlotStyle->{Black, Blue}];

Show[{tlustyspec, gp}]

]

(*Secondary graybody comparison function where extinction probability is evaluated using the opacity is outputted from TLUSTY itself. *)
Options[GraybodyCompare2]={optsize->Large}
GraybodyCompare2[dir_?StringQ, OptionsPattern[]]:=Module[{dir2, o, atm,  dens , \[Chi], \[Chi]all, \[Kappa]all, \[Kappa]r, \[Tau], near, pos, \[Kappa], \[Sigma], ne,\[Epsilon], \[Mu]e, T,\[Nu],z, nef, Tf, \[Chi]f, densf, m, tau0, keep, \[Delta]z, spec, planck, tlustyspec, gb, gp, gb1, gp1, params, t, q, Teff, Qg, size},
size=OptionValue[optsize];

dir2 =StringReplace[dir, "/fort"~~__->""];
params=ParseFile[dir2];
t=params[[1]];
q=params[[3]];
Teff=10^t;
Qg=10^q;

o=Opac[dir2<> "/fort.85"];
atm=ParseAtm[dir<>"/fort.7"];
tau0=2/3;

T=atm[[All,2]];
ne=atm[[All,3]];
dens=atm[[All,4]];
z=atm[[All,5]];

\[Delta]z=z//Differences;
o=o//Reverse;

\[Chi]=o[[All,2]];
\[Nu]=o[[All,1]];
\[Tau]=-#[[2;;]] \[Delta]z dens[[2;;]]&/@\[Chi];
keep=Position[\[Tau], x_/;Min[x]<tau0, 1];
\[Nu]=Extract[\[Nu], keep];
\[Tau]= Extract[\[Tau], keep];
\[Chi]=Extract[\[Chi], keep];

T=(Transpose[{#, T[[2;;]]}]&/@\[Tau]);
Tf=Interpolation/@T;
T=Table[Tf[[i]][tau0], {i, 1, Length[\[Nu]]}];

dens=(Transpose[{#, dens[[2;;]]}]&/@\[Tau]);
densf=Interpolation/@dens;
dens=Table[densf[[i]][tau0], {i, 1, Length[\[Nu]]}];

ne=(Transpose[{#, ne[[2;;]]}]&/@\[Tau]);
nef=Interpolation/@ne;
ne=Table[nef[[i]][tau0], {i, 1, Length[\[Nu]]}];

\[Chi]=Table[Transpose[{\[Tau][[i]],\[Chi][[i, 2;;]]}], {i, 1, \[Nu]//Length}];
\[Chi]f=Interpolation/@\[Chi];
\[Chi]=Table[\[Chi]f[[i]][tau0], {i, 1, Length[\[Nu]]}];

\[Mu]e=ne mp/dens;
\[Sigma]=\[Mu]e 0.4;
\[Kappa]=\[Chi]-\[Sigma];
\[Epsilon]=\[Kappa]/(\[Kappa]+\[Sigma]);

gb1=Table[{\[Nu][[i]], 2\[Pi] \[Nu][[i]](B[T[[i]]]/.{\[Nu]1->\[Nu][[i]]})   (\[Epsilon][[i]]^(1/2))/(1+\[Epsilon][[i]]^(1/2))}, {i, 1, Length[\[Nu]]}];

tlustyspec=PlotF[dir<> "/fort.13", optrange->{0.1 \[Nu]peak[Teff], 10 \[Nu]peak[Teff]}, optsize->size, optcol->Red];
gp1=ListLogLogPlot[gb1, Joined->True];
Show[tlustyspec, gp1]
]

(*Approximate color corrected blackbody based in part on the prescription in: 2012MNRAS .420.1848D*)
Options[ColCorrected]={optsize->Large}
ColCorrected[dir_?StringQ, OptionsPattern[]]:=Module[{dir2, params, Teff, Qg, t, q, tlustyspec,fcol, bb, size},
size=OptionValue[optsize];

dir2 =StringReplace[dir, "/fort"~~__->""];
params=ParseFile[dir2];
t=params[[1]];
q=params[[3]];
Teff=10^t;
Qg=10^q;
fcol=If[Teff>=10^5,(72 *keV/(kb*Teff))^(1/9), (Teff/(3 10^4))^0.82];

tlustyspec=PlotF[dir2<> "/fort.13",optrange-> {0.05 \[Nu]peak[Teff], 10 \[Nu]peak[Teff]}, optsize->size];

bb=LogLogPlot[{ \[Pi] \[Nu]1 B[Teff, fcol]}, {\[Nu]1, 0.05 \[Nu]peak[Teff], 10 \[Nu]peak[Teff]}];
Show[tlustyspec, bb]

]

(*Options[GraybodyCompare3]={ optsize->Large}
GraybodyCompare3[dir_?StringQ, OptionsPattern[]]:=
Module[{dir2, Tp, \[Epsilon], bb, gb, bp, tlustyspec, gp,params, t, q, m, Teff, Qg, Dmtot, size}, 
size=OptionValue[optsize];


dir2=StringReplace[dir, "/fort"~~__->""];
params=ParseFile[dir2];
t=params[[1]];
m=params[[2]];
q=params[[3]];

Teff=10^t;
Qg=10^q;
Dmtot=10^m;

gb=Graybody3[Teff, Qg, Dmtot];


bb=\[Pi] B[Teff] \[Nu]1;

bp=Plot[bb, {\[Nu]1, 0.1 \[Nu]peak[Teff], 10 \[Nu]peak[Teff]}];
gp=ListLogLogPlot[gb, Joined->True];
tlustyspec=PlotF[dir2<> "/fort.13"(*,Green*),optrange-> {0.1\[Nu]peak[Teff], 10 \[Nu]peak[Teff]}, optsize->size];

Show[tlustyspec, bp, gp]

]*)



