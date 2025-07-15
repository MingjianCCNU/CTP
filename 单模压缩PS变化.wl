(* ::Package:: *)

Clear["Global`*"];

r = 0;
Sn  = 4;

CFsqz[x_,p_]=Exp[-1/2(Exp[2r]x^2+Exp[-2r]p^2)]*LaguerreL[Sn, Exp[2r]x^2+Exp[-2r]p^2]
CFsqz[0,0]

CFsqzD2x[x_] = -D[CFsqz[x,0],{x,2}];
CFsqzD2p[p_] = -D[CFsqz[0,p],{p,2}];

CM = {{CFsqzD2p[0],0},{0,CFsqzD2x[0]}}//MatrixForm



Clear["Global`*"];
Import["ToMatlab.m", "Package"];
r = .;
Sn  = 8;

CM = {{Exp[-2r],0},{0,Exp[2r]}}//MatrixForm;

CFsqz[x_,p_] = Exp[-1/2(Exp[2r]x^2+Exp[-2r]p^2)];
CFsqzl[a_,ac_] = CFsqz[(a+ac)/2,(a-ac)/2/I];

CFsqzsub[a_,ac_] = -Exp[-1/2(a*ac)]D[Exp[1/2(a*ac)]CFsqzl[a,ac],{a,Sn},{ac,Sn}];
CFsqzsubxp[x_,p_] = FullSimplify[CFsqzsub[x+I*p,x-I*p]/CFsqzsub[0,0]];

CFsqzD2x[x_] = FullSimplify[-D[CFsqzsubxp[x,0],{x,2}]];
CFsqzD2p[p_] = FullSimplify[-D[CFsqzsubxp[0,p],{p,2}]];

CM = FullSimplify[{{CFsqzD2p[0],0},{0,CFsqzD2x[0]}}]//MatrixForm
ToMatlab[FullSimplify[CFsqzD2x[0]]]
(*Solve[CFsqzD2x[0]-\[ExponentialE]^(2 r)==0,r,Reals]
Solve[CFsqzD2p[0]-\[ExponentialE]^(-2 r)==0,r,Reals]*)

Plot[{CFsqzD2x[0],E^(2 r)},{r,-2,0}]


Log[3]/2.0
ArcTanh[0.5]



