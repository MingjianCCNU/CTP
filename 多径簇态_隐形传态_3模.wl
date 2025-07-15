(* ::Package:: *)

(*\:591a\:5f84\:7c07\:6001\:6216\:7ebf\:6027\:7c07\:6001N=3*)
Clear["Global`*"];
(*Import["ToMatlab.m", "Package"];*)
Nmode=3;

T=1;
rdB = -5;
r=rdB/8.69;


IN = IdentityMatrix[Nmode];

\[CapitalOmega]0 = {{0,1},{-1,0}};
\[CapitalOmega] = KroneckerProduct[IN,\[CapitalOmega]0];

\[CapitalSigma]sqz = {{Exp[-2r],0},{0,Exp[2r]}};
\[CapitalSigma]0 = KroneckerProduct[IN,\[CapitalSigma]sqz];

(*Linear cluster*)
Uadj = SparseArray[{{i_, j_} /; 
(i==j-1)||(i==j+1) -> 1}
, {Nmode, Nmode},0] // Normal;

(*Linear optical cluster *)
U = (IN+I*Uadj) . Inverse[MatrixPower[Uadj . Uadj+IN,1/2]];
U //MatrixForm;

(*Reshape matrix*)
Ure = SparseArray[{{i_, j_} /; (Mod[i,2] == 1 && i==2j-1) || (Mod[i,2] == 0 && i-1==2(j-Nmode)-1) -> 1}, {2*Nmode, 2*Nmode},0] // Normal;
Ure //MatrixForm;

S = Join[Join[Re[U],-Im[U],2],Join[Im[U],Re[U],2]];
S //MatrixForm;

(*Reshaping matrix*)
Sre = Ure . S . Transpose[Ure];
Round[Sre,0.01] //MatrixForm;
\[CapitalSigma] = Sre . \[CapitalSigma]0 . Transpose[Sre];
\[CapitalSigma] //MatrixForm;

(*(*Canonical Cluster*)
Sqnd = IdentityMatrix[2*Nmode];
For[n=1,n<=Nmode,n++,
For[k=n,k<=Nmode,k++,
Sqnd=Sqnd . Normal[SparseArray[{\
{i_,i_} -> 1,\
{i_,j_} /;  Uadj[[n,k]]==1 && ((i==2n && j==2k-1) || (i==2k && j==2n-1)) -> 1
},{2*Nmode, 2*Nmode}]];
]];
Sqnd//MatrixForm;
\[CapitalSigma] = Sqnd . \[CapitalSigma]0 . Transpose[Sqnd];*)

(*(*Through the channel*)
\[CapitalSigma] = ArrayFlatten[({
{\[CapitalSigma], 0},
{0, IdentityMatrix[2]}})];

B = ArrayFlatten[({
{IdentityMatrix[2Nmode-2], 0},
{0,{{Sqrt[T],0,Sqrt[1-T],0},{0,Sqrt[T],0,Sqrt[1-T]},{-Sqrt[1-T],0,Sqrt[T],0},{0,-Sqrt[1-T],0,Sqrt[T]}}}})];

\[CapitalSigma] = B . \[CapitalSigma] . Transpose[B];
\[CapitalSigma] = \[CapitalSigma][[1;;2Nmode, 1;;2Nmode]];

\[CapitalSigma] //MatrixForm;*)

(*CF after the channel*)

Chi[\[CapitalXi]_] = Exp[-(1/2)Transpose[\[CapitalXi]] . \[CapitalOmega] . \[CapitalSigma] . Transpose[\[CapitalOmega]] . \[CapitalXi]];
CF[a_,ac_,b_,bc_,c_,cc_]=Chi[{{a+ac}/2,{(a-ac)/2/I},{b+bc}/2,{(b-bc)/2/I},{c+cc}/2,{(c-cc)/2/I}}];
(*CF[a_,ac_,b_,bc_,c_,cc_]=Chi[{{a+ac},{(a-ac)/I},{b+bc},{(b-bc)/I},{c+cc},{(c-cc)/I}}];*)

(*ChiPS[a_,ac_,b_,bc_,c_,cc_]=Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CF[a,ac,b,bc,c,cc],a,ac,{b,1},{bc,1},c,cc];*)
ChiPS[a_,ac_,b_,bc_,c_,cc_]=Simplify[Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CF[a,ac,b,bc,c,cc],
(*{a,1},{ac,2},{b,1},{bc,0},{c,0},{cc,0},*)
{a,1},{ac,1},{b,0},{bc,0},{c,1},{cc,1}]]
(*ChiPS[a_,ac_,b_,bc_,c_,cc_]=Exp[-((a*ac)/2)]*D[Exp[(a*ac)/2]*CF[a,ac,b,bc,c,cc],a,ac];*)
(*ChiPS[a_,ac_,b_,bc_,c_,cc_]=Exp[-((b*bc)/2)]*D[Exp[(b*bc)/2]*CF[a,ac,b,bc,c,cc],{b,3},{bc,3}];*)
(*ChiPS[a_,ac_,b_,bc_,c_,cc_]=Exp[-((a*ac)/2)]*D[Exp[(a*ac)/2]*CF[a,ac,b,bc,c,cc],{a,2},{ac,2}];*)

(*ChiPS[a_,ac_,b_,bc_,c_,cc_]=Exp[-(a*ac)/2]*D[Exp[(a*ac)]*D[Exp[-(a*ac)/2]*CF[a,ac,b,bc,c,cc],a,ac],a,ac];*)

ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_] = ChiPS[(x1+I*p1),(x1-I*p1),(x2+I*p2),(x2-I*p2),(x3+I*p3),(x3-I*p3)];
(*ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_] = ChiPS[(x1+I*p1)/2,(x1-I*p1)/2,(x2+I*p2)/2,(x2-I*p2)/2,(x3+I*p3)/2,(x3-I*p3)/2];*)

ChiPSxp[0,0,0,0,0,0]

ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_] = Simplify[ChiPSxp[x1,p1,x2,p2,x3,p3]/ChiPSxp[0,0,0,0,0,0]];

Fmod[x_,p_]=Simplify[Chi[{{-p},{-x},{-x},{0},{p},{-x}}]]
FPSmod[x_,p_]=Simplify[ChiPSxp[-p,-x,-x,0,p,-x]];

RangePlot=2;
Plot3D[Re[FPSmod[x,p]]/Re[Fmod[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},PlotPoints->200,PlotRange->All,RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2]]
Plot[Re[FPSmod[0,x]]/Re[Fmod[0,x]],{x,-RangePlot,RangePlot}]

Re[Fmod[1,3]]
Re[FPSmod[1,3]]




(*PS\:5728U\:524d\:7b49\:6548\:6027\:6d4b\:8bd5*)
Clear["Global`*"];
(*Import["ToMatlab.m", "Package"];*)
Nmode=3;

rdB = -5;
r=rdB/8.69;

IN = IdentityMatrix[Nmode];

\[CapitalOmega]0 = {{0,1},{-1,0}};
\[CapitalOmega] = KroneckerProduct[IN,\[CapitalOmega]0];

\[CapitalSigma]sqz = {{Exp[-2r],0},{0,Exp[2r]}};
\[CapitalSigma]0 = KroneckerProduct[IN,\[CapitalSigma]sqz];

(*Linear cluster*)
Uadj = SparseArray[{{i_, j_} /; 
(i==j-1)||(i==j+1) -> 1}
, {Nmode, Nmode},0] // Normal;

(*Linear optical cluster *)
U = (IN+I*Uadj) . Inverse[MatrixPower[Uadj . Uadj+IN,1/2]];
U //MatrixForm

(*Reshape matrix*)
Ure = SparseArray[{{i_, j_} /; (Mod[i,2] == 1 && i==2j-1) || (Mod[i,2] == 0 && i-1==2(j-Nmode)-1) -> 1}, {2*Nmode, 2*Nmode},0] // Normal;
Ure //MatrixForm;

S = Join[Join[Re[U],-Im[U],2],Join[Im[U],Re[U],2]];
S //MatrixForm;

(*Reshaping matrix*)
Sre = Ure . S . Transpose[Ure];
Round[Sre,0.01] //MatrixForm;
\[CapitalSigma] = Sre . \[CapitalSigma]0 . Transpose[Sre];
\[CapitalSigma] //MatrixForm;

Chi[\[CapitalXi]_] = Exp[-(1/2)Transpose[\[CapitalXi]] . \[CapitalOmega] . \[CapitalSigma] . Transpose[\[CapitalOmega]] . \[CapitalXi]];
CF[a_,ac_,b_,bc_,c_,cc_]=Chi[{{a+ac}/2,{(a-ac)/2/I},{b+bc}/2,{(b-bc)/2/I},{c+cc}/2,{(c-cc)/2/I}}];
ChiPS[a_,ac_,b_,bc_,c_,cc_]=Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CF[a,ac,b,bc,c,cc],a,ac];
ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_] = ChiPS[(x1+I*p1),(x1-I*p1),(x2+I*p2),(x2-I*p2),(x3+I*p3),(x3-I*p3)];
ChiPSxp[0,0,0,0,0,0]
ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_] = Simplify[ChiPSxp[x1,p1,x2,p2,x3,p3]/ChiPSxp[0,0,0,0,0,0]];

ChiPSxpMatin[{{x1_},{p1_},{x2_},{p2_},{x3_},{p3_}}]=ChiPSxp[x1,p1,x2,p2,x3,p3];

FullSimplify[ChiPSxpMatin[Sre . {{x1},{p1},{x2},{p2},{x3},{p3}}]]

(*\:5728\:7ecf\:8fc7\:7ebf\:6027\:7f51\:7edc\:4e4b\:524d\:5bf9\:6a21\:65bd\:52a0\:76f8\:5e72\:7684PS*)
ChiSqz[\[CapitalXi]_] = Exp[-(1/2)Transpose[\[CapitalXi]] . \[CapitalOmega] . \[CapitalSigma]0 . Transpose[\[CapitalOmega]] . \[CapitalXi]];
U1 = 1/2+1/(2 Sqrt[3]);
(*U1 = 1/3 (3/2+Sqrt[3]/2);*)
U1c = Conjugate[U1];
U2 = -(I/Sqrt[3]);
(*U2 = \[ImaginaryI]/Sqrt[3];*)
U2c = Conjugate[U2];
U3 = -(1/2)+1/(2 Sqrt[3]);
(*U3 = 1/3 (-(3/2)+Sqrt[3]/2);*)
U3c = Conjugate[U3];
CFSqz[a_,ac_,b_,bc_,c_,cc_]=ChiSqz[{{a+ac}/2,{(a-ac)/2/I},{b+bc}/2,{(b-bc)/2/I},{c+cc}/2,{(c-cc)/2/I}}];
CFSqzps[a_,ac_,b_,bc_,c_,cc_] = 
U1*U1c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],a,ac]+
U2*U2c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],b,bc]+
U3*U3c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],c,cc]+
U1*U2c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],b,ac]+
U2*U1c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],a,bc]+
U1*U3c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],c,ac]+
U3*U1c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],a,cc]+
U2*U3c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],c,bc]+
U3*U2c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],b,cc];

CFSqzps[a_,ac_,b_,bc_,c_,cc_] = 
U1*U1c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],a,ac]+
U2*U2c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],b,bc]+
U3*U3c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],c,cc]+
U1*U2c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],a,bc]+
U2*U1c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],b,ac]+
U1*U3c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],a,cc]+
U3*U1c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],c,ac]+
U2*U3c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],b,cc]+
U3*U2c*Exp[-((a*ac+b*bc+c*cc)/2)]*D[Exp[(a*ac+b*bc+c*cc)/2]*CFSqz[a,ac,b,bc,c,cc],c,bc];
CFSqzps[0,0,0,0,0,0]
CFSqzpsxp[x1_,p1_,x2_,p2_,x3_,p3_] = FullSimplify[CFSqzps[(x1+I*p1),(x1-I*p1),(x2+I*p2),(x2-I*p2),(x3+I*p3),(x3-I*p3)]/CFSqzps[0,0,0,0,0,0]]


Expand[-6I(U[[1,1]]a1+U[[1,2]]a2+U[[1,3]]a3)(U[[2,1]]a1+U[[2,2]]a2+U[[2,3]]a3)]

Expand[-6(U[[1,1]]a1+U[[1,2]]a2+U[[1,3]]a3)(U[[3,1]]a1+U[[3,2]]a2+U[[3,3]]a3)]


 


Round[Abs[1I/6+1I/2Sqrt[3]],0.01]
