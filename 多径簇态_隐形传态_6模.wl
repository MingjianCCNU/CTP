(* ::Package:: *)

(*\:591a\:5f84\:7c07\:6001N=6*)
Clear["Global`*"];
Import["ToMatlab.m", "Package"];
Nmode=6;
T=1; (*Channel loss*)
rdB = -3;
r=rdB/8.69;

IN = IdentityMatrix[Nmode];

\[CapitalOmega]0 = {{0,1},{-1,0}};
\[CapitalOmega] = KroneckerProduct[IN,\[CapitalOmega]0];

\[CapitalSigma]sqz = {{Exp[-2r],0},{0,Exp[2r]}};
\[CapitalSigma]0 = KroneckerProduct[IN,\[CapitalSigma]sqz];

(*Multi-rail cluster*)
Uadj = {{0,1,1,1,1,0},{1,0,0,0,0,1},{1,0,0,0,0,1},{1,0,0,0,0,1},{1,0,0,0,0,1},{0,1,1,1,1,0}};

(*Multi-rail optical cluster *)
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

(*Through the channel*)
\[CapitalSigma] = ArrayFlatten[({
{\[CapitalSigma], 0},
{0, IdentityMatrix[2]}})];

B = ArrayFlatten[({
{IdentityMatrix[2Nmode-2], 0},
{0,{{Sqrt[T],0,Sqrt[1-T],0},{0,Sqrt[T],0,Sqrt[1-T]},{-Sqrt[1-T],0,Sqrt[T],0},{0,-Sqrt[1-T],0,Sqrt[T]}}}})];

\[CapitalSigma] = B . \[CapitalSigma] . Transpose[B];
\[CapitalSigma] = \[CapitalSigma][[1;;2Nmode, 1;;2Nmode]];

\[CapitalSigma] //MatrixForm

(*CF after the channel*)
nSubtraction = 2;
Chi[\[CapitalXi]_] = Exp[-(1/2)Transpose[\[CapitalXi]] . \[CapitalOmega] . \[CapitalSigma] . Transpose[\[CapitalOmega]] . \[CapitalXi]];
CF[a_,ac_,b_,bc_,c_,cc_,d_,dc_,e_,ec_,f_,fc_]=Chi[{{a+ac}/2,{(a-ac)/2/I},{b+bc}/2,{(b-bc)/2/I},{c+cc}/2,{(c-cc)/2/I},{d+dc}/2,{(d-dc)/2/I},{e+ec}/2,{(e-ec)/2/I},{f+fc}/2,{(f-fc)/2/I}}](*\:6b64\:5904\:628ax\:548cp\:5206\:522b\:5f53\:4f5c\:5b9e\:90e8\:548c\:865a\:90e8*)
ChiPS[a_,ac_,b_,bc_,c_,cc_,d_,dc_,e_,ec_,f_,fc_]=Simplify[Exp[-((a*ac+b*bc+c*cc+d*dc+e*ec+f*fc)/2)]*D[Exp[(a*ac+b*bc+c*cc+d*dc+e*ec+f*fc)/2]*CF[a,ac,b,bc,c,cc,d,dc,e,ec,f,fc],{a,nSubtraction},{ac,nSubtraction},{f,nSubtraction},{fc,nSubtraction}]];

ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_,x5_,p5_,x6_,p6_] = ChiPS[(x1+I*p1),(x1-I*p1),(x2+I*p2),(x2-I*p2),(x3+I*p3),(x3-I*p3),(x4+I*p4),(x4-I*p4),(x5+I*p5),(x5-I*p5),(x6+I*p6),(x6-I*p6)];
ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_,x5_,p5_,x6_,p6_] = ChiPSxp[x1,p1,x2,p2,x3,p3,x4,p4,x5,p5,x6,p6]/ChiPSxp[0,0,0,0,0,0,0,0,0,0,0,0];

Fmod[x_,p_]=Simplify[Chi[{{p},{x},{x/4},{0},{x/4},{0},{x/4},{0},{x/4},{0},{-p},{x}}]]
FPSmod[x_,p_]=Simplify[ChiPSxp[p,x,x/4,0,x/4,0,x/4,0,x/4,0,-p,x]]
(*Pi/Integrate[FPSmod[x,p],{x,-Infinity,Infinity},{p,-Infinity,Infinity}]*)
(*Collect[FPSmod[x,p]/Fmod[x,p],x^2,p^2]*)
(*ToMatlab[Fmod[x,p]];
ToMatlab[FPSmod[x,p]];*)
RangePlot=4;
Plot3D[Re[FPSmod[x,p]]/Re[Fmod[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},PlotPoints->200,PlotRange->All,RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2]]
Plot[{Re[FPSmod[x,x]],Re[Fmod[x,x]]},{x,-RangePlot,RangePlot}]




