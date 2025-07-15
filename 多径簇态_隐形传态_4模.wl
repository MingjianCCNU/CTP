(* ::Package:: *)

(*\:591a\:5f84\:7c07\:6001N=4*)
Clear["Global`*"];
SetDirectory[NotebookDirectory[]];
Import["ToMatlab.m", "Package"];
Nmode = 4;
T=1; (*Channel loss*)
rdB = -8;
r=rdB/8.69;
nSubtraction = 1;
IN = IdentityMatrix[Nmode];

\[CapitalOmega]0 = {{0,1},{-1,0}};
\[CapitalOmega] = KroneckerProduct[IN,\[CapitalOmega]0];

\[CapitalSigma]sqz = {{Exp[-2r],0},{0,Exp[2r]}};
\[CapitalSigma]0 = KroneckerProduct[IN,\[CapitalSigma]sqz];

(*Multi-rail cluster*)
Uadj = {{0,1,1,0},{1,0,0,1},{1,0,0,1},{0,1,1,0}};

(*Multi-rail optical cluster *)
U = (IN+I*Uadj) . Inverse[MatrixPower[Uadj . Uadj+IN,1/2]];
Round[U,0.01] //MatrixForm;

(*Reshape matrix*)
Ure = SparseArray[{{i_, j_} /; (Mod[i,2] == 1 && i==2j-1) || (Mod[i,2] == 0 && i-1==2(j-Nmode)-1) -> 1}, {2*Nmode, 2*Nmode},0] // Normal;
Ure //MatrixForm;

S = Join[Join[Re[U],-Im[U],2],Join[Im[U],Re[U],2]];
(*Reshaping matrix*)
Sre = Ure . S . Transpose[Ure];
Round[Sre,0.01] //MatrixForm;
\[CapitalSigma] = Sre . \[CapitalSigma]0 . Transpose[Sre];
\[CapitalSigma] //MatrixForm

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

(*CF after the channel*)

Chi[\[CapitalXi]_] = Exp[-(1/2)Transpose[\[CapitalXi]] . \[CapitalOmega] . \[CapitalSigma] . Transpose[\[CapitalOmega]] . \[CapitalXi]];
CF[a_,ac_,b_,bc_,c_,cc_,d_,dc_]=Simplify[Chi[{{a+ac}/2,{(a-ac)/2/I},{b+bc}/2,{(b-bc)/2/I},{c+cc}/2,{(c-cc)/2/I},{d+dc}/2,{(d-dc)/2/I}}]](*\:6b64\:5904\:628ax\:548cp\:5206\:522b\:5f53\:4f5c\:5b9e\:90e8\:548c\:865a\:90e8*)
(*CF[a_,ac_,b_,bc_,c_,cc_,d_,dc_]=Simplify[Chi[{{a+ac},{(a-ac)/I},{b+bc},{(b-bc)/I},{c+cc},{(c-cc)/I},{d+dc},{(d-dc)/I}}]](*\:6b64\:5904\:628ax\:548cp\:5206\:522b\:5f53\:4f5c\:5b9e\:90e8\:548c\:865a\:90e8*)*)
(*ChiPS[a_,ac_,b_,bc_,c_,cc_,d_,dc_]=Exp[-((a*ac+b*bc+c*cc+d*dc)/2)]*D[Exp[(a*ac+b*bc+c*cc+d*dc)/2]*CF[a,ac,b,bc,c,cc,d,dc],{a,10},{ac,10},b,bc,c,cc,d,dc];*)
(*ChiPS[a_,ac_,b_,bc_,c_,cc_,d_,dc_]=Exp[((a*ac+b*bc+c*cc+d*dc)/2)]*D[Exp[-(a*ac+b*bc+c*cc+d*dc)/2]*CF[a,ac,b,bc,c,cc,d,dc],b,bc,c,cc];*)
ChiPS[a_,ac_,b_,bc_,c_,cc_,d_,dc_]=Exp[-((a*ac+b*bc+c*cc+d*dc)/2)]*D[Exp[(a*ac+b*bc+c*cc+d*dc)/2]*CF[a,ac,b,bc,c,cc,d,dc],{a,nSubtraction},{ac,nSubtraction},{d,nSubtraction},{dc,nSubtraction}];

ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_] = ChiPS[(x1+I*p1),(x1-I*p1),(x2+I*p2),(x2-I*p2),(x3+I*p3),(x3-I*p3),(x4+I*p4),(x4-I*p4)];
(*ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_] = ChiPS[(x1+I*p1)/2,(x1-I*p1)/2,(x2+I*p2)/2,(x2-I*p2)/2,(x3+I*p3)/2,(x3-I*p3)/2,(x4+I*p4)/2,(x4-I*p4)/2];*)
Prob = ChiPSxp[0,0,0,0,0,0,0,0]
ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_] = ChiPSxp[x1,p1,x2,p2,x3,p3,x4,p4]/ChiPSxp[0,0,0,0,0,0,0,0];

Fmod[x_,p_]=Simplify[Chi[{{p},{x},{x/2},{0},{x/2},{0},{-p},{x}}]]
FPSmod[x_,p_]=Simplify[ChiPSxp[p,x,x/2,0,x/2,0,-p,x]]
(*Pi/Integrate[FPSmod[x,p],{x,-Infinity,Infinity},{p,-Infinity,Infinity}]*)
(*Collect[FPSmod[x,p]/Fmod[x,p],x^2,p^2]*)
(*ToMatlab[Fmod[x,p]];
ToMatlab[FPSmod[x,p]];*)
RangePlot=2;
(*Plot3D[Re[FPSmod[x,p]]/Re[Fmod[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},PlotPoints->200,PlotRange->All,RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2]]*)
surfacePlot = Plot3D[Re[FPSmod[x,p]]/Re[Fmod[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},PlotPoints->200,PlotRange->All]

data = Table[Re[FPSmod[x,p]]/Re[Fmod[x,p]],{x,-2,2,0.04},{p,-2,2,0.04}];
Export["data.csv", data, "CSV"];




CFcoh[x_,p_] = Exp[-1/2*(x^2+p^2)];
CFcoh[x,p]CFcoh[x,p]
Integrate[CFcoh[x,p]CFcoh[x,p],{x,-2.0,2.0},{p,-2.0,2.0}]/Pi

