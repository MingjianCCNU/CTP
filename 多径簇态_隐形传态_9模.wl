(* ::Package:: *)

(*\:591a\:5f84\:7c07\:6001N=9*)
Clear["Global`*"];
Import["ToMatlab.m", "Package"];
Nmode = 9;
T = 1; (*Channel loss*)
rdB = -5;
r=rdB/8.69;
(*r = .;*)

Fcoh[x_,p_] = Exp[-1/2(x^2+p^2)];
IN = IdentityMatrix[Nmode];

\[CapitalOmega]0 = {{0,1},{-1,0}};
\[CapitalOmega] = KroneckerProduct[IN,\[CapitalOmega]0];

\[CapitalSigma]sqz = {{Exp[-2r],0},{0,Exp[2r]}};
\[CapitalSigma]0 = KroneckerProduct[IN,\[CapitalSigma]sqz];

(*Multi-rail cluster*)
K = Nmode; (* \:793a\:4f8b\:7ef4\:5ea6\:ff0c\:53ef\:4fee\:6539\:4e3a\:4efb\:610f\:6b63\:6574\:6570 *)
Uadj = SparseArray[
  {{1, 1} -> 0, {1, K} -> 0, {K, 1} -> 0, {K, K} -> 0} ~Join~
  Table[{i, 1} -> 1, {i, 2, K - 1}] ~Join~
  Table[{i, K} -> 1, {i, 2, K - 1}] ~Join~
  Table[{1, j} -> 1, {j, 2, K - 1}] ~Join~
  Table[{K, j} -> 1, {j, 2, K - 1}],
  {K, K}, 0
];

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

(*CF after the channel*)
nSubtraction = 1;
Chi[\[CapitalXi]_] = Exp[-(1/2)Transpose[\[CapitalXi]] . \[CapitalOmega] . \[CapitalSigma] . Transpose[\[CapitalOmega]] . \[CapitalXi]];
ChiInv[\[CapitalXi]_] = Chi[Sre . \[CapitalXi]];

ChiSqz[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_,x5_,p5_,x6_,p6_,x7_,p7_,x8_,p8_,x9_,p9_] = ChiInv[{{x1},{p1},{x2},{p2},{x3},{p3},{x4},{p4},{x5},{p5},{x6},{p6},{x7},{p7},{x8},{p8},{x9},{p9}}];

CF[a_,ac_,b_,bc_,c_,cc_,d_,dc_,e_,ec_,f_,fc_,g_,gc_,h_,hc_,i_,ic_]=Chi[{
{a+ac}/2,{(a-ac)/2/I},{b+bc}/2,{(b-bc)/2/I},{c+cc}/2,{(c-cc)/2/I},
{d+dc}/2,{(d-dc)/2/I},{e+ec}/2,{(e-ec)/2/I},{f+fc}/2,{(f-fc)/2/I},
{g+gc}/2,{(g-gc)/2/I},{h+hc}/2,{(h-hc)/2/I},{i+ic}/2,{(i-ic)/2/I}}](*\:6b64\:5904\:628ax\:548cp\:5206\:522b\:5f53\:4f5c\:5b9e\:90e8\:548c\:865a\:90e8*)
ChiPS[a_,ac_,b_,bc_,c_,cc_,d_,dc_,e_,ec_,f_,fc_,g_,gc_,h_,hc_,i_,ic_]=Exp[-((a*ac+b*bc+c*cc+d*dc+e*ec+f*fc+g*gc+h*hc+i*ic)/2)]*
D[Exp[(a*ac+b*bc+c*cc+d*dc+e*ec+f*fc+g*gc+h*hc+i*ic)/2]*CF[a,ac,b,bc,c,cc,d,dc,e,ec,f,fc,g,gc,h,hc,i,ic],{a,nSubtraction},{ac,nSubtraction},{i,nSubtraction},{ic,nSubtraction}];

ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_,x5_,p5_,x6_,p6_,x7_,p7_,x8_,p8_,x9_,p9_] = ChiPS[
(x1+I*p1),(x1-I*p1),(x2+I*p2),(x2-I*p2),(x3+I*p3),(x3-I*p3),(x4+I*p4),(x4-I*p4),(x5+I*p5),(x5-I*p5),(x6+I*p6),(x6-I*p6),(x7+I*p7),(x7-I*p7),(x8+I*p8),(x8-I*p8),(x9+I*p9),(x9-I*p9)];
ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_,x5_,p5_,x6_,p6_,x7_,p7_,x8_,p8_,x9_,p9_] = ChiPSxp[
x1,p1,x2,p2,x3,p3,x4,p4,x5,p5,x6,p6,x7,p7,x8,p8,x9,p9]/ChiPSxp[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

Fmod[x_,p_]=FullSimplify[Chi[{{p},{x},{x/7},{0},{x/7},{0},{x/7},{0},{x/7},{0},{x/7},{0},{x/7},{0},{x/7},{0},{-p},{x}}]]
FPSmod[x_,p_]=ChiPSxp[p,x,x/7,0,x/7,0,x/7,0,x/7,0,x/7,0,x/7,0,x/7,0,-p,x];

ChiPSxpMat[{{x1_},{p1_},{x2_},{p2_},{x3_},{p3_},{x4_},{p4_},{x5_},{p5_},{x6_},{p6_},{x7_},{p7_},{x8_},{p8_},{x9_},{p9_}}] = ChiPSxp[x1,p1,x2,p2,x3,p3,x4,p4,x5,p5,x6,p6,x7,p7,x8,p8,x9,p9];
ChiPsxpInv[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_,x5_,p5_,x6_,p6_,x7_,p7_,x8_,p8_,x9_,p9_] = ChiPSxpMat[Sre . {{x1},{p1},{x2},{p2},{x3},{p3},{x4},{p4},{x5},{p5},{x6},{p6},{x7},{p7},{x8},{p8},{x9},{p9}}];

VarX[x_] = -D[ChiSqz[x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],{x,2}];
VarX[x_] = -D[ChiSqz[0,0,x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],{x,2}];
ToMatlab[FullSimplify[VarX[0]]]
(*\:7b2c\:4e00\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],{x,2}];
ToMatlab[FullSimplify[VarXPS[0]]]
(*\:7b2c\:4e8c\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[0,0,x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],{x,2}];
ToMatlab[FullSimplify[VarXPS[0]]]
(*\:7b2c\:4e09\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[0,0,0,0,x,0,0,0,0,0,0,0,0,0,0,0,0,0],{x,2}];
ToMatlab[FullSimplify[VarXPS[0]]]
(*\:7b2c\:56db\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[0,0,0,0,0,0,x,0,0,0,0,0,0,0,0,0,0,0],{x,2}];
ToMatlab[FullSimplify[VarXPS[0]]]
(*\:7b2c\:4e94\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,0,0,0],{x,2}];
ToMatlab[FullSimplify[VarXPS[0]]]

RangePlot = 2.2;
Plotsim[x_,p_]=FullSimplify[Re[FPSmod[x,p]]/Re[Fmod[x,p]]];
P1=Plot3D[Plotsim[x,p],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
PlotPoints->300,PlotRange->All,
PlotStyle->Directive[RGBColor["#ff8066"],Opacity[0.99]],
MeshFunctions->{#3 &},Mesh->{{1}},MeshStyle->Dashed,
ColorFunction -> (Blend[{RGBColor["#9b89b3"], RGBColor["#845ec2"]}, #3] &),
RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2],
AxesLabel->{x,p,Subscript[\[Eta], PS]/\[Eta]}];

P2=Plot3D[1,{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
PlotPoints->200,PlotRange->All,
PlotStyle->Directive[RGBColor["#b0a8b9"],Opacity[0.5]],MeshFunctions->{#3 &},Mesh->{{1}},MeshStyle->Dashed,
AxesLabel->Automatic,
BoundaryStyle->None];

Pout = Show[P1,P2]

F = Integrate[Fcoh[-x,-p]*Fcoh[x,p]*Fmod[x,p],{x,-10,10},{p,-10,10}]/Pi
Fs = Integrate[Fcoh[-x,-p]*Fcoh[x,p]*FPSmod[x,p],{x,-10,10},{p,-10,10}]/Pi



