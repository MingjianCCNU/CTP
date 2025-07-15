(* ::Package:: *)

(*\:591a\:5f84\:7c07\:6001N=5*)
Clear["Global`*"];
(*Import["ToMatlab.m", "Package"];*)
Nmode = 5;

rdB = -2.5;
r = rdB/8.69;
(*r = .;*)

nSubtraction = 1;
nSubtraction2 = 1;

IN = IdentityMatrix[Nmode];

\[CapitalOmega]0 = {{0,1},{-1,0}};
\[CapitalOmega] = KroneckerProduct[IN,\[CapitalOmega]0];

\[CapitalSigma]sqz = {{Exp[-2r],0},{0,Exp[2r]}};

(*Print["Initial CM:",
\[CapitalSigma]sqz //MatrixForm]*)

\[CapitalSigma]0 = KroneckerProduct[IN,\[CapitalSigma]sqz];

(*Multi-rail cluster*)
Uadj = {{0,1,1,1,0},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{0,1,1,1,0}};

(*Multi-rail optical cluster *)
U = (IN+I*Uadj) . Inverse[MatrixPower[Uadj . Uadj+IN,1/2]];

(*Reshape matrix*)
Ure = SparseArray[{{i_, j_} /; (Mod[i,2] == 1 && i==2j-1) || (Mod[i,2] == 0 && i-1==2(j-Nmode)-1) -> 1}, {2*Nmode, 2*Nmode},0] // Normal;
Ure //MatrixForm;

S = Join[Join[Re[U],-Im[U],2],Join[Im[U],Re[U],2]];

(*Print[
"Symplectic Transformation:",
S //MatrixForm]*)

(*Reshaping matrix*)
Sre = Ure . S . Transpose[Ure];
Round[Sre,0.01] //MatrixForm;
\[CapitalSigma] = Sre . \[CapitalSigma]0 . Transpose[Sre];

(*Print[
"Covariance Matrix:",
\[CapitalSigma] //MatrixForm]*)

Chi[\[CapitalXi]_] = Exp[-(1/2)Transpose[\[CapitalXi]] . \[CapitalOmega] . \[CapitalSigma] . Transpose[\[CapitalOmega]] . \[CapitalXi]];
ChiInv[\[CapitalXi]_] = Chi[Sre . \[CapitalXi]];

ChiSqz[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_,x5_,p5_] = ChiInv[{{x1},{p1},{x2},{p2},{x3},{p3},{x4},{p4},{x5},{p5}}];

CF[a_,ac_,b_,bc_,c_,cc_,d_,dc_,e_,ec_]=Chi[{{a+ac}/2,{(a-ac)/2/I},{b+bc}/2,{(b-bc)/2/I},{c+cc}/2,{(c-cc)/2/I},{d+dc}/2,{(d-dc)/2/I},{e+ec}/2,{(e-ec)/2/I}}];
ChiPS[a_,ac_,b_,bc_,c_,cc_,d_,dc_,e_,ec_]=Exp[-((a*ac+b*bc+c*cc+d*dc+e*ec)/2)]*D[Exp[(a*ac+b*bc+c*cc+d*dc+e*ec)/2]*CF[a,ac,b,bc,c,cc,d,dc,e,ec],
{a,nSubtraction},{ac,nSubtraction},{e,nSubtraction2},{ec,nSubtraction2}];
(*ChiPS[a_,ac_,b_,bc_,c_,cc_,d_,dc_,e_,ec_]=Exp[-((a*ac+b*bc+c*cc+d*dc+e*ec)/2)]*D[Exp[(a*ac+b*bc+c*cc+d*dc+e*ec)/2]*CF[a,ac,b,bc,c,cc,d,dc,e,ec],{b,nSubtraction},{bc,nSubtraction},{c,nSubtraction2},{cc,nSubtraction2},{d,nSubtraction2},{dc,nSubtraction2}];*)
ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_,x5_,p5_] = ChiPS[(x1+I*p1),(x1-I*p1),(x2+I*p2),(x2-I*p2),(x3+I*p3),(x3-I*p3),(x4+I*p4),(x4-I*p4),(x5+I*p5),(x5-I*p5)]/ChiPS[0,0,0,0,0,0,0,0,0,0];
ChiPS[0,0,0,0,0,0,0,0,0,0]
ChiPSxpMat[{{x1_},{p1_},{x2_},{p2_},{x3_},{p3_},{x4_},{p4_},{x5_},{p5_}}] = ChiPSxp[x1,p1,x2,p2,x3,p3,x4,p4,x5,p5];
ChiPsxpInv[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_,x5_,p5_] = ChiPSxpMat[Sre . {{x1},{p1},{x2},{p2},{x3},{p3},{x4},{p4},{x5},{p5}}];

VarX[x_] = -D[ChiSqz[x,0,0,0,0,0,0,0,0,0],{x,2}];
VarX[x_] = -D[ChiSqz[0,0,x,0,0,0,0,0,0,0],{x,2}];
ToMatlab[FullSimplify[VarX[0]]]
(*\:7b2c\:4e00\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[x,0,0,0,0,0,0,0,0,0],{x,2}];
ToMatlab[FullSimplify[VarXPS[0]]]
(*\:7b2c\:4e8c\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[0,0,x,0,0,0,0,0,0,0],{x,2}];
ToMatlab[FullSimplify[VarXPS[0]]]
(*\:7b2c\:4e09\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[0,0,0,0,x,0,0,0,0,0],{x,2}];
ToMatlab[FullSimplify[VarXPS[0]]]
(*\:7b2c\:56db\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[0,0,0,0,0,0,x,0,0,0],{x,2}];
ToMatlab[FullSimplify[VarXPS[0]]]
(*\:7b2c\:4e94\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[0,0,0,0,0,0,0,0,x,0],{x,2}];
ToMatlab[FullSimplify[VarXPS[0]]]
(*\:9a8c\:8bc1\:53d1\:73b0\:5b58\:5728\:5bf9\:79f0\:6027\:ff0c\:4e24\:7aef\:6a21\:7684\:65b9\:5dee\:76f8\:540c\:ff0c\:4e2d\:95f4\:6a21\:7684\:65b9\:5dee\:76f8\:540c*)
VarP[p_] = -D[ChiSqz[0,p,0,0,0,0,0,0,0,0],{p,2}];
VarP[p_] = -D[ChiSqz[0,0,0,p,0,0,0,0,0,0],{p,2}];
VarPPS[p_] = -D[ChiPsxpInv[0,p,0,0,0,0,0,0,0,0],{p,2}];
FullSimplify[VarPPS[0]]
VarPPS[p_] = -D[ChiPsxpInv[0,0,0,p,0,0,0,0,0,0],{p,2}];
FullSimplify[VarPPS[0]]
VarPPS[p_] = -D[ChiPsxpInv[0,0,0,0,0,p,0,0,0,0],{p,2}];
FullSimplify[VarPPS[0]]
VarPPS[p_] = -D[ChiPsxpInv[0,0,0,0,0,0,0,p,0,0],{p,2}];
FullSimplify[VarPPS[0]]
VarPPS[p_] = -D[ChiPsxpInv[0,0,0,0,0,0,0,0,0,p],{p,2}];
FullSimplify[VarPPS[0]]

Fmod[x_,p_]=Simplify[Chi[{-{p},-{x},{-x},{0},{p},{0},{x},{0},{-p},{x}}]];
FPSmod[x_,p_]=Simplify[ChiPSxp[-p,-x,-x,0,p,0,x,0,-p,x]];
RangePlot = 2.2;
DarkBlue = "#0089ba";
LightGrey = "#b0a8b9";
Orangee = "#ff6f91";
Yelloww = "#ffc75f";


PlotData[x_,p_] = Simplify[ChiPsxpInv[x,p,0,0,0,0,0,0,0,0]];

RangePlot = 3;
DarkBlue = "#0089ba";
LightGrey = "#b0a8b9";
Orangee = "#ff6f91";
Yelloww = "#ffc75f";
P1 = Plot3D[Re[PlotData[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
PlotTheme->"Minimal",
PlotPoints->300,PlotRange->All,
PlotStyle->Directive[Opacity[0.99]],
Lighting -> Automatic,
AxesEdge -> Automatic,
MeshFunctions->{#3 &},Mesh->{{1}},MeshStyle->Dashed,
RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2],
ColorFunction -> (Blend[{RGBColor[LightGrey], RGBColor[DarkBlue]}, #3] &),
AxesLabel->{x,p,Subscript[\[Chi], coh]},
Boxed->False,PlotRangePadding->None,
Axes->False,
ViewPoint-> 10{3,5,5}]
Export["C:\\Users\\PS\\Desktop\\endMode.png",P1,ImageResolution -> 300, Background->None];


P2=Plot3D[1,{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
PlotPoints->200,PlotRange->All,
PlotStyle->Directive[RGBColor["#b0a8b9"],Opacity[0.5]],MeshFunctions->{#3 &},Mesh->{{1}},MeshStyle->Dashed,
AxesLabel->Automatic,
BoundaryStyle->None];
(*ContourPlot[Re[FPSmod[x,p]]/Re[Fmod[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},PlotLegends->Automatic,ContourLabels->True,ColorFunction->(If[#<1,Lighter[Blue,#],Lighter[Blue,#]]&),ColorFunctionScaling->False]*)

Pout = Show[P1,P2]
(*Export["r1.svg",Show[P1,P2]]*)
(*Export["C:\\Users\\PS\\Desktop\\ratior3.png",Pout,ImageResolution -> 600];*)

