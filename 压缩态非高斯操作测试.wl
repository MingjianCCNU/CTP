(* ::Package:: *)

Clear["Global`*"];
SetDirectory[NotebookDirectory[]];
Import["ToMatlab.m", "Package"];

r = -1;
Nmode = 1;

IN = IdentityMatrix[Nmode];

J0 = {{1,I},{1,-I}}/2;
Z0 = {{1,0},{0,-1}};
\[CapitalOmega]0 = {{0,1},{-1,0}};
J = KroneckerProduct[IN,J0];
Z = KroneckerProduct[IN,Z0];
\[CapitalOmega] = KroneckerProduct[IN,\[CapitalOmega]0];

\[CapitalSigma] = {{Exp[-2r],0},{0,Exp[2r]}}

(*{{ac,a}} . Z . J
ConjugateTranspose[J] . Z . {{a},{ac}}*)

nsub = 2;
CF[a_,ac_] = Simplify[Exp[-{{ac,a}} . Z . J . \[CapitalSigma] . ConjugateTranspose[J] . Z . {{a},{ac}}/2]];
CFxp[x_,p_] = Simplify[CF[x+I*p,x-I*p]]
CFPs[a_,ac_] = Simplify[Exp[-a*ac/2]*D[Exp[a*ac/2]*CF[a,ac],{a,nsub},{ac,nsub}]];
CFPsxp[x_,p_] = CFPs[x+I*p,x-I*p]/CFPs[0,0];

RangePlot = 3;
DarkBlue = "#0089ba";
LightGrey = "#b0a8b9";
Orangee = "#ff6f91";
Yelloww = "#ffc75f";
P1 = Plot3D[Re[CFxp[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
PlotTheme->"Minimal",
PlotPoints->300,PlotRange->All,
PlotStyle->Directive[Opacity[0.99]],
Lighting -> Automatic,
AxesEdge -> Automatic,
MeshFunctions->{#3 &},Mesh->{{1}},MeshStyle->Dashed,
(*ColorFunction -> "BlueGreenYellow",*)
RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2],
ColorFunction -> (Blend[{RGBColor[Yelloww], RGBColor[Orangee]}, #3] &),
AxesLabel->{x,p,Subscript[\[Chi], coh]},
Boxed->False,PlotRangePadding->None,
AxesOrigin->{0,0,0},Axes->{True, True, False},AxesLabel->{x,p},AxesStyle -> Arrowheads[{0.03,0.03}],
ViewPoint-> 10{3,5,5}]
Export["C:\\Users\\PS\\Desktop\\sqzin.png",P1,ImageResolution -> 600];

(*\:751f\:6210\:6700\:7b80\:56fe\:7684\:4ee3\:7801\:6bb5*)
(*RangePlot = 3;
DarkBlue = "#0089ba";
LightGrey = "#b0a8b9";
Orangee = "#c34a36";
P1 = Plot3D[Re[CFxp[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
PlotTheme->"Minimal",
PlotPoints->300,PlotRange->All,
PlotStyle->Directive[Opacity[0.99]],
Lighting -> Automatic,
AxesEdge -> Automatic,
MeshFunctions->{#3 &},Mesh->{{1}},MeshStyle->Dashed,
RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2],
ColorFunction -> (Blend[{RGBColor["#ffc75f"], RGBColor["#ff6f91"]}, #3] &),
ViewPoint-> 10{3,5,5},Boxed->False,Axes->False,PlotRangePadding->None
(*AxesOrigin->{0,0,0},Axes->{True, True, False},AxesLabel->{x,p},AxesStyle -> Arrowheads[{0.03,0.03}]*)]*)
(*Export["C:\\Users\\PS\\Desktop\\cohIni.png",P1,ImageResolution -> 600];*)


CFxp[0,p]

(*RangePlot = 2;
Plot3D[Re[CFxp[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},PlotPoints->200,PlotRange->All,RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2]]*)
(*\:9a8c\:8bc1\:65b9\:5dee\:662f\:5426\:7b26\:5408\:6587\:732e[PhysRevA.111.043704]\:4e2d\:7684\:516c\:5f0f*)
Varx[x_] = -D[CFPsxp[x,0],{x,2}];
FullSimplify[Varx[0]]
(*Round[Varx[0],0.01]*)
Varp[p_] = -D[CFPsxp[0,p],{p,2}];

FullSimplify[Varp[0]]
(*Round[Varp[0],0.01]*)
(*Plot[{Varx[0]/Exp[2r]},{r,-2,2}]
Plot[{Varp[0]/Exp[-2r]},{r,-2,2}]*)
Vx = Round[Exp[2r]*(1+4Sinh[r]^2*(2Sinh[r]^2+Cosh[r]Sinh[r])/(2Sinh[r]^4+(Cosh[r]Sinh[r])^2)),0.01]
Vy = Round[Exp[-2r]*(1+4Sinh[r]^2*(2Sinh[r]^2-Cosh[r]Sinh[r])/(2Sinh[r]^4+(Cosh[r]Sinh[r])^2)),0.01]





(*TMSV N=2*)
Clear["Global`*"];
(*Import["ToMatlab.m", "Package"];*)
Nmode=2;

rdB = -15;
r = rdB/8.69;
r =  -.04;
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
Uadj = {{0,1},{1,0}};

(*Multi-rail optical cluster *)
U = (IN+I*Uadj) . Inverse[MatrixPower[Uadj . Uadj+IN,1/2]];

(*Reshape matrix*)
Ure = SparseArray[{{i_, j_} /; (Mod[i,2] == 1 && i==2j-1) || (Mod[i,2] == 0 && i-1==2(j-Nmode)-1) -> 1}, {2*Nmode, 2*Nmode},0] // Normal;
Ure //MatrixForm;

S = Join[Join[Re[U],-Im[U],2],Join[Im[U],Re[U],2]];

(*Reshaping matrix*)
Sre = Ure . S . Transpose[Ure];

Sre = KroneckerProduct

(*Print[
"Reordered Symplectic Transformation:",
Sre //MatrixForm]*)

Round[Sre,0.01] //MatrixForm;
\[CapitalSigma] = Sre . \[CapitalSigma]0 . Transpose[Sre];

Print[
"Covariance Matrix:",
\[CapitalSigma] //MatrixForm]

Chi[\[CapitalXi]_] = Exp[-(1/2)Transpose[\[CapitalXi]] . \[CapitalOmega] . \[CapitalSigma] . Transpose[\[CapitalOmega]] . \[CapitalXi]];
ChiInv[\[CapitalXi]_] = Chi[Sre . \[CapitalXi]];

ChiSqz[x1_,p1_,x2_,p2_] = ChiInv[{{x1},{p1},{x2},{p2}}];

CF[a_,ac_,b_,bc_]=Chi[{{a+ac}/2,{(a-ac)/2/I},{b+bc}/2,{(b-bc)/2/I}}];
ChiPS[a_,ac_,b_,bc_]=Exp[-((a*ac+b*bc)/2)]*D[Exp[(a*ac+b*bc)/2]*CF[a,ac,b,bc],{a,nSubtraction},{ac,nSubtraction},{b,nSubtraction2},{bc,nSubtraction2}];
ChiPSxp[x1_,p1_,x2_,p2_] = ChiPS[(x1+I*p1),(x1-I*p1),(x2+I*p2),(x2-I*p2)]/ChiPS[0,0,0,0];

ChiPSxpMat[{{x1_},{p1_},{x2_},{p2_}}] = ChiPSxp[x1,p1,x2,p2];
ChiPsxpInv[x1_,p1_,x2_,p2_] = ChiPSxpMat[Sre . {{x1},{p1},{x2},{p2}}];

(*Varac[a_,ac_,b_,bc_] = FullSimplify[-D[ChiPS[a,ac,b,bc]/ChiPS[0,0,0,0],{a,2}]];
Varac[0,0,0,0]
Varac[a_,ac_,b_,bc_] = FullSimplify[-D[ChiPS[a,ac,b,bc]/ChiPS[0,0,0,0],{ac,2}]];
Varac[0,0,0,0]*)

VarX[x_] = -D[ChiSqz[x,0,0,0],{x,2}];
VarX[x_] = -D[ChiSqz[0,0,x,0],{x,2}];
FullSimplify[VarX[0]]

(*\:7b2c\:4e00\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[x,0,0,0],{x,2}];
FullSimplify[VarXPS[0]]
(*\:7b2c\:4e8c\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[0,0,x,0],{x,2}];
FullSimplify[VarXPS[0]]
(*\:9a8c\:8bc1\:53d1\:73b0\:5b58\:5728\:5bf9\:79f0\:6027\:ff0c\:4e24\:7aef\:6a21\:7684\:65b9\:5dee\:76f8\:540c\:ff0c\:4e2d\:95f4\:6a21\:7684\:65b9\:5dee\:76f8\:540c*)

VarP[p_] = -D[ChiSqz[0,p,0,0],{p,2}];
VarP[p_] = -D[ChiSqz[0,0,0,p],{p,2}];
FullSimplify[VarP[0]]
VarPPS[p_] = -D[ChiPsxpInv[0,p,0,0],{p,2}];
FullSimplify[VarPPS[0]]
VarPPS[p_] = -D[ChiPsxpInv[0,0,0,p],{p,2}];
FullSimplify[VarPPS[0]]
Fmod[x_,p_]=Simplify[Chi[{{p},{x},{x},{p}}]]
FPSmod[x_,p_]=Simplify[ChiPSxp[p,x,x,p]];
RangePlot=0.5;
(*Plot3D[Re[FPSmod[x,p]]/Re[Fmod[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},PlotPoints->200,PlotRange->All,RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2]]*)
Plot3D[Re[FPSmod[x,p]]/Re[Fmod[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},PlotPoints->200,PlotRange->All]

(*\:9a8c\:8bc1\:65b9\:5dee\:662f\:5426\:7b26\:5408\:6587\:732e[PhysRevA.111.043704]\:4e2d\:7684\:516c\:5f0f*)
Vx = Exp[2r] *(1+2 Exp[r]Sinh[r]/Cosh[2r])
Vy = Exp[-2r](1-2 Exp[-r]Sinh[r]/Cosh[2r])

Plot[{Exp[2r]*(1+4Sinh[r]^2*(2Sinh[r]^2+Cosh[r]Sinh[r])/(2Sinh[r]^4+(Cosh[r]Sinh[r])^2)), Exp[2r] *(1+2 Exp[r]Sinh[r]/Cosh[2r])}, {r,-1,0}]


(*N=2, 2 independent sqz states with !!orthogonal!! sqz direction *)
Clear["Global`*"];
(*Import["ToMatlab.m", "Package"];*)
Nmode=2;

r = 0.4;
nSubtraction = 2;
nSubtraction2 = 2;

IN = IdentityMatrix[Nmode];

Sbs = KroneckerProduct[{{1,1},{1,-1}}/Sqrt[2],IN];

\[CapitalOmega]0 = {{0,1},{-1,0}};
\[CapitalOmega] = KroneckerProduct[IN,\[CapitalOmega]0];

\[CapitalSigma]sqz = {{Exp[-2r],0},{0,Exp[2r]}};
\[CapitalSigma]sqz2 = {{Exp[2r],0},{0,Exp[-2r]}};

(*Print["Initial CM:",
\[CapitalSigma]sqz //MatrixForm]*)

\[CapitalSigma]0 = KroneckerProduct[IN,\[CapitalSigma]sqz];

\[CapitalSigma] = \[CapitalSigma]0;
\[CapitalSigma] = ArrayFlatten[{{\[CapitalSigma]sqz,0},{0,\[CapitalSigma]sqz2}}];

Print[
"Covariance Matrix:",
\[CapitalSigma] //MatrixForm]

Chi[\[CapitalXi]_] = Exp[-(1/2)Transpose[\[CapitalXi]] . \[CapitalOmega] . \[CapitalSigma] . Transpose[\[CapitalOmega]] . \[CapitalXi]];

ChiSqz[x1_,p1_,x2_,p2_] = Chi[{{x1},{p1},{x2},{p2}}];

CF[a_,ac_,b_,bc_]=Chi[{{a+ac}/2,{(a-ac)/2/I},{b+bc}/2,{(b-bc)/2/I}}];
ChiPS[a_,ac_,b_,bc_]=
Exp[-((a*ac+b*bc)/2)]*D[Exp[(a*ac+b*bc)/2]*CF[a,ac,b,bc],{a,nSubtraction},{ac,nSubtraction},{b,0},{bc,0}]+
Exp[-((a*ac+b*bc)/2)]*D[Exp[(a*ac+b*bc)/2]*CF[a,ac,b,bc],{a,0},{ac,0},{b,nSubtraction},{bc,nSubtraction}]+
-Exp[-((a*ac+b*bc)/2)]*D[Exp[(a*ac+b*bc)/2]*CF[a,ac,b,bc],{a,nSubtraction},{ac,0},{b,0},{bc,nSubtraction}]+
-Exp[-((a*ac+b*bc)/2)]*D[Exp[(a*ac+b*bc)/2]*CF[a,ac,b,bc],{a,0},{ac,nSubtraction},{b,nSubtraction},{bc,0}];

(*ChiPS[a_,ac_,b_,bc_]=
Exp[-((a*ac+b*bc)/2)]*D[Exp[(a*ac+b*bc)/2]*CF[a,ac,b,bc],{a,nSubtraction},{ac,nSubtraction},{b,0},{bc,0}];*)

    
ChiPSxp[x1_,p1_,x2_,p2_] = FullSimplify[ChiPS[(x1+I*p1),(x1-I*p1),(x2+I*p2),(x2-I*p2)]/ChiPS[0,0,0,0]]
(*ChiPSxpIn[{{x1_},{p1_},{x2_},{p2_}}] = ChiPSxp[x1,p1,x2,p2];
ChiPSxpBs[x1_,p1_,x2_,p2_] = Simplify[ChiPSxpIn[Sbs.{{x1},{p1},{x2},{p2}}]]*)


VarX[x_] = -D[ChiSqz[x,0,0,0],{x,2}];
VarX[x_] = -D[ChiSqz[0,0,x,0],{x,2}];
(*FullSimplify[VarX[0]]*)

(*\:7b2c\:4e00\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPSxp[x,0,0,0],{x,2}];
FullSimplify[VarXPS[0]]
(*\:7b2c\:4e8c\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPSxp[0,0,x,0],{x,2}];
FullSimplify[VarXPS[0]]

VarP[p_] = -D[ChiPSxp[0,p,0,0],{p,2}];
FullSimplify[VarP[0]]
VarP[p_] = -D[ChiPSxp[0,0,0,p],{p,2}];
FullSimplify[VarP[0]]






(*\:591a\:5f84\:7c07\:6001N=4*)
Clear["Global`*"];
(*Import["ToMatlab.m", "Package"];*)
Nmode=4;

rdB = -9;
r=rdB/8.69;
nSubtraction = 1;
nSubtraction2 = 0;

IN = IdentityMatrix[Nmode];

\[CapitalOmega]0 = {{0,1},{-1,0}};
\[CapitalOmega] = KroneckerProduct[IN,\[CapitalOmega]0];

\[CapitalSigma]sqz = {{Exp[-2r],0},{0,Exp[2r]}};

(*Print["Initial CM:",
\[CapitalSigma]sqz //MatrixForm]*)

\[CapitalSigma]0 = KroneckerProduct[IN,\[CapitalSigma]sqz];

(*Multi-rail cluster*)
Uadj = {{0,1,1,0},{1,0,0,1},{1,0,0,1},{0,1,1,0}};

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

ChiSqz[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_] = ChiInv[{{x1},{p1},{x2},{p2},{x3},{p3},{x4},{p4}}];

CF[a_,ac_,b_,bc_,c_,cc_,d_,dc_]=Chi[{{a+ac}/2,{(a-ac)/2/I},{b+bc}/2,{(b-bc)/2/I},{c+cc}/2,{(c-cc)/2/I},{d+dc}/2,{(d-dc)/2/I}}];
ChiPS[a_,ac_,b_,bc_,c_,cc_,d_,dc_]=Exp[-((a*ac+b*bc+c*cc+d*dc)/2)]*D[Exp[(a*ac+b*bc+c*cc+d*dc)/2]*CF[a,ac,b,bc,c,cc,d,dc],{a,nSubtraction},{ac,nSubtraction},{d,nSubtraction2},{dc,nSubtraction2}];
ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_] = ChiPS[(x1+I*p1),(x1-I*p1),(x2+I*p2),(x2-I*p2),(x3+I*p3),(x3-I*p3),(x4+I*p4),(x4-I*p4)]/ChiPS[0,0,0,0,0,0,0,0];

ChiPSxpMat[{{x1_},{p1_},{x2_},{p2_},{x3_},{p3_},{x4_},{p4_}}] = ChiPSxp[x1,p1,x2,p2,x3,p3,x4,p4];
ChiPsxpInv[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_] = ChiPSxpMat[Sre . {{x1},{p1},{x2},{p2},{x3},{p3},{x4},{p4}}];


VarX[x_] = -D[ChiSqz[x,0,0,0,0,0,0,0],{x,2}];
VarX[x_] = -D[ChiSqz[0,0,x,0,0,0,0,0],{x,2}];
(*\:7b2c\:4e00\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[x,0,0,0,0,0,0,0],{x,2}];
(*\:7b2c\:4e8c\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[0,0,x,0,0,0,0,0],{x,2}];
(*\:7b2c\:4e09\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[0,0,0,0,x,0,0,0],{x,2}];
(*\:7b2c\:56db\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[0,0,0,0,0,0,x,0],{x,2}];
(*\:9a8c\:8bc1\:53d1\:73b0\:5b58\:5728\:5bf9\:79f0\:6027\:ff0c\:4e24\:7aef\:6a21\:7684\:65b9\:5dee\:76f8\:540c\:ff0c\:4e2d\:95f4\:6a21\:7684\:65b9\:5dee\:76f8\:540c*)
FullSimplify[VarX[0]]
FullSimplify[VarXPS[0]]

VarP[p_] = -D[ChiSqz[0,p,0,0,0,0,0,0],{p,2}];
VarP[p_] = -D[ChiSqz[0,0,0,p,0,0,0,0],{p,2}];
VarPPS[p_] = -D[ChiPsxpInv[0,p,0,0,0,0,0,0],{p,2}];
VarPPS[p_] = -D[ChiPsxpInv[0,0,0,p,0,0,0,0],{p,2}];
VarPPS[p_] = -D[ChiPsxpInv[0,0,0,0,0,p,0,0],{p,2}];
VarPPS[p_] = -D[ChiPsxpInv[0,0,0,0,0,0,0,p],{p,2}];
FullSimplify[VarP[0]]
FullSimplify[VarPPS[0]]

Fmod[x_,p_]=Simplify[Chi[{{p},{x},{x/2},{0},{x/2},{0},{-p},{x}}]]
FPSmod[x_,p_]=Simplify[ChiPSxp[p,x,x/2,0,x/2,0,-p,x]];
RangePlot=0.5;
(*Plot3D[Re[FPSmod[x,p]]/Re[Fmod[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},PlotPoints->200,PlotRange->All,RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2]]*)
Plot3D[Re[FPSmod[x,p]]/Re[Fmod[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},PlotPoints->200,PlotRange->All]


(*\:591a\:5f84\:7c07\:6001N=5*)
Clear["Global`*"];
(*Import["ToMatlab.m", "Package"];*)
Nmode = 5;

rdB = -12;
r = rdB/8.69;
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
FullSimplify[VarX[0]]
(*\:7b2c\:4e00\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[x,0,0,0,0,0,0,0,0,0],{x,2}];
FullSimplify[VarXPS[0]]
(*\:7b2c\:4e8c\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[0,0,x,0,0,0,0,0,0,0],{x,2}];
FullSimplify[VarXPS[0]]
(*\:7b2c\:4e09\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[0,0,0,0,x,0,0,0,0,0],{x,2}];
FullSimplify[VarXPS[0]]
(*\:7b2c\:56db\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[0,0,0,0,0,0,x,0,0,0],{x,2}];
FullSimplify[VarXPS[0]]
(*\:7b2c\:4e94\:4e2a\:6a21\:7684x\:573a\:65b9\:5dee*)
VarXPS[x_] = -D[ChiPsxpInv[0,0,0,0,0,0,0,0,x,0],{x,2}];
FullSimplify[VarXPS[0]]
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

Fmod[x_,p_]=Simplify[Chi[{{p},{x},{x/3},{0},{x/3},{0},{x/3},{0},{-p},{x}}]]
FPSmod[x_,p_]=Simplify[ChiPSxp[p,x,x/3,0,x/3,0,x/3,0,-p,x]]
RangePlot = 2.2;
(*Plot3D[{Re[FPSmod[x,p]]/Re[Fmod[x,p]],1},{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},PlotPoints->200,PlotRange->All,RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2]]*)
P1=Plot3D[Re[FPSmod[x,p]]/Re[Fmod[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
PlotPoints->300,PlotRange->All,
PlotStyle->Directive[RGBColor["#ff8066"],Opacity[0.99]],
MeshFunctions->{#3 &},Mesh->{{1}},MeshStyle->Dashed,
ColorFunction -> (Blend[{RGBColor["#9b89b3"], RGBColor["#845ec2"]}, #3] &),
RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2],
AxesLabel->{x,p,Subscript[\[Eta], PS]/\[Eta]},ViewPoint-> {3,5,5}]

P2=Plot3D[1,{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
PlotPoints->200,PlotRange->All,
PlotStyle->Directive[RGBColor["#b0a8b9"],Opacity[0.5]],MeshFunctions->{#3 &},Mesh->{{1}},MeshStyle->Dashed,
AxesLabel->Automatic,
BoundaryStyle->None,ViewPoint-> {3,5,5}];
(*ContourPlot[Re[FPSmod[x,p]]/Re[Fmod[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},PlotLegends->Automatic,ContourLabels->True,ColorFunction->(If[#<1,Lighter[Blue,#],Lighter[Blue,#]]&),ColorFunctionScaling->False]*)

Integrate[FPSmod[x,p]*Exp[-(x^2+p^2)]/Pi,{x,-100,100},{p,-100,100}]
Integrate[Fmod[x,p]*Exp[-(x^2+p^2)]/Pi,{x,-100,100},{p,-100,100}]
Pout = Show[P1,P2]
(*Export["r1.svg",Show[P1,P2]]*)
(*Export["C:\\Users\\PS\\Desktop\\ratior3t.png",Pout,ImageResolution -> 600];*)


Clear["Global`*"];
SetDirectory[NotebookDirectory[]];
Import["ToMatlab.m", "Package"];

r = -0/8.69;
r = .;
Nmode = 1;

IN = IdentityMatrix[Nmode];

J0 = {{1,I},{1,-I}}/2;
Z0 = {{1,0},{0,-1}};
\[CapitalOmega]0 = {{0,1},{-1,0}};
J = KroneckerProduct[IN,J0];
Z = KroneckerProduct[IN,Z0];
\[CapitalOmega] = KroneckerProduct[IN,\[CapitalOmega]0];

\[CapitalSigma] = {{Exp[-2r],0},{0,Exp[2r]}};

(*{{ac,a}} . Z . J
ConjugateTranspose[J] . Z . {{a},{ac}}*)


CF[a_,ac_] = FullSimplify[Exp[-(1/2){{ac,a}} . Z . J . \[CapitalSigma] . ConjugateTranspose[J] . Z . {{a},{ac}}]];
CFxp[x_,p_] = FullSimplify[CF[(x+I*p),(x-I*p)]];
VarX[x_]=-D[CFxp[x,0],{x,2}];
VarX[0];
VarP[x_]=-D[CFxp[0,x],{x,2}];
VarP[0];

CFsqz[x1_,p1_,x2_,p2_] = FullSimplify[CFxp[x1,p1]*CFxp[x2,p2]]
CFtmsv[x1_,p1_,x2_,p2_] = FullSimplify[CFsqz[x1,p1-x2,x2,p2-x1]]
FullSimplify[CFtmsv[p,x,x,p]];

(*CFtmsv[x1_,p1_,x2_,p2_] = Simplify[
CFxp[(x1+x2)/Sqrt[2],(p1+p2)/Sqrt[2]]*
CFxp[(p1-p2)/Sqrt[2],(x1-x2)/Sqrt[2]]];
FullSimplify[CFtmsv[x,-p,x,p]]*)

Vx1
VarX[x_]=-D[CFtmsv[x,0,0,0],{x,2}];
FullSimplify[VarX[0]]

Vq1
VarX[x_]=-D[CFtmsv[0,x,0,0],{x,2}];
FullSimplify[VarX[0]]

Covq1x2
CovX[x1_,x2_]=-D[CFtmsv[0,x1,x2,0],x1,x2];
FullSimplify[CovX[0,0]]

Covx1q2
CovX[x1_,x2_]=-D[CFtmsv[x1,0,0,x2],x1,x2];
FullSimplify[CovX[0,0]]

Vx2
VarP[x_]=-D[CFtmsv[0,0,x,0],{x,2}];
FullSimplify[VarP[0]]

Vq2
VarP[x_]=-D[CFtmsv[0,0,0,x],{x,2}];
FullSimplify[VarP[0]]

