(* ::Package:: *)

Clear["Global`*"];
SetDirectory[NotebookDirectory[]];
Import["ToMatlab.m", "Package"];
Nmode=2;
T=1;
r=-0.01/8.69;

IN = IdentityMatrix[Nmode];

J0 = {{1,I},{1,-I}}/2;
Z0 = {{1,0},{0,-1}};
\[CapitalOmega]0 = {{0,1},{-1,0}};
J = KroneckerProduct[IN,J0];
Z = KroneckerProduct[IN,Z0];
\[CapitalOmega] = KroneckerProduct[IN,\[CapitalOmega]0];
X = KroneckerProduct[IN, {{0,1},{1,0}}];

\[CapitalSigma]sqz = {{Exp[-2r],0},{0,Exp[2r]}};
\[CapitalSigma]0 = KroneckerProduct[IN,\[CapitalSigma]sqz];

X = KroneckerProduct[IN, {{0,1},{1,0}}];

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
\[CapitalSigma] = Sqnd . \[CapitalSigma]0 . Transpose[Sqnd];
\[CapitalSigma] //MatrixForm*)

(*Through the channel*)
\[CapitalSigma] = ArrayFlatten[({
{\[CapitalSigma], 0},
{0, IdentityMatrix[2]}})];

B = ArrayFlatten[({
{IdentityMatrix[2Nmode-2], 0},
{0,{{Sqrt[T],0,Sqrt[1-T],0},{0,Sqrt[T],0,Sqrt[1-T]},{-Sqrt[1-T],0,Sqrt[T],0},{0,-Sqrt[1-T],0,Sqrt[T]}}}})];

\[CapitalSigma] = B . \[CapitalSigma] . Transpose[B];
\[CapitalSigma] = \[CapitalSigma][[1;;2Nmode, 1;;2Nmode]];

FullSimplify[MatrixForm[\[CapitalSigma]]]

FullSimplify[MatrixForm[X . (Z . J . \[CapitalSigma] . ConjugateTranspose[J] . Z-DiagonalMatrix[{1/2,1/2,1/2,1/2}])]]

(*CF after the channel*)

Chi[\[CapitalXi]_] = Exp[-(1/2)Transpose[\[CapitalXi]] . \[CapitalOmega] . \[CapitalSigma] . Transpose[\[CapitalOmega]] . \[CapitalXi]];
CF[a_,ac_,b_,bc_]=Simplify[Chi[{(a+ac)/2,(a-ac)/2/I,(b+bc)/2,(b-bc)/2/I}]]
(*CF[a_,ac_,b_,bc_]=Simplify[Chi[{(a+ac),(a-ac)/I,(b+bc),(b-bc)/I}]]*)

(*Photon subtraction*)
ChiPS[a_,ac_,b_,bc_]=Exp[-((a*ac+b*bc)/2)]*D[Exp[(a*ac+b*bc)/2]*CF[a,ac,b,bc],
{a,1},{ac,1}];

(*(*delocalized operations*)
ChiPS[a_,ac_,b_,bc_]=Simplify[
Exp[-((a*ac+b*bc)/2)]*D[Exp[(a*ac+b*bc)/2]*CF[a,ac,b,bc],{a,2},{ac,2}]+
Exp[-((a*ac+b*bc)/2)]*D[Exp[(a*ac+b*bc)/2]*CF[a,ac,b,bc],{b,2},{bc,2}]+
Exp[-((a*ac+b*bc)/2)]*D[Exp[(a*ac+b*bc)/2]*CF[a,ac,b,bc],{ac,2},{b,2}]+
Exp[-((a*ac+b*bc)/2)]*D[Exp[(a*ac+b*bc)/2]*CF[a,ac,b,bc],{bc,2},{a,2}]]*)

ChiPS[0,0,0,0]

ChiPSxp[x1_,p1_,x2_,p2_] = ChiPS[(x1+I*p1),(x1-I*p1),(x2+I*p2),(x2-I*p2)];
(*ChiPSxp[x1_,p1_,x2_,p2_] = ChiPS[(x1+I*p1)/2,(x1-I*p1)/2,(x2+I*p2)/2,(x2-I*p2)/2];*)

ChiPSxp[x1_,p1_,x2_,p2_] = Simplify[ChiPSxp[x1,p1,x2,p2]/ChiPSxp[0,0,0,0]];

Fmod[x_,p_]=Simplify[Chi[{-p,-x,-x,-p}]]
FPSmod[x_,p_]=Simplify[ChiPSxp[p,x,x,p]]

(*RangePlot=5;
Plot3D[Re[FPSmod[x,p]]/Re[Fmod[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},PlotPoints->200,PlotRange->All,RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2]]
Plot3D[Re[FPSmod[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},PlotPoints->200,PlotRange->All,RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2]]

FPSmod[1,3]
*)

(*TMSV test case*)
\[CapitalSigma]tmsv = {{Cosh[2 r],0,Sinh[2 r],0},{0,Cosh[2 r],0,-Sinh[2 r]},{Sinh[2 r],0,Cosh[2 r],0},{0,-Sinh[2 r],0,Cosh[2 r]}};
MatrixForm[\[CapitalSigma]tmsv]
MatrixForm[X . (Z . J . \[CapitalSigma]tmsv . ConjugateTranspose[J] . Z-DiagonalMatrix[{1/2,1/2,1/2,1/2}])]


FullSimplify[(U[[1,1]]*a1+U[[1,2]]*a2)(U[[2,1]]*a1+U[[2,2]]*a2)]


(*\:7ebf\:6027\:7c07\:6001*)
Clear["Global`*"];
Import["ToMatlab.m", "Package"];
Nmode = 4;
T = 1; (*Channel loss*)
SymbolicVector = Table[r[i], {i, 1, Nmode}];

IN = IdentityMatrix[Nmode];

J0 = {{1,I},{1,-I}}/2;
Z0 = {{1,0},{0,-1}};
\[CapitalOmega]0 = {{0,1},{-1,0}};
J = KroneckerProduct[IN,J0];
Z = KroneckerProduct[IN,Z0];
\[CapitalOmega] = KroneckerProduct[IN,\[CapitalOmega]0];

rdB = -1;
r = rdB/8.69;
\[CapitalSigma]sqz = {{Exp[-2r],0},{0,Exp[2r]}};
\[CapitalSigma]0 = KroneckerProduct[IN,\[CapitalSigma]sqz];

(*(*\:8003\:8651\:4e0d\:5747\:8861squeeze\:7684\:60c5\:51b5*)
r1 = -12/8.69;
r2 = -5/8.69;
r4 = -12/8.69;
r3 = -5/8.69;

r2 = -12/8.69;
r1 = -5/8.69;
r3 = -12/8.69;
r4 = -5/8.69;

r2 = -12/8.69;
r1 = -5/8.69;
r4 = -12/8.69;
r3 = -5/8.69;

r1 = .;
r2 = .;
r3 = .;
r4 = .;

\[CapitalSigma]0 = DiagonalMatrix[{Exp[-2r1], Exp[2r1], Exp[-2r2], Exp[2r2], Exp[-2r3], Exp[2r3] ,Exp[-2r4], Exp[2r4]}]*)

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

FullSimplify[\[CapitalSigma]] //MatrixForm

(*CF after the channel*)

Chi[\[CapitalXi]_] = Exp[-(1/2)Transpose[\[CapitalXi]] . \[CapitalOmega] . \[CapitalSigma] . Transpose[\[CapitalOmega]] . \[CapitalXi]];

(*\:6b64\:5904\:628ax\:548cp\:5206\:522b\:5f53\:4f5c\:5b9e\:90e8\:548c\:865a\:90e8*)
CF[a_,ac_,b_,bc_,c_,cc_,d_,dc_]=Simplify[Chi[{{a+ac}/2,{(a-ac)/2/I},{b+bc}/2,{(b-bc)/2/I},{c+cc}/2,{(c-cc)/2/I},{d+dc}/2,{(d-dc)/2/I}}]];

(*\:6d4b\:8bd5\:4e00\:4e0bchi\:7684\:4e24\:79cd\:8868\:8fbe\:5f0f\:7ed3\:679c\:662f\:5426\:4e00\:81f4*)
(*CF2[a_,ac_,b_,bc_,c_,cc_,d_,dc_]=Simplify[Exp[-0.5*{{ac,a,bc,b,cc,c,dc,d}}.Z.J.\[CapitalSigma].ConjugateTranspose[J].Z.{{a},{ac},{b},{bc},{c},{cc},{d},{dc}}]]
Simplify[CF[a,ac,b,bc,c,cc,d,dc]/CF2[a,ac,b,bc,c,cc,d,dc]]*)

(*ChiPS[a_,ac_,b_,bc_,c_,cc_,d_,dc_]=Exp[-((a*ac+b*bc+c*cc+d*dc)/2)]*D[Exp[(a*ac+b*bc+c*cc+d*dc)/2]*CF[a,ac,b,bc,c,cc,d,dc],a,ac,d,dc];*)
ChiPS[a_,ac_,b_,bc_,c_,cc_,d_,dc_]=Exp[-((a*ac+b*bc+c*cc+d*dc)/2)]*D[Exp[(a*ac+b*bc+c*cc+d*dc)/2]*CF[a,ac,b,bc,c,cc,d,dc],a,ac,b,bc,c,cc,d,dc];

ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_] = ChiPS[(x1+I*p1),(x1-I*p1),(x2+I*p2),(x2-I*p2),(x3+I*p3),(x3-I*p3),(x4+I*p4),(x4-I*p4)];
ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_] = ChiPSxp[x1,p1,x2,p2,x3,p3,x4,p4]/ChiPSxp[0,0,0,0,0,0,0,0];

Fmod[x_,p_]=Simplify[Chi[{-{p},-{x},{-x},{0},{p},{0},{x},{p}}]]
FPSmod[x_,p_]=Simplify[ChiPSxp[-p,-x,-x,0,p,0,x,p]]
(*Pi/Integrate[FPSmod[x,p],{x,-Infinity,Infinity},{p,-Infinity,Infinity}]*)
(*Collect[FPSmod[x,p]/Fmod[x,p],x^2,p^2]*)
(*ToMatlab[Fmod[x,p]];
ToMatlab[FPSmod[x,p]];*)
RangePlot = 2;
Plot3D[Re[FPSmod[x,p]]/Re[Fmod[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},PlotPoints->200,PlotRange->All,RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2]]
Fmod[x_,p_]=Simplify[Chi[{-p,-x,-x,0,p,0,x,p}]]


(*\:7ebf\:6027\:7c07\:6001*)
Clear["Global`*"];
Import["ToMatlab.m", "Package"];
Nmode = 5;
T = 1; (*Channel loss*)
nSubtraction = 1;
nSubtraction2 = 1;

SymbolicVector = Table[r[i], {i, 1, Nmode}];

IN = IdentityMatrix[Nmode];

J0 = {{1,I},{1,-I}}/2;
Z0 = {{1,0},{0,-1}};
\[CapitalOmega]0 = {{0,1},{-1,0}};
J = KroneckerProduct[IN,J0];
Z = KroneckerProduct[IN,Z0];
\[CapitalOmega] = KroneckerProduct[IN,\[CapitalOmega]0];
X = KroneckerProduct[IN, {{0,1},{1,0}}];

MatrixForm[X];

rdB = -6;
r = rdB/8.69;
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

FullSimplify[\[CapitalSigma]] //MatrixForm;

MatrixForm[(-0.5)X . (Z . J . \[CapitalSigma] . ConjugateTranspose[J] . Z-DiagonalMatrix[{0,0,1/2,1/2,0,0,0,0,1/2,1/2}])]

(*CF after the channel*)

Chi[\[CapitalXi]_] = Exp[-(1/2)Transpose[\[CapitalXi]] . \[CapitalOmega] . \[CapitalSigma] . Transpose[\[CapitalOmega]] . \[CapitalXi]];

(*\:6b64\:5904\:628ax\:548cp\:5206\:522b\:5f53\:4f5c\:5b9e\:90e8\:548c\:865a\:90e8*)
CF[a_,ac_,b_,bc_,c_,cc_,d_,dc_,e_,ec_]=Chi[{{a+ac}/2,{(a-ac)/2/I},{b+bc}/2,{(b-bc)/2/I},{c+cc}/2,{(c-cc)/2/I},{d+dc}/2,{(d-dc)/2/I},{e+ec}/2,{(e-ec)/2/I}}](*\:6b64\:5904\:628ax\:548cp\:5206\:522b\:5f53\:4f5c\:5b9e\:90e8\:548c\:865a\:90e8*)

(*\:6d4b\:8bd5\:4e00\:4e0bchi\:7684\:4e24\:79cd\:8868\:8fbe\:5f0f\:7ed3\:679c\:662f\:5426\:4e00\:81f4*)
(*CF2[a_,ac_,b_,bc_,c_,cc_,d_,dc_]=Simplify[Exp[-0.5*{{ac,a,bc,b,cc,c,dc,d}}.Z.J.\[CapitalSigma].ConjugateTranspose[J].Z.{{a},{ac},{b},{bc},{c},{cc},{d},{dc}}]]
Simplify[CF[a,ac,b,bc,c,cc,d,dc]/CF2[a,ac,b,bc,c,cc,d,dc]]*)

(*ChiPS[a_,ac_,b_,bc_,c_,cc_,d_,dc_,e_,ec_]=Exp[-((a*ac+b*bc+c*cc+d*dc+e*ec)/2)]*D[Exp[(a*ac+b*bc+c*cc+d*dc+e*ec)/2]*CF[a,ac,b,bc,c,cc,d,dc,e,ec],
{a,nSubtraction},{ac,nSubtraction},{e,nSubtraction2},{ec,nSubtraction2},{b,nSubtraction},{bc,nSubtraction},{d,nSubtraction2},{dc,nSubtraction2}];*)
ChiPS[a_,ac_,b_,bc_,c_,cc_,d_,dc_,e_,ec_]=Exp[-((a*ac+b*bc+c*cc+d*dc+e*ec)/2)]*D[Exp[(a*ac+b*bc+c*cc+d*dc+e*ec)/2]*CF[a,ac,b,bc,c,cc,d,dc,e,ec],
(*{a,nSubtraction},{ac,nSubtraction},{e,nSubtraction2},{ec,nSubtraction2},*)
{b,1},{bc,1},{d,1},{dc,1},
{a,1},{ac,1},{e,1},{ec,1}];

ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_,x5_,p5_] = ChiPS[(x1+I*p1),(x1-I*p1),(x2+I*p2),(x2-I*p2),(x3+I*p3),(x3-I*p3),(x4+I*p4),(x4-I*p4),(x5+I*p5),(x5-I*p5)];
Prob = ChiPSxp[0,0,0,0,0,0,0,0,0,0]
ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_,x5_,p5_] = ChiPSxp[x1,p1,x2,p2,x3,p3,x4,p4,x5,p5]/ChiPSxp[0,0,0,0,0,0,0,0,0,0];

Fmod[x_,p_]=Simplify[Chi[{-{p},-{x},{-x},{0},{p},{0},{x},{0},{-p},{x}}]];
FPSmod[x_,p_]=Simplify[ChiPSxp[-p,-x,-x,0,p,0,x,0,-p,x]];

RangePlot = 2.2;
Plotsim[x_,p_]=FullSimplify[Re[FPSmod[x,p]]/Re[Fmod[x,p]]];
P1=Plot3D[Plotsim[x,p],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
PlotPoints->300,PlotRange->All,
PlotStyle->Directive[RGBColor["#ff8066"],Opacity[0.99]],
MeshFunctions->{#3 &},Mesh->{{1}},MeshStyle->Dashed,
ColorFunction -> (Blend[{RGBColor["#9b89b3"], RGBColor["#845ec2"]}, #3] &),
RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2],
AxesLabel->{x,p,Subscript[\[Eta], PS]/\[Eta]},
ViewPoint-> {3,5,5}];

P2=Plot3D[1,{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
PlotPoints->200,PlotRange->All,
PlotStyle->Directive[RGBColor["#b0a8b9"],Opacity[0.5]],MeshFunctions->{#3 &},Mesh->{{1}},MeshStyle->Dashed,
AxesLabel->Automatic,
BoundaryStyle->None,ViewPoint-> {3,5,5}];

Pout = Show[P1,P2]

Export["C:\\Users\\PS\\Desktop\\lineareta.png",Pout,ImageResolution -> 600];

Fcoh[x_,p_] = Exp[-1/2(x^2+p^2)];

F = Integrate[Fcoh[-x,-p]*Fcoh[x,p]*Fmod[x,p],{x,-Infinity,Infinity},{p,-Infinity,Infinity}]/Pi
Fs = Integrate[Fcoh[-x,-p]*Fcoh[x,p]*FPSmod[x,p],{x,-Infinity,Infinity},{p,-Infinity,Infinity}]/Pi

data = Table[Re[FPSmod[x, y]]/Re[Fmod[x, y]], {x, -4, 4, 0.16}, {y, -4, 4, 0.16}];
Export["data.csv", ArrayReshape[Flatten[data],{51,51}], "CSV"];


(*\:7ebf\:6027\:7c07\:6001*)
Clear["Global`*"];
Import["ToMatlab.m", "Package"];
Nmode = 6;
T = 1; (*Channel loss*)
nSubtraction = 1;
nSubtraction2 = 1;

SymbolicVector = Table[r[i], {i, 1, Nmode}];

IN = IdentityMatrix[Nmode];

J0 = {{1,I},{1,-I}}/2;
Z0 = {{1,0},{0,-1}};
\[CapitalOmega]0 = {{0,1},{-1,0}};
J = KroneckerProduct[IN,J0];
Z = KroneckerProduct[IN,Z0];
\[CapitalOmega] = KroneckerProduct[IN,\[CapitalOmega]0];

rdB = -1;
r = rdB/8.69;
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

FullSimplify[\[CapitalSigma]];

(*CF after the channel*)

Chi[\[CapitalXi]_] = Exp[-(1/2)Transpose[\[CapitalXi]] . \[CapitalOmega] . \[CapitalSigma] . Transpose[\[CapitalOmega]] . \[CapitalXi]];

(*\:6b64\:5904\:628ax\:548cp\:5206\:522b\:5f53\:4f5c\:5b9e\:90e8\:548c\:865a\:90e8*)
CF[a1_,a1c_,a2_,a2c_,a3_,a3c_,a4_,a4c_,a5_,a5c_,a6_,a6c_]=Chi[{{a1+a1c}/2,{(a1-a1c)/2/I},{a2+a2c}/2,{(a2-a2c)/2/I},
{a3+a3c}/2,{(a3-a3c)/2/I},{a4+a4c}/2,{(a4-a4c)/2/I},{a5+a5c}/2,{(a5-a5c)/2/I},{a6+a6c}/2,{(a6-a6c)/2/I}}](*\:6b64\:5904\:628ax\:548cp\:5206\:522b\:5f53\:4f5c\:5b9e\:90e8\:548c\:865a\:90e8*)

ChiPS[a1_,a1c_,a2_,a2c_,a3_,a3c_,a4_,a4c_,a5_,a5c_,a6_,a6c_]=Exp[-((a1*a1c+a2*a2c+a3*a3c+a4*a4c+a5*a5c+a6*a6c)/2)]*
D[Exp[(a1*a1c+a2*a2c+a3*a3c+a4*a4c+a5*a5c+a6*a6c)/2]*CF[a1,a1c,a2,a2c,a3,a3c,a4,a4c,a5,a5c,a6,a6c],
{a1,nSubtraction},{a1c,nSubtraction},{a6,nSubtraction},{a6c,nSubtraction},
{a2,nSubtraction},{a2c,nSubtraction},{a5,nSubtraction},{a5c,nSubtraction},
{a3,nSubtraction},{a3c,nSubtraction},{a4,nSubtraction},{a4c,nSubtraction}];

ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_,x5_,p5_,x6_,p6_] = ChiPS[(x1+I*p1),(x1-I*p1),(x2+I*p2),(x2-I*p2),(x3+I*p3),(x3-I*p3),(x4+I*p4),(x4-I*p4),(x5+I*p5),(x5-I*p5),(x6+I*p6),(x6-I*p6)];
ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_,x5_,p5_,x6_,p6_] = ChiPSxp[x1,p1,x2,p2,x3,p3,x4,p4,x5,p5,x6,p6]/ChiPSxp[0,0,0,0,0,0,0,0,0,0,0,0];

Fmod[x_,p_]=Simplify[Chi[{-{p},-{x},{-x},{0},{p},{0},{x},{0},{-p},{0},{-x},{-p}}]]
FPSmod[x_,p_]=Simplify[ChiPSxp[-p,-x,-x,0,p,0,x,0,-p,0,-x,-p]]

RangePlot = 2;
Plot3D[Re[FPSmod[x,p]]/Re[Fmod[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},PlotPoints->200,PlotRange->All,RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2]]


(*\:7ebf\:6027\:7c07\:6001N=7*)
Clear["Global`*"];
Import["ToMatlab.m", "Package"];
Nmode = 7;
T = 1; (*Channel loss*)
rdB = -4.57041;

r = rdB/8.69;

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
nSubtraction = 1;
Chi[\[CapitalXi]_] = Exp[-(1/2)Transpose[\[CapitalXi]] . \[CapitalOmega] . \[CapitalSigma] . Transpose[\[CapitalOmega]] . \[CapitalXi]];
CF[a_,ac_,b_,bc_,c_,cc_,d_,dc_,e_,ec_,f_,fc_,g_,gc_]=Chi[{{a+ac}/2,{(a-ac)/2/I},{b+bc}/2,{(b-bc)/2/I},{c+cc}/2,{(c-cc)/2/I},{d+dc}/2,{(d-dc)/2/I},{e+ec}/2,{(e-ec)/2/I},{f+fc}/2,{(f-fc)/2/I},{g+gc}/2,{(g-gc)/2/I}}](*\:6b64\:5904\:628ax\:548cp\:5206\:522b\:5f53\:4f5c\:5b9e\:90e8\:548c\:865a\:90e8*)
(*ChiPS[a_,ac_,b_,bc_,c_,cc_,d_,dc_,e_,ec_,f_,fc_,g_,gc_]=Simplify[Exp[-((a*ac+b*bc+c*cc+d*dc+e*ec+f*fc+g*gc)/2)]*D[Exp[(a*ac+b*bc+c*cc+d*dc+e*ec+f*fc+g*gc)/2]*CF[a,ac,b,bc,c,cc,d,dc,e,ec,f,fc,g,gc],{b,nSubtraction},{bc,nSubtraction},{f,nSubtraction},{fc,nSubtraction}]];*)
ChiPS[a_,ac_,b_,bc_,c_,cc_,d_,dc_,e_,ec_,f_,fc_,g_,gc_]=Exp[-((a*ac+b*bc+c*cc+d*dc+e*ec+f*fc+g*gc)/2)]*D[Exp[(a*ac+b*bc+c*cc+d*dc+e*ec+f*fc+g*gc)/2]*CF[a,ac,b,bc,c,cc,d,dc,e,ec,f,fc,g,gc],
{a,nSubtraction},{ac,nSubtraction},{g,nSubtraction},{gc,nSubtraction},
{b,nSubtraction},{bc,nSubtraction},{f,nSubtraction},{fc,nSubtraction}];
ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_,x5_,p5_,x6_,p6_,x7_,p7_] = ChiPS[(x1+I*p1),(x1-I*p1),(x2+I*p2),(x2-I*p2),(x3+I*p3),(x3-I*p3),(x4+I*p4),(x4-I*p4),(x5+I*p5),(x5-I*p5),(x6+I*p6),(x6-I*p6),(x7+I*p7),(x7-I*p7)];
ChiPSxp[x1_,p1_,x2_,p2_,x3_,p3_,x4_,p4_,x5_,p5_,x6_,p6_,x7_,p7_] = ChiPSxp[x1,p1,x2,p2,x3,p3,x4,p4,x5,p5,x6,p6,x7,p7]/ChiPSxp[0,0,0,0,0,0,0,0,0,0,0,0,0,0];

Fmod[x_,p_]=Simplify[Chi[{{-p},{-x},{-x},{0},{p},{0},{x},{0},{-p},{0},{-x},{0},{p},{-x}}]]
FPSmod[x_,p_]=Simplify[ChiPSxp[-p,-x,-x,0,p,0,x,0,-p,0,-x,0,p,-x]];

RangePlot = 2.2;
(*Plot3D[{Re[FPSmod[x,p]]/Re[Fmod[x,p]],1},{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},PlotPoints->200,PlotRange->All,RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2]]*)
(*P1=Plot3D[Re[FPSmod[x,p]]/Re[Fmod[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
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

Pout = Show[P1,P2]*)
Fcoh[x_,p_] = Exp[-1/2(x^2+p^2)];
F = Integrate[Fcoh[-x,-p]*Fcoh[x,p]*Fmod[x,p],{x,-Infinity,Infinity},{p,-Infinity,Infinity}]/Pi
Fs = Integrate[Fcoh[-x,-p]*Fcoh[x,p]*FPSmod[x,p],{x,-Infinity,Infinity},{p,-Infinity,Infinity}]/Pi


Fcoh[-x,-p]*Fcoh[x,p]
