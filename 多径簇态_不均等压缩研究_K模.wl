(* ::Package:: *)

(*\:591a\:5f84\:7c07\:6001K\:53d6\:53d8\:91cf*)
Clear["Global`*"];
(* \:83b7\:53d6\:5f53\:524d\:7b14\:8bb0\:672c\:7684\:8def\:5f84 *)
currentNotebookPath = DirectoryName[NotebookFileName[]];
(* \:8bbe\:7f6e\:5f53\:524d\:5de5\:4f5c\:76ee\:5f55\:4e3a\:8be5\:8def\:5f84 *)
SetDirectory[currentNotebookPath];

Import["ToMatlab.m", "Package"];
r = .;
Nmode = 9;
T = 1;
IN = IdentityMatrix[Nmode];

\[CapitalOmega]0 = {{0,1},{-1,0}};
\[CapitalOmega] = KroneckerProduct[IN,\[CapitalOmega]0];

(*\:8003\:8651\:4e0d\:5747\:8861squeeze\:7684\:60c5\:51b5*)

K = Nmode; (* \:8bbe\:7f6e\:5411\:91cf\:7ef4\:5ea6 *)

$Assumptions = {r1k \[Element] Reals, r2 \[Element] Reals};
r[i_] := r; (* \:5b9a\:4e49r\:51fd\:6570\:ff0c\:8fd9\:91cc\:91c7\:7528\:7ebf\:6027\:589e\:957f\:793a\:4f8b *)
(*r[1] = -1.25;
r[K] = -1.25;*)
\[CapitalSigma]0 = DiagonalMatrix[Flatten[Table[{Exp[-2 r[i]], Exp[2 r[i]]}, {i, K}]]];

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
MatrixForm[Uadj]; (* \:663e\:793a\:77e9\:9635 *)

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
\[CapitalSigma] //MatrixForm;

(*CF after the channel*)
Chi[\[CapitalXi]_] = Exp[-(1/2)Transpose[\[CapitalXi]] . \[CapitalOmega] . \[CapitalSigma] . Transpose[\[CapitalOmega]] . \[CapitalXi]];
Fmod[x_,p_] = FullSimplify[Chi[Join[
  {-p, -x},
  Flatten[Table[{-x/(K-2), 0}, {i, (K-2)}]],
  {p, -x}
]]]

(*Pi/Integrate[FPSmod[x,p],{x,-Infinity,Infinity},{p,-Infinity,Infinity}]*)
(*Collect[FPSmod[x,p]/Fmod[x,p],x^2,p^2]*)
(*ToMatlab[Fmod[x,p]];
ToMatlab[FPSmod[x,p]];*)
Fmod[x,p]*Exp[-1/2(x^2+p^2)]

Cohin[x_,p_] = Exp[-1/2(x^2+p^2)]
Cohout[x_,p_] = Cohin[x,p]*Fmod[x,p]
Fidelity1 = FullSimplify[Integrate[Cohin[x,p]*Cohout[-x,-p],{x,-Infinity,Infinity},{p,-Infinity,Infinity}]/Pi]



