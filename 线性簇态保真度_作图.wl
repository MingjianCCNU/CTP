(* ::Package:: *)

Clear["Global`*"];
SetDirectory[NotebookDirectory[]];
Import["ToMatlab.m", "Package"];

(*\:8f93\:5165\:6a21\:7684r*)

$Assumptions = {r \[Element] Reals};
r = 0;
Nmode = 1;

IN = IdentityMatrix[Nmode];

J0 = {{1,I},{1,-I}}/2;
Z0 = {{1,0},{0,-1}};
\[CapitalOmega]0 = {{0,1},{-1,0}};
J = KroneckerProduct[IN,J0];
Z = KroneckerProduct[IN,Z0];
\[CapitalOmega] = KroneckerProduct[IN,\[CapitalOmega]0];

\[CapitalSigma] = {{Exp[-2r],0},{0,Exp[2r]}}

nsub = 2;
CF[a_,ac_] = Simplify[Exp[-{{ac,a}} . Z . J . \[CapitalSigma] . ConjugateTranspose[J] . Z . {{a},{ac}}/2]];
CFxp[x_,p_] = Simplify[CF[x+I*p,x-I*p]];

CFchannel1[x_,p_] = E^(-0.33462247687924274` p^2-0.33462247687924374` x^2);
CFchannel2[x_,p_] = E^(-0.23474492661905244` p^2-0.23474492661905189` x^2);
CFchannel3[x_,p_] = E^(-0.4746000270036079` p^2-0.09476737649468603` x^2);

CFout1[x_,p_] = CFxp[x,p]*CFchannel1[x,p];
CFout2[x_,p_] = CFxp[x,p]*CFchannel2[x,p];
CFout3[x_,p_] = CFxp[x,p]*CFchannel3[x,p];

(*r14\:5927\:60c5\:51b5*)
Fidelity1 = ToMatlab[Integrate[CFxp[x,p]E^(-0.33462247687924274` p^2-0.33462247687924374` x^2)*CFxp[-x,-p],{x,-Infinity,Infinity},{p,-Infinity,Infinity}]/Pi]
(*r23\:5927\:60c5\:51b5*)
Fidelity2 = ToMatlab[Integrate[CFxp[x,p]E^(-0.23474492661905244` p^2-0.23474492661905189` x^2)*CFxp[-x,-p],{x,-Infinity,Infinity},{p,-Infinity,Infinity}]/Pi]
(*r24\:5927\:60c5\:51b5*)
Fidelity3 = ToMatlab[Integrate[CFxp[x,p]E^(-0.4746000270036079` p^2-0.09476737649468603` x^2)*CFxp[-x,-p],{x,-Infinity,Infinity},{p,-Infinity,Infinity}]/Pi]

(*Plot axes*)
RangePlot = 3;
P1 = Plot3D[0,{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
PlotTheme->"Minimal",
PlotPoints->300,PlotRange->All,
PlotStyle->None,
MeshStyle->None,
Boxed->False,
PlotRangePadding->None,
RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2],
AxesOrigin->{0,0,0},Axes->{True, True, False},
AxesLabel->{x,p},
AxesStyle -> Directive[Thickness[0.01]],
ViewPoint-> 10{3,5,5},
LabelStyle -> Directive[FontSize -> 40],
ImagePadding->None]
Export["C:\\Users\\PS\\Desktop\\axes.png",P1,ImageResolution -> 300, Background->None,ImageMargins->0];

(*img = Import["C:\\Users\\PS\\Desktop\\axes.png"];
croppedImg = ImageCrop[img];

Export["C:\\Users\\PS\\Desktop\\axes_cropped.png", croppedImg];*)


(*Plot input mode*)
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
Axes->False,
ViewPoint-> 10{3,5,5},
LabelStyle -> Directive[FontSize -> 40],ImagePadding->None]
Export["C:\\Users\\PS\\Desktop\\sqzin.png",P1,ImageResolution -> 300, Background->None,ImageMargins->0];

(*Plot output mode1*)
RangePlot = 3;
DarkBlue = "#0089ba";
LightGrey = "#b0a8b9";
Orangee = "#ff6f91";
Yelloww = "#ffc75f";
P1 = Plot3D[Re[CFout1[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
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
Axes->False,
ViewPoint-> 10{3,5,5}]
Export["C:\\Users\\PS\\Desktop\\sqzout1.png",P1,ImageResolution -> 300, Background->None];

(*Plot output mode2*)
RangePlot = 3;
DarkBlue = "#0089ba";
LightGrey = "#b0a8b9";
Orangee = "#ff6f91";
Yelloww = "#ffc75f";
P1 = Plot3D[Re[CFout2[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
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
Axes->False,
ViewPoint-> 10{3,5,5}]
Export["C:\\Users\\PS\\Desktop\\sqzout2.png",P1,ImageResolution -> 300, Background->None];

(*Plot output mode*)
RangePlot = 3;
DarkBlue = "#0089ba";
LightGrey = "#b0a8b9";
Orangee = "#ff6f91";
Yelloww = "#ffc75f";
P1 = Plot3D[Re[CFout3[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
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
Axes->False,
ViewPoint-> 10{3,5,5}]
Export["C:\\Users\\PS\\Desktop\\sqzout3.png",P1,ImageResolution -> 300, Background->None];

(*Plot channel response 1*)
RangePlot = 3;
DarkBlue = "#0089ba";
LightGrey = "#b0a8b9";
Orangee = "#ff6f91";
Yelloww = "#ffc75f";
P1 = Plot3D[Re[CFchannel1[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
PlotTheme->"Minimal",
PlotPoints->300,PlotRange->All,
PlotStyle->Directive[Opacity[0.99]],
Lighting -> Automatic,
AxesEdge -> Automatic,
MeshFunctions->{#3 &},Mesh->{{1}},MeshStyle->Dashed,
(*ColorFunction -> "BlueGreenYellow",*)
RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2],
ColorFunction -> (Blend[{RGBColor[LightGrey], RGBColor[DarkBlue]}, #3] &),
AxesLabel->{x,p,Subscript[\[Chi], coh]},
Boxed->False,PlotRangePadding->None,
Axes->False,
ViewPoint-> 10{3,5,5}]
Export["C:\\Users\\PS\\Desktop\\ch1.png",P1,ImageResolution -> 300, Background->None];

(*Plot channel response 2*)
RangePlot = 3;
DarkBlue = "#0089ba";
LightGrey = "#b0a8b9";
Orangee = "#ff6f91";
Yelloww = "#ffc75f";
P1 = Plot3D[Re[CFchannel2[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
PlotTheme->"Minimal",
PlotPoints->300,PlotRange->All,
PlotStyle->Directive[Opacity[0.99]],
Lighting -> Automatic,
AxesEdge -> Automatic,
MeshFunctions->{#3 &},Mesh->{{1}},MeshStyle->Dashed,
(*ColorFunction -> "BlueGreenYellow",*)
RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2],
ColorFunction -> (Blend[{RGBColor[LightGrey], RGBColor[DarkBlue]}, #3] &),
AxesLabel->{x,p,Subscript[\[Chi], coh]},
Boxed->False,PlotRangePadding->None,
Axes->False,
ViewPoint-> 10{3,5,5}]
Export["C:\\Users\\PS\\Desktop\\ch2.png",P1,ImageResolution -> 300, Background->None];

(*Plot channel response 3*)
RangePlot = 3;
DarkBlue = "#0089ba";
LightGrey = "#b0a8b9";
Orangee = "#ff6f91";
Yelloww = "#ffc75f";
P1 = Plot3D[Re[CFchannel3[x,p]],{x,-RangePlot,RangePlot},{p,-RangePlot,RangePlot},
PlotTheme->"Minimal",
PlotPoints->300,PlotRange->All,
PlotStyle->Directive[Opacity[0.99]],
Lighting -> Automatic,
AxesEdge -> Automatic,
MeshFunctions->{#3 &},Mesh->{{1}},MeshStyle->Dashed,
(*ColorFunction -> "BlueGreenYellow",*)
RegionFunction -> Function[{x, y, z}, x^2 + y^2 <= RangePlot^2],
ColorFunction -> (Blend[{RGBColor[LightGrey], RGBColor[DarkBlue]}, #3] &),
AxesLabel->{x,p,Subscript[\[Chi], coh]},
Boxed->False,PlotRangePadding->None,
Axes->False,
ViewPoint-> 10{3,5,5}]
Export["C:\\Users\\PS\\Desktop\\ch3.png",P1,ImageResolution -> 300, Background->None];


ContourPlot3D[
  x^2 + y^2 + z^2 == -1,  (* \:65e0\:89e3\:7684\:65b9\:7a0b\:ff0c\:4e0d\:751f\:6210\:5b9e\:9645\:56fe\:5f62 *)
  {x, -1, 1}, {y, -1, 1}, {z, -1, 1},
  Axes -> True,
  AxesLabel -> {"X", "Y", "Z"},
  Boxed -> True,
  ContourStyle -> None    (* \:4e0d\:663e\:793a\:4efb\:4f55\:8f6e\:5ed3 *)
]
