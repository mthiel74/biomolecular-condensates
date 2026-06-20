#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Flory-Huggins thermodynamics: free energy, spinodal, binodal,
   critical point, and phase diagrams for biomolecular condensates.

   Output: docs/images/fh_free_energy.png
           docs/images/fh_phase_diagram.png
           docs/images/fh_chain_length.png
*)

repoRoot = ParentDirectory @ DirectoryName[$InputFileName];
imgDir   = FileNameJoin[{repoRoot, "docs", "images"}];
If[!DirectoryQ[imgDir], CreateDirectory[imgDir, CreateIntermediateDirectories -> True]];

(* ================================================================
   Flory-Huggins free energy density
   f(\[Phi]) = \[Phi] ln(\[Phi])/N1 + (1-\[Phi]) ln(1-\[Phi])/N2 + \[Chi] \[Phi](1-\[Phi])
   ================================================================ *)

fFH[\[Phi]_, \[Chi]_, N1_, N2_] :=
   \[Phi] Log[\[Phi]]/N1 + (1 - \[Phi]) Log[1 - \[Phi]]/N2 + \[Chi] \[Phi] (1 - \[Phi]);

muFH[\[Phi]_, \[Chi]_, N1_, N2_] :=
   (Log[\[Phi]] + 1)/N1 - (Log[1 - \[Phi]] + 1)/N2 + \[Chi] (1 - 2 \[Phi]);

d2fFH[\[Phi]_, \[Chi]_, N1_, N2_] :=
   1/(N1 \[Phi]) + 1/(N2 (1 - \[Phi])) - 2 \[Chi];

piOsmotic[\[Phi]_, \[Chi]_, N1_, N2_] :=
   fFH[\[Phi], \[Chi], N1, N2] - \[Phi] muFH[\[Phi], \[Chi], N1, N2];

(* Analytic critical point *)
criticalPoint[N1_, N2_] := Module[{sqr = Sqrt[N2/N1]},
   <|"\[Phi]c" -> sqr/(1 + sqr),
     "\[Chi]c" -> (1/Sqrt[N1] + 1/Sqrt[N2])^2/2|>];

(* Spinodal curve: locus of d2f/d\[Phi]2 = 0 *)
spinodalCurve[N1_, N2_, npts_: 300] := Module[
   {cp = criticalPoint[N1, N2], \[Phi]vals, \[Chi]vals},
   \[Phi]vals = Subdivide[0.001, 0.999, npts];
   \[Chi]vals = (1/(N1 #) + 1/(N2 (1 - #)))/2 & /@ \[Phi]vals;
   Select[Transpose[{\[Phi]vals, \[Chi]vals}],
      #[[2]] >= cp["\[Chi]c"] &]
];

(* Binodal via common tangent *)
binodalCurve[N1_, N2_, npts_: 60] := Module[
   {cp = criticalPoint[N1, N2], \[Chi]max, \[Chi]vals, results = {}},
   \[Chi]max = cp["\[Chi]c"] * 3;
   \[Chi]vals = Subdivide[cp["\[Chi]c"] + 0.001, \[Chi]max, npts];
   Do[
     Quiet @ Check[
       Module[{sol, \[Phi]L, \[Phi]R},
         sol = FindRoot[
           {muFH[\[Phi]L, \[Chi], N1, N2] == muFH[\[Phi]R, \[Chi], N1, N2],
            piOsmotic[\[Phi]L, \[Chi], N1, N2] == piOsmotic[\[Phi]R, \[Chi], N1, N2]},
           {\[Phi]L, cp["\[Phi]c"]/3, 0.001, cp["\[Phi]c"]},
           {\[Phi]R, 1 - cp["\[Phi]c"]/3, cp["\[Phi]c"], 0.999}];
         If[0 < (\[Phi]L /. sol) < (\[Phi]R /. sol) < 1,
           AppendTo[results, {\[Chi], \[Phi]L /. sol, \[Phi]R /. sol}]]],
       Null],
     {\[Chi], \[Chi]vals}];
   results
];


(* ================================================================
   Figure 1: Free energy curves at several \[Chi]
   ================================================================ *)

Print["Generating Flory-Huggins figures..."];

Module[{N1 = 100, N2 = 1, cp, \[Chi]vals, \[Phi]range, plots},
   cp = criticalPoint[N1, N2];
   \[Chi]vals = {0.3, 0.5, 1.0, 1.5, 2.5} cp["\[Chi]c"];
   \[Phi]range = Subdivide[0.001, 0.999, 500];

   plots = Table[
     ListLinePlot[
       Transpose[{\[Phi]range, fFH[#, \[Chi], N1, N2] & /@ \[Phi]range}],
       PlotStyle -> {Thick},
       PlotLegends -> None],
     {\[Chi], \[Chi]vals}];

   Export[FileNameJoin[{imgDir, "fh_free_energy.png"}],
     Plot[
       Evaluate[Table[fFH[\[Phi], \[Chi], N1, N2], {\[Chi], \[Chi]vals}]],
       {\[Phi], 0.001, 0.999},
       PlotStyle -> (Directive[Thick, #] & /@ {
          Darker[Blue, 0.3], Blue, Darker[Green],
          Darker[Orange], Darker[Red]}),
       PlotLegends -> (Row[{"\[Chi] = ", NumberForm[#, {3, 3}]}] & /@ \[Chi]vals),
       AxesLabel -> {"\[Phi]", "f(\[Phi])/kT"},
       PlotLabel -> Style[Row[{"Flory\[Dash]Huggins free energy  (", Subscript["N", "1"],
          " = ", N1, ", ", Subscript["N", "2"], " = ", N2, ")"}], 14],
       PlotRange -> {Automatic, {-0.15, 0.02}},
       GridLines -> {{}, {0}},
       GridLinesStyle -> Directive[Gray, Dashed],
       ImageSize -> 640, Frame -> True,
       FrameLabel -> {"Volume fraction \[Phi]",
          "Free energy density  f(\[Phi])/kT"}],
     ImageResolution -> 200];
   Print["  saved fh_free_energy.png"];
];


(* ================================================================
   Figure 2: Phase diagram (\[Chi] vs \[Phi])
   ================================================================ *)

Module[{N1 = 100, N2 = 1, cp, sp, bn, pSpin, pBin, pCrit},
   cp = criticalPoint[N1, N2];
   sp = spinodalCurve[N1, N2];
   bn = binodalCurve[N1, N2];

   pSpin = ListLinePlot[sp,
      PlotStyle -> Directive[Red, Dashed, Thick],
      PlotLegends -> {"Spinodal"}];

   pBin = If[Length[bn] > 0,
     ListLinePlot[
       {bn[[All, {2, 1}]], bn[[All, {3, 1}]]},
       PlotStyle -> Directive[Blue, Thick],
       PlotLegends -> {"Binodal", None}],
     Graphics[]];

   pCrit = ListPlot[{{cp["\[Phi]c"], cp["\[Chi]c"]}},
      PlotStyle -> {Black, PointSize[Large]},
      PlotLegends -> {"Critical point"}];

   Export[FileNameJoin[{imgDir, "fh_phase_diagram.png"}],
     Show[pBin, pSpin, pCrit,
       PlotRange -> {{0, 1}, {0, cp["\[Chi]c"] * 3}},
       Frame -> True,
       FrameLabel -> {"Volume fraction \[Phi]",
          "Flory\[Dash]Huggins parameter \[Chi]  (\[Proportional] 1/T)"},
       PlotLabel -> Style[Row[{"Phase diagram  (", Subscript["N", "1"],
          " = 100, ", Subscript["N", "2"], " = 1)"}], 14],
       ImageSize -> 640,
       GridLines -> Automatic, GridLinesStyle -> Directive[LightGray]],
     ImageResolution -> 200];
   Print["  saved fh_phase_diagram.png"];
];


(* ================================================================
   Figure 3: Chain length effect on phase diagram
   ================================================================ *)

Module[{Nvals = {1, 5, 20, 100, 500}, colors, plots},
   colors = ColorData["DarkRainbow"] /@ Subdivide[0, 1, Length[Nvals] - 1];

   plots = MapThread[
     Module[{N1 = #1, cp, sp, T\[Phi]},
       cp = criticalPoint[N1, 1];
       sp = spinodalCurve[N1, 1, 200];
       T\[Phi] = {#[[1]], 1/#[[2]]} & /@ sp;
       ListLinePlot[T\[Phi],
         PlotStyle -> Directive[Thick, #2],
         PlotLegends -> {Row[{"N = ", N1}]}]
     ] &,
     {Nvals, colors}];

   Export[FileNameJoin[{imgDir, "fh_chain_length.png"}],
     Show[plots,
       PlotRange -> {{0, 1}, {0, Automatic}},
       Frame -> True,
       FrameLabel -> {"Volume fraction \[Phi]", "Reduced temperature T*"},
       PlotLabel -> Style["Chain length widens the two-phase region", 14],
       ImageSize -> 640,
       GridLines -> Automatic, GridLinesStyle -> Directive[LightGray]],
     ImageResolution -> 200];
   Print["  saved fh_chain_length.png"];
];

Print["Phase diagram figures complete."];
