#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Liquid \[Rule] gel \[Rule] aggregate transition relevant to
   neurodegenerative disease (FUS, TDP-43, tau).

   Output: docs/images/disease_transition.png
           docs/images/disease_phase_diagram.png
*)

repoRoot = ParentDirectory @ DirectoryName[$InputFileName];
imgDir   = FileNameJoin[{repoRoot, "docs", "images"}];
If[!DirectoryQ[imgDir], CreateDirectory[imgDir, CreateIntermediateDirectories -> True]];

Print["Disease transition figures..."];


(* ================================================================
   Three-state kinetic model: liquid \[Rule] gel \[Rule] aggregate
   d[liquid]/dt  = -kGel \[CenterDot] [liquid]
   d[gel]/dt     = kGel \[CenterDot] [liquid] - kAgg \[CenterDot] [gel]
   d[aggregate]/dt = kAgg \[CenterDot] [gel]
   ================================================================ *)

agingModel[kGel_, kAgg_, tMax_: 100.] :=
   NDSolve[{
     l'[t] == -kGel l[t],
     g'[t] == kGel l[t] - kAgg g[t],
     a'[t] == kAgg g[t],
     l[0] == 1, g[0] == 0, a[0] == 0},
     {l, g, a}, {t, 0, tMax}][[1]];


(* Three scenarios *)
solSlow     = agingModel[0.02, 0.005];
solModerate = agingModel[0.05, 0.01];
solFast     = agingModel[0.1, 0.02];

makePanelPlot[sol_, title_] := Plot[
   Evaluate[{l[t], g[t], a[t]} /. sol],
   {t, 0, 100},
   PlotStyle -> {
     Directive[Thick, Blue],
     Directive[Thick, Orange],
     Directive[Thick, Red]},
   Filling -> {1 -> Axis, 2 -> Axis, 3 -> Axis},
   FillingStyle -> {
     Directive[Opacity[0.15], LightBlue],
     Directive[Opacity[0.15], LightOrange],
     Directive[Opacity[0.15], LightRed]},
   PlotLegends -> Placed[{"Liquid", "Gel", "Aggregate"}, {0.8, 0.5}],
   Frame -> True,
   FrameLabel -> {"Time (a.u.)", "Fraction"},
   PlotLabel -> Style[title, 12],
   PlotRange -> {{0, 100}, {0, 1.05}},
   GridLines -> Automatic, GridLinesStyle -> Directive[LightGray],
   ImageSize -> 350];

Export[FileNameJoin[{imgDir, "disease_transition.png"}],
   Labeled[
     GraphicsGrid[{{
       makePanelPlot[solSlow, "Slow maturation\n(wild type)"],
       makePanelPlot[solModerate, "Moderate\n(risk variant)"],
       makePanelPlot[solFast, "Fast (disease mutant)\n(FUS/TDP-43)"]
     }}, ImageSize -> 1050, Spacings -> 15],
     Style["Condensate aging: liquid \[Rule] gel \[Rule] aggregate \
(neurodegeneration)", 14, Bold], Top],
   ImageResolution -> 200];
Print["  saved disease_transition.png"];


(* ================================================================
   Phase diagram of condensate material states
   ================================================================ *)

Export[FileNameJoin[{imgDir, "disease_phase_diagram.png"}],
   Module[{boundLL, boundGel, boundAgg, c},
     c = Subdivide[0., 10., 200];
     boundLL  = 3.0 + 2.0 Exp[-0.5 #] & /@ c;
     boundGel = 2.0 + 1.0 Exp[-0.3 #] & /@ c;
     boundAgg = 1.0 + 0.5 Exp[-0.2 #] & /@ c;

     Show[
       (* Regions *)
       ListLinePlot[{
         Transpose[{c, boundLL}],
         Transpose[{c, boundGel}],
         Transpose[{c, boundAgg}]},
         PlotStyle -> {
           Directive[Thick, Blue],
           Directive[Thick, Orange],
           Directive[Thick, Red]},
         Filling -> {
           1 -> {10, Directive[Opacity[0.1], LightBlue]},
           2 -> {{1}, Directive[Opacity[0.1], LightOrange]},
           3 -> {{2}, Directive[Opacity[0.1], LightRed]}},
         PlotLegends -> {"Liquid condensate", "Gel transition", "Solid aggregate"}],

       (* Labels *)
       Graphics[{
         Text[Style["Mixed\n(one phase)", 13, Bold, Blue], {5, 8}],
         Text[Style["Liquid condensate", 12, Bold, Darker[Orange]], {7, 4.2}],
         Text[Style["Gel", 12, Bold, Darker[Red, 0.3]], {7, 2.5}],
         Text[Style["Solid aggregate", 11, Bold, Darker[Red]], {7, 0.7}],
         {Red, Thick, Arrowheads[0.04],
          Arrow[{{2, 4.5}, {4, 3}}]},
         Text[Style["Disease mutations\nshift boundaries", 10, Italic, Red],
           {1.5, 5.2}]
       }],

       Frame -> True,
       FrameLabel -> {"Protein concentration (a.u.)",
          "Temperature / Interaction strength (a.u.)"},
       PlotLabel -> Style["Condensate material state phase diagram", 14],
       PlotRange -> {{0, 10}, {0, 10}},
       GridLines -> Automatic, GridLinesStyle -> Directive[LightGray],
       ImageSize -> 640]
   ],
   ImageResolution -> 200];
Print["  saved disease_phase_diagram.png"];

Print["Disease transition figures complete."];
