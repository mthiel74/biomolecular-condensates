#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Concentration enhancement and reaction kinetics inside condensates.
   The origin-of-life argument: partition coefficient K amplifies
   bimolecular reaction rates by K^2.

   Output: docs/images/kin_rate_enhancement.png
           docs/images/kin_ribozyme.png
           docs/images/kin_prebiotic_space.png
*)

repoRoot = ParentDirectory @ DirectoryName[$InputFileName];
imgDir   = FileNameJoin[{repoRoot, "docs", "images"}];
If[!DirectoryQ[imgDir], CreateDirectory[imgDir, CreateIntermediateDirectories -> True]];

Print["Condensate kinetics figures..."];


(* ================================================================
   Figure 1: Rate enhancement vs partition coefficient
   ================================================================ *)

Export[FileNameJoin[{imgDir, "kin_rate_enhancement.png"}],
   LogLogPlot[
     {K, K^2, K^3},
     {K, 1, 1000},
     PlotStyle -> {
       Directive[Thick, Darker[Green]],
       Directive[Thick, Blue],
       Directive[Thick, Dashed, Red]},
     PlotLegends -> Placed[{
       "Unimolecular (A \[Rule] B): K",
       "Bimolecular (A + B \[Rule] C): K\:00b2",
       "Trimolecular: K\:00b3"}, {0.35, 0.85}],
     Frame -> True,
     FrameLabel -> {"Partition coefficient K",
        "Reaction rate enhancement factor"},
     PlotLabel -> Style["Rate enhancement inside biomolecular condensates", 14],
     GridLines -> {{100}, {}},
     GridLinesStyle -> Directive[Gray, Dashed],
     Epilog -> {
       Text[Style["Typical condensate\nK \[TildeTilde] 100", 10, Italic, Gray],
         {Log[100], Log[10]}, {0, -1.5}]},
     ImageSize -> 640,
     PlotRange -> {{1, 1000}, {1, 10^9}}],
   ImageResolution -> 200];
Print["  saved kin_rate_enhancement.png"];


(* ================================================================
   Figure 2: Ribozyme kinetics — replication vs degradation
   ================================================================ *)

Module[{Vmax = 1.0, Km = 10.0, kDeg = 0.05, Sfixed = 1.0,
        Krange, vRep, vDeg},

   Krange = Subdivide[1., 200., 300];

   vRep = Table[
     Module[{Seff = K Sfixed, Veff = K Vmax},
       Veff Seff / (Km + Seff)],
     {K, Krange}];

   vDeg = kDeg # Sfixed & /@ Krange;

   Export[FileNameJoin[{imgDir, "kin_ribozyme.png"}],
     Show[
       ListLinePlot[
         {Transpose[{Krange, vRep}], Transpose[{Krange, vDeg}]},
         PlotStyle -> {Directive[Thick, Blue], Directive[Thick, Dashed, Red]},
         PlotLegends -> {"Replication rate", "Degradation rate"},
         Filling -> {1 -> {{2}, Directive[Opacity[0.15], Green]}}],
       Frame -> True,
       FrameLabel -> {"Partition coefficient K", "Rate (\[Mu]M/min)"},
       PlotLabel -> Style["Replication vs degradation: the origin-of-life window", 14],
       GridLines -> Automatic, GridLinesStyle -> Directive[LightGray],
       ImageSize -> 640],
     ImageResolution -> 200];
   Print["  saved kin_ribozyme.png"];
];


(* ================================================================
   Figure 3: Prebiotic phase space (K vs concentration)
   ================================================================ *)

Module[{Vmax = 1.0, Km = 10.0, kDeg = 1.0,
        Krange, Crange, netRate},

   (* kDeg = 1.0 gives K_critical ~ 10 at low c (the biologically
      meaningful regime). At high c replication saturates while
      degradation grows linearly, creating a nontrivial boundary. *)
   Krange = Subdivide[1., 200., 150];
   Crange = 10^Subdivide[-3., 1., 150];

   netRate = Table[
     Module[{Seff = K c, Veff = K Vmax, rep, deg},
       rep = Veff Seff / (Km + Seff);
       deg = kDeg K c;
       rep - deg],
     {K, Krange}, {c, Crange}];

   Export[FileNameJoin[{imgDir, "kin_prebiotic_space.png"}],
     Show[
       ListContourPlot[
         Flatten[Table[
           {Krange[[i]], Log10[Crange[[j]]], netRate[[i, j]]},
           {i, Length[Krange]}, {j, Length[Crange]}], 1],
         Contours -> 20,
         ColorFunction -> "RedBlueTones",
         ContourLabels -> None,
         PlotLegends -> Placed[BarLegend[{"RedBlueTones", Automatic},
           LegendLabel -> "Net rate"], Right],
         FrameLabel -> {"Partition coefficient K",
            Row[{"log", Subscript["", "10"], "([Substrate]", Subscript["", "bulk"], " / \[Mu]M)"}]},
         PlotLabel -> Style["Net replication: replication \[Minus] degradation", 14]],
       ListContourPlot[
         Flatten[Table[
           {Krange[[i]], Log10[Crange[[j]]], netRate[[i, j]]},
           {i, Length[Krange]}, {j, Length[Crange]}], 1],
         Contours -> {0},
         ContourStyle -> Directive[Black, Thick],
         ContourLabels -> None],
       Frame -> True,
       ImageSize -> 640],
     ImageResolution -> 200];
   Print["  saved kin_prebiotic_space.png"];
];

Print["Kinetics figures complete."];
