#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Protocell dynamics: Lifshitz-Slyozov coarsening, proto-Darwinian
   selection, and growth/division/competition.

   Output: docs/images/proto_ostwald.png
           docs/images/proto_selection.png
           docs/images/proto_growth_division.png
*)

repoRoot = ParentDirectory @ DirectoryName[$InputFileName];
imgDir   = FileNameJoin[{repoRoot, "docs", "images"}];
If[!DirectoryQ[imgDir], CreateDirectory[imgDir, CreateIntermediateDirectories -> True]];

Print["Protocell dynamics figures..."];


(* ================================================================
   Lifshitz-Slyozov-Wagner coarsening simulation
   ================================================================ *)

ostwaldRipening[nDroplets_, tMax_, dt_, seed_] := Module[
   {R, Rstar, dRdt, alive, t = 0., lc = 1.0, D0 = 1.0,
    history = {}, rng},

   SeedRandom[seed];
   R = Exp[RandomVariate[NormalDistribution[0, 0.5], nDroplets]];
   R = Clip[R, {0.01, \[Infinity]}];
   AppendTo[history, {0., R}];

   While[t < tMax && Count[R, x_ /; x > 0.01] > 2,
     alive = UnitStep[R - 0.01];
     Rstar = Total[R alive] / Max[1, Total[alive]];

     dRdt = alive D0 lc (1/Rstar - 1/Clip[R, {0.01, \[Infinity]}]) /
            Clip[R, {0.01, \[Infinity]}];
     R += dRdt dt;
     R = Clip[R, {0, \[Infinity]}];
     R = R UnitStep[R - 0.01]; (* dissolve tiny *)

     t += dt;
     If[Mod[Round[t/dt], Round[1/dt]] == 0,
       AppendTo[history, {t, R}]];
   ];
   history
];


Print["  Running Ostwald ripening (100 droplets)..."];
history = ostwaldRipening[100, 100., 0.01, 42];

(* Mean and max radius over time *)
radiusData = {#[[1]],
   Mean[Select[#[[2]], # > 0.01 &]],
   Max[#[[2]]],
   Count[#[[2]], x_ /; x > 0.01]} & /@ history;

Export[FileNameJoin[{imgDir, "proto_ostwald.png"}],
   GraphicsGrid[{{
     ListLinePlot[
       {radiusData[[All, {1, 2}]], radiusData[[All, {1, 3}]]},
       PlotStyle -> {Directive[Thick, Blue], Directive[Thick, Dashed, Red]},
       PlotLegends -> {"Mean R", "Max R"},
       Frame -> True,
       FrameLabel -> {"Time", "Radius"},
       PlotLabel -> Style["Droplet coarsening", 13],
       GridLines -> Automatic, GridLinesStyle -> Directive[LightGray],
       ImageSize -> 400],
     ListLinePlot[radiusData[[All, {1, 4}]],
       PlotStyle -> Directive[Thick, Darker[Green]],
       Frame -> True,
       FrameLabel -> {"Time", "Surviving droplets"},
       PlotLabel -> Style["Ostwald ripening", 13],
       GridLines -> Automatic, GridLinesStyle -> Directive[LightGray],
       ImageSize -> 400]
   }}, ImageSize -> 800],
   ImageResolution -> 200];
Print["  saved proto_ostwald.png"];


(* ================================================================
   Proto-Darwinian selection simulation
   ================================================================ *)

protocellSimulation[nInit_, tMax_, dt_, cBulk_, seed_] := Module[
   {cells, t = 0., history = {}, maxPop = 300, rng, newCells, fitness},

   SeedRandom[seed];
   (* Each cell: {radius, K, kRep, age, generation} *)
   cells = Table[
     {RandomReal[{0.5, 2.0}],
      RandomReal[{10, 200}],
      RandomReal[{0.01, 0.2}],
      0., 0},
     {nInit}];

   While[t < tMax && Length[cells] > 0,
     (* Growth *)
     cells = Map[
       Module[{r = #[[1]], K = #[[2]], kR = #[[3]], age = #[[4]], gen = #[[5]],
               cEff, netRate, dV, vol},
         cEff = K cBulk;
         netRate = kR cEff - 0.05 cEff;
         vol = 4/3 Pi r^3;
         If[netRate > 0, vol += netRate vol dt 0.01];
         r = (3 vol / (4 Pi))^(1/3);
         {r, K, kR, age + dt, gen}] &,
       cells];

     (* Division *)
     newCells = {};
     Do[
       If[cell[[1]] > 3.0,
         Module[{rNew = cell[[1]] / 2^(1/3), K = cell[[2]], kR = cell[[3]],
                 gen = cell[[5]], Kn, kn},
           Kn = 1 + 0.05 RandomReal[NormalDistribution[]];
           kn = 1 + 0.05 RandomReal[NormalDistribution[]];
           AppendTo[newCells, {rNew, K Kn, kR kn, 0., gen + 1}];
           AppendTo[newCells, {rNew, K/Kn, kR/kn, 0., gen + 1}]],
         AppendTo[newCells, cell]],
       {cell, cells}];
     cells = newCells;

     (* Death *)
     cells = Select[cells, #[[1]] > 0.1 &];

     (* Resource limit *)
     If[Length[cells] > maxPop,
       fitness = #[[2]] #[[3]] & /@ cells;
       cells = cells[[Ordering[fitness, -maxPop]]]];

     (* Record *)
     If[Mod[Round[t/dt], 10] == 0 && Length[cells] > 0,
       AppendTo[history,
         {t, Length[cells],
          Mean[cells[[All, 2]]],
          Mean[cells[[All, 3]]],
          Max[cells[[All, 5]]]}]];
     t += dt];

   history
];


Print["  Running protocell selection (50 cells, 500 time units)..."];
protoHist = protocellSimulation[50, 500., 0.5, 0.1, 42];

If[Length[protoHist] > 0,
   Export[FileNameJoin[{imgDir, "proto_selection.png"}],
     GraphicsGrid[{{
       ListLinePlot[protoHist[[All, {1, 3}]],
         PlotStyle -> Directive[Thick, Red],
         Frame -> True,
         FrameLabel -> {"Time", "Mean K"},
         PlotLabel -> Style["Selection for higher K", 13],
         GridLines -> Automatic, GridLinesStyle -> Directive[LightGray],
         ImageSize -> 350],
       ListLinePlot[protoHist[[All, {1, 4}]],
         PlotStyle -> Directive[Thick, Darker[Green]],
         Frame -> True,
         FrameLabel -> {"Time", "Mean replication rate"},
         PlotLabel -> Style["Selection for faster replication", 13],
         GridLines -> Automatic, GridLinesStyle -> Directive[LightGray],
         ImageSize -> 350],
       ListLinePlot[protoHist[[All, {1, 5}]],
         PlotStyle -> Directive[Thick, Purple],
         Frame -> True,
         FrameLabel -> {"Time", "Max generation"},
         PlotLabel -> Style["Lineage depth", 13],
         GridLines -> Automatic, GridLinesStyle -> Directive[LightGray],
         ImageSize -> 350]
     }}, ImageSize -> 1000],
     ImageResolution -> 200];
   Print["  saved proto_selection.png"],
   Print["  WARN: no protocell history to plot"]
];


(* ================================================================
   Growth-division schematic
   ================================================================ *)

Export[FileNameJoin[{imgDir, "proto_growth_division.png"}],
   Graphics[{
     (* Phase 1: small droplet *)
     {EdgeForm[Directive[Thick, Blue]], FaceForm[Directive[Opacity[0.3], LightBlue]],
      Disk[{0, 0}, 0.5]},
     Text[Style["Small\ncondensate", 10], {0, -0.9}],

     (* Arrow *)
     {Thick, Darker[Green], Arrowheads[0.03],
      Arrow[{{0.7, 0}, {2.3, 0}}]},
     Text[Style["Growth\n(accretion)", 9, Italic, Darker[Green]], {1.5, 0.3}],

     (* Phase 2: large droplet *)
     {EdgeForm[Directive[Thick, Blue]], FaceForm[Directive[Opacity[0.3], LightBlue]],
      Disk[{3.5, 0}, 1.0]},
     Text[Style["Large\ncondensate", 10], {3.5, -1.4}],

     (* Arrow *)
     {Thick, Red, Arrowheads[0.03],
      Arrow[{{4.7, 0}, {6.3, 0}}]},
     Text[Style["Division", 9, Italic, Red], {5.5, 0.3}],

     (* Phase 3: two daughters *)
     {EdgeForm[Directive[Thick, Blue]], FaceForm[Directive[Opacity[0.3], LightBlue]],
      Disk[{7.5, 0.5}, 0.6]},
     {EdgeForm[Directive[Thick, Blue]], FaceForm[Directive[Opacity[0.3], LightBlue]],
      Disk[{7.5, -0.5}, 0.6]},
     Text[Style["Daughters\n(with variation)", 10], {7.5, -1.4}],

     (* Competition arrow *)
     {Thick, Orange, Arrowheads[0.03],
      Arrow[{{8.3, 0.5}, {9.8, 0.5}}]},
     {Thick, Gray, Dashed, Arrowheads[0.03],
      Arrow[{{8.3, -0.5}, {9.0, -1.2}}]},
     Text[Style["Fitter\nsurvives", 9, Italic, Orange], {9.3, 0.8}],
     Text[Style["Weaker\ndissolves", 9, Italic, Gray], {9.5, -1.0}],

     (* Winner *)
     {EdgeForm[Directive[Thick, Darker[Green]]],
      FaceForm[Directive[Opacity[0.3], LightGreen]],
      Disk[{10.5, 0.5}, 0.55]}
   },
   PlotRange -> {{-1.2, 11.5}, {-2, 1.8}},
   ImageSize -> 800,
   PlotLabel -> Style["Protocell life cycle: growth \[RightArrow] \
division \[RightArrow] selection", 14]],
   ImageResolution -> 200];
Print["  saved proto_growth_division.png"];

Print["Protocell figures complete."];
