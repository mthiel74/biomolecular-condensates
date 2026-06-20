#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Cahn-Hilliard equation: spinodal decomposition.

   We use Wolfram's NDSolve for a 1D Cahn-Hilliard PDE, which
   WL handles natively with Method of Lines. The 2D simulation
   is in the Python implementation; here we focus on producing
   stunning Wolfram-native visualisations of the 1D dynamics and
   the structure factor evolution.

   Output: docs/images/ch_1d_evolution.png
           docs/images/ch_free_energy_landscape.png
*)

repoRoot = ParentDirectory @ DirectoryName[$InputFileName];
imgDir   = FileNameJoin[{repoRoot, "docs", "images"}];
If[!DirectoryQ[imgDir], CreateDirectory[imgDir, CreateIntermediateDirectories -> True]];

Print["Cahn-Hilliard figures..."];


(* ================================================================
   1D Cahn-Hilliard via NDSolve

   \[PartialD]\[Phi]/\[PartialD]t = M \[PartialD]^2/\[PartialD]x^2
       [ \[PartialD]f/\[PartialD]\[Phi] - \[Kappa] \[PartialD]^2\[Phi]/\[PartialD]x^2 ]

   with Flory-Huggins f(\[Phi]) and periodic BCs.
   ================================================================ *)

Module[{L = 20, \[Chi] = 1.5, N1 = 100, N2 = 1, \[Kappa] = 0.5,
        mob = 1.0, \[Phi]0 = 0.3, noise = 0.02, sol, eps = 10.^-6},

   (* Chemical potential from FH *)
   mu[\[Phi]_] := (Log[Clip[\[Phi], {eps, 1 - eps}]] + 1)/N1
                - (Log[Clip[1 - \[Phi], {eps, 1 - eps}]] + 1)/N2
                + \[Chi] (1 - 2 \[Phi]);

   (* Initial condition: mean + small perturbation *)
   SeedRandom[42];
   nPts = 200;
   initData = Table[\[Phi]0 + noise Sin[2 Pi 3 x/L] +
      noise Sin[2 Pi 7 x/L] + noise/2 Sin[2 Pi 13 x/L],
      {x, Subdivide[0., L, nPts]}];
   initFunc = Interpolation[
      Transpose[{Subdivide[0., L, nPts], initData}],
      InterpolationOrder -> 3, PeriodicInterpolation -> True];

   Print["  Solving 1D Cahn-Hilliard PDE via NDSolve..."];
   sol = Quiet @ NDSolve[{
      D[\[Phi][x, t], t] == mob D[mu[\[Phi][x, t]] - \[Kappa] D[\[Phi][x, t], x, x], x, x],
      \[Phi][x, 0] == initFunc[x],
      \[Phi][0, t] == \[Phi][L, t],
      Derivative[1, 0][\[Phi]][0, t] == Derivative[1, 0][\[Phi]][L, t]},
      \[Phi], {x, 0, L}, {t, 0, 50},
      Method -> {"MethodOfLines",
        "SpatialDiscretization" -> {"TensorProductGrid",
          "MinPoints" -> 200, "MaxPoints" -> 400}},
      MaxSteps -> 100000,
      AccuracyGoal -> 4, PrecisionGoal -> 4];

   If[MatchQ[sol, {{_Rule}}],
     Print["  NDSolve succeeded."];

     (* Figure 1: space-time evolution *)
     Export[FileNameJoin[{imgDir, "ch_1d_evolution.png"}],
       GraphicsColumn[{
         DensityPlot[Evaluate[\[Phi][x, t] /. sol[[1]]],
           {x, 0, L}, {t, 0, 50},
           PlotRange -> {0, 0.8},
           ColorFunction -> "TemperatureMap",
           PlotLegends -> Placed[BarLegend[{"TemperatureMap", {0, 0.8}},
             LegendLabel -> "\[Phi]"], Right],
           FrameLabel -> {"Position x", "Time t"},
           PlotLabel -> Style["1D spinodal decomposition (space-time)", 14],
           ImageSize -> 600, AspectRatio -> 0.4],
         Plot[
           Evaluate[Table[\[Phi][x, t] /. sol[[1]], {t, {0, 2, 10, 30, 50}}]],
           {x, 0, L},
           PlotStyle -> (Directive[Thick, #] & /@ {
             Darker[Blue, 0.3], Blue, Darker[Green], Orange, Red}),
           PlotLegends -> {"t = 0", "t = 2", "t = 10", "t = 30", "t = 50"},
           Frame -> True,
           FrameLabel -> {"Position x", "Volume fraction \[Phi]"},
           PlotLabel -> Style["Composition profiles at selected times", 13],
           PlotRange -> {0, 1},
           GridLines -> Automatic, GridLinesStyle -> Directive[LightGray],
           ImageSize -> 600]
       }, Spacings -> 15],
       ImageResolution -> 200];
     Print["  saved ch_1d_evolution.png"],

     Print["  WARN: NDSolve did not converge. Skipping CH figures."]
   ];
];


(* ================================================================
   Free energy landscape figure (analytical, no simulation needed)
   ================================================================ *)

Module[{N1 = 100, N2 = 1},
   Export[FileNameJoin[{imgDir, "ch_free_energy_landscape.png"}],
     Show[
       Plot3D[
         \[Phi] Log[Clip[\[Phi], {0.001, 0.999}]]/N1 +
         (1 - \[Phi]) Log[Clip[1 - \[Phi], {0.001, 0.999}]]/N2 +
         \[Chi] \[Phi] (1 - \[Phi]),
         {\[Phi], 0.01, 0.99}, {\[Chi], 0.3, 2.0},
         PlotRange -> All,
         ColorFunction -> "TemperatureMap",
         AxesLabel -> {"\[Phi]", "\[Chi]", "f(\[Phi])"},
         PlotLabel -> Style["Flory\[Dash]Huggins free energy landscape", 14],
         MeshFunctions -> {#3 &},
         MeshStyle -> Directive[Gray, Thin],
         ImageSize -> 640,
         ViewPoint -> {-2.5, -1.5, 1.5},
         Lighting -> "Neutral"],
       ImageResolution -> 200]
   ];
   Print["  saved ch_free_energy_landscape.png"];
];

Print["Cahn-Hilliard figures complete."];
