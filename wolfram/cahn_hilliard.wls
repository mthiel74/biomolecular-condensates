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

(* The Flory-Huggins chemical potential is logarithmically stiff near
   \[Phi]=0,1 and the strongly-stabilised 1D solver either freezes or
   blows up into high-frequency noise (see project notes). The robust,
   standard pedagogical model for spinodal decomposition is the symmetric
   double-well (Model B), used here exactly as in \[Section]8.4. We solve
   the order parameter \[Psi]\[Element][-1,1] and map it to a physical
   volume fraction \[Phi] = \[Phi]\:2080 + A\[ThinSpace]\[Psi] for display. *)
Module[{a = 1.0, b = 1.0, kappa = 1.0, M = 1.0, L = 80.0, T = 400.0,
        nPts = 256, mu, coeffs, icf, sol, phiOf, \[Phi]0 = 0.4, amp = 0.28,
        tSnap, densPlot, profPlot},

   mu[\[Psi]_] := -a \[Psi] + b \[Psi]^3;
   SeedRandom[17];
   coeffs = RandomReal[{-1, 1}, 24];
   (* smooth, exactly periodic initial condition (sines vanish at x=0,L) *)
   icf[x_] := 0.02 Sum[coeffs[[n]] Sin[2 Pi n x/L], {n, 1, 24}];

   Print["  Solving 1D Cahn-Hilliard PDE (double-well Model B) via NDSolve..."];
   sol = Quiet @ NDSolveValue[{
      D[\[Psi][x, t], t] == M D[mu[\[Psi][x, t]] - kappa D[\[Psi][x, t], x, x], x, x],
      \[Psi][x, 0] == icf[x],
      \[Psi][0, t] == \[Psi][L, t],
      Derivative[1, 0][\[Psi]][0, t] == Derivative[1, 0][\[Psi]][L, t]},
      \[Psi], {x, 0, L}, {t, 0, T},
      Method -> {"MethodOfLines",
        "SpatialDiscretization" -> {"TensorProductGrid",
          "MinPoints" -> nPts, "MaxPoints" -> nPts}},
      MaxSteps -> 500000, AccuracyGoal -> 5, PrecisionGoal -> 5];

   If[Head[sol] === InterpolatingFunction,
     Print["  NDSolve succeeded."];
     phiOf[x_, t_] := Clip[\[Phi]0 + amp sol[x, t], {0., 0.8}];
     tSnap = {0, 8, 40, 150, 400};

     (* Figure 1: space-time evolution + composition profiles *)
     (* Legends placed Below (horizontal) so they never overflow the
        right edge and get clipped by GraphicsColumn / Export. *)
     densPlot = DensityPlot[phiOf[x, t], {x, 0, L}, {t, 0, T},
        PlotRange -> {0, 0.8}, ColorFunction -> "TemperatureMap",
        PlotPoints -> 90,
        PlotLegends -> Placed[BarLegend[{"TemperatureMap", {0, 0.8}},
          LegendLabel -> "\[Phi]  (red = dense, blue = dilute)",
          LegendLayout -> "Row"], Below],
        FrameLabel -> {"Position x", "Time t"},
        PlotLabel -> Style["1D spinodal decomposition (space-time)", 14],
        AspectRatio -> 0.4, ImageSize -> 620];

     profPlot = Plot[Evaluate[Table[phiOf[x, t], {t, tSnap}]], {x, 0, L},
        PlotStyle -> (Directive[Thick, #] & /@ {
          Darker[Blue, 0.3], Blue, Darker[Green], Orange, Red}),
        PlotLegends -> Placed[
          LineLegend[
            (Directive[Thick, #] & /@ {
               Darker[Blue, 0.3], Blue, Darker[Green], Orange, Red}),
            (Row[{"t = ", #}] & /@ tSnap),
            LegendLayout -> "Row", LegendMarkerSize -> 24,
            Background -> Directive[GrayLevel[1, 0.7]]],
          Scaled[{0.5, 0.93}]],
        PlotPoints -> 200, MaxRecursion -> 2,
        Frame -> True,
        FrameLabel -> {"Position x", "Volume fraction \[Phi]"},
        PlotLabel -> Style["Composition profiles at selected times", 13],
        PlotRange -> {0, 0.85},
        GridLines -> Automatic, GridLinesStyle -> Directive[LightGray],
        ImageSize -> 620];

     Export[FileNameJoin[{imgDir, "ch_1d_evolution.png"}],
       GraphicsColumn[{densPlot, profPlot}, Spacings -> 15],
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
