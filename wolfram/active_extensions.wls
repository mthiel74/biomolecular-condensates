#!/usr/bin/env wolframscript
(* ::Package:: *)

(*  Active / chemically-driven extensions of the condensate model.

    Pure Wolfram Language companion to the Phase-2 Python modules
    (src/active_model_b.py, src/fuel_driven_droplets.py,
    src/thermal_fission.py, src/stickers_spacers.py).  WL is used where it is
    strongest: the analytic linear-stability of the fuel-driven model (finite
    unstable band, dissolution threshold, k^-1/2 arrest length), the
    sticker-modified phase diagram, and a 1D reaction-Cahn-Hilliard solved
    natively with NDSolve that exhibits arrested coarsening.

    Output (docs/images/):
      act_wl_fuel_dispersion.png   omega(k) family: finite band + dissolution
      act_wl_arrest_length.png     arrest length lambda_max(k) ~ k^-1/2
      act_wl_stickers.png          sticker valence widens the binodal
      act_wl_reaction_1d.png       1D CH: passive coarsens vs reaction arrests
*)

repoRoot = ParentDirectory @ DirectoryName[$InputFileName];
imgDir   = FileNameJoin[{repoRoot, "docs", "images"}];
If[!DirectoryQ[imgDir], CreateDirectory[imgDir, CreateIntermediateDirectories -> True]];

Print["Active-extension figures (Wolfram Language)..."];

accent  = RGBColor[0.78, 0.12, 0.36];
accent2 = RGBColor[0.12, 0.44, 0.70];
gold    = RGBColor[0.91, 0.64, 0.24];
ink     = GrayLevel[0.1];

(* ==================================================================
   1.  Fuel-driven dispersion (double well, f'' = -a):
       omega(k) = -M k^2 (f'' + kappa k^2) - kReact
                = M k^2 (a - kappa k^2) - kReact
       Finite unstable band; dissolution when kReact > M a^2/(4 kappa).
   ================================================================== *)
Module[{a = 1.0, kappa = 1.0, M = 1.0, krs, kdiss, curves},
   kdiss = M a^2/(4 kappa);                       (* = 0.25 *)
   krs = {0.0, 0.05, 0.15, kdiss};
   curves = Table[M k^2 (a - kappa k^2) - kr, {kr, krs}];
   Export[FileNameJoin[{imgDir, "act_wl_fuel_dispersion.png"}],
     Show[
       Plot[Evaluate@curves, {k, 0, 1.3},
         PlotStyle -> {Directive[ink, Thick],
                       Directive[accent2, Thick],
                       Directive[accent, Thick],
                       Directive[gold, Thick, Dashed]},
         Frame -> True, GridLines -> Automatic,
         GridLinesStyle -> Directive[LightGray, Dotted],
         FrameLabel -> {"wavenumber  k", "growth rate  \[Omega](k)"},
         PlotLegends -> Placed[{"passive  k=0", "k=0.05", "k=0.15",
            "\!\(\*SubscriptBox[\(k\), \(diss\)]\)=0.25"}, {0.22, 0.28}],
         PlotLabel -> Style[
            "Fuel-driven dispersion  \[Omega] = M \!\(\*SuperscriptBox[\(k\), \(2\)]\)(a \[Minus] \[Kappa] \!\(\*SuperscriptBox[\(k\), \(2\)]\)) \[Minus] \!\(\*SubscriptBox[\(k\), \(react\)]\):  reaction stabilises low k",
            13],
         PlotRange -> {{0, 1.3}, {-0.3, 0.3}}, ImageSize -> 640],
       Graphics[{Gray, Line[{{0, 0}, {1.3, 0}}]}]],
     ImageResolution -> 200];
   Print["  saved act_wl_fuel_dispersion.png"];
];

(* ==================================================================
   2.  Arrest length.  Band edges solve M kappa k^4 - M a k^2 + kReact = 0:
         kmin^2 = (M a - Sqrt[(M a)^2 - 4 M kappa kReact])/(2 M kappa).
       lambda_max = 2 pi / kmin  ~  k^-1/2  (Zwicker R* ~ sqrt(D/k)),
       floored by the fastest passive mode lambda* = 2 pi / sqrt(a/2kappa).
   ================================================================== *)
Module[{a = 1.0, kappa = 1.0, M = 1.0, kdiss, kmin, lamMax, lamStar},
   kdiss = M a^2/(4 kappa);
   kmin[kr_] := Sqrt[(M a - Sqrt[(M a)^2 - 4 M kappa kr])/(2 M kappa)];
   lamMax[kr_] := 2 Pi/kmin[kr];
   lamStar = 2 Pi/Sqrt[a/(2 kappa)];
   Export[FileNameJoin[{imgDir, "act_wl_arrest_length.png"}],
     Show[
       LogLogPlot[{lamMax[kr], lamStar, lamMax[0.05] Sqrt[0.05/kr]},
         {kr, 0.004, kdiss 0.999},
         PlotStyle -> {Directive[accent2, Thick],
                       Directive[ink, Thick, Dashed],
                       Directive[gold, Thick, Dotted]},
         Frame -> True, GridLines -> Automatic,
         GridLinesStyle -> Directive[LightGray, Dotted],
         FrameLabel -> {"reaction rate  \!\(\*SubscriptBox[\(k\), \(react\)]\)",
                        "arrested length scale"},
         PlotLegends -> Placed[{
            "\[Lambda]\!\(\*SubscriptBox[\(\), \(max\)]\) = 2\[Pi]/\!\(\*SubscriptBox[\(k\), \(min\)]\)",
            "\[Lambda]* (fastest mode floor)",
            "\!\(\*SuperscriptBox[\(k\), \(-1/2\)]\) (Zwicker)"}, {0.7, 0.78}],
         PlotLabel -> Style[
            "Reaction selects a length: \[Lambda]\!\(\*SubscriptBox[\(\), \(max\)]\)~\!\(\*SuperscriptBox[\(k\), \(-1/2\)]\), merges with floor at \!\(\*SubscriptBox[\(k\), \(diss\)]\)",
            13],
         ImageSize -> 640]],
     ImageResolution -> 200];
   Print["  saved act_wl_arrest_length.png"];
];

(* ==================================================================
   3.  Stickers-and-spacers: association free energy widens the binodal.
       f(phi) = phi Log[phi]/N1 + (1-phi) Log[1-phi]/N2 + chi phi(1-phi)
                - rho_s Log[1 + K rho_s],   rho_s = (valence/N1) phi.
       Binodal by the lower convex hull (parameter-free double tangent).
   ================================================================== *)
Module[{nN1 = 100, nN2 = 1, K = 80.0, vals = {0, 2, 4, 8}, chiC, phiGrid,
        fOf, binodal, curves, cols},
   chiC = 0.5 (1/Sqrt[nN1] + 1/Sqrt[nN2])^2;
   phiGrid = Subdivide[0.0008, 0.9992, 1200];
   fOf[phi_, chi_, v_] := phi Log[phi]/nN1 + (1 - phi) Log[1 - phi]/nN2 +
      chi phi (1 - phi) -
      (v/nN1 phi) Log[1 + K (v/nN1) phi];

   (* tie-line endpoints from the lower convex hull of {phi, f} *)
   binodal[v_] := Module[{out = {}, chi, pts, hull, seg},
      Do[
        pts = Transpose[{phiGrid, fOf[#, chi, v] & /@ phiGrid}];
        (* lower hull: keep points below the chord of their neighbours *)
        hull = Module[{h = {}},
           Do[
             While[Length[h] >= 2 &&
                (h[[-1, 2]] - h[[-2, 2]]) (pts[[i, 1]] - h[[-2, 1]]) >=
                (pts[[i, 2]] - h[[-2, 2]]) (h[[-1, 1]] - h[[-2, 1]]),
                h = Most[h]];
             h = Append[h, pts[[i]]], {i, Length[pts]}];
           h];
        (* widest hull edge that skips interior points = the tie line *)
        seg = MaximalBy[
           Select[Partition[hull, 2, 1], (#[[2, 1]] - #[[1, 1]]) > 0.02 &],
           (#[[2, 1]] - #[[1, 1]]) &, 1];
        If[seg =!= {}, AppendTo[out, {seg[[1, 1, 1]], seg[[1, 2, 1]], chi}]],
        {chi, chiC*0.45, chiC*3, (chiC*3 - chiC*0.45)/120}];
      out];

   cols = {ink, RGBColor[0.55, 0.2, 0.55], accent, gold};
   curves = Table[
      With[{b = binodal[vals[[j]]]},
        If[b === {}, Nothing,
          {ListLinePlot[{#[[1]], #[[3]]} & /@ b,
              PlotStyle -> Directive[cols[[j]], Thick]],
           ListLinePlot[{#[[2]], #[[3]]} & /@ b,
              PlotStyle -> Directive[cols[[j]], Thick]]}]],
      {j, Length[vals]}];

   Export[FileNameJoin[{imgDir, "act_wl_stickers.png"}],
     Show[Flatten[curves],
       Frame -> True, GridLines -> Automatic,
       GridLinesStyle -> Directive[LightGray, Dotted],
       FrameLabel -> {"volume fraction  \[Phi]", "\[Chi]  (\[Proportional] 1/T)"},
       PlotLabel -> Style[
          "Stickers-and-spacers: valence (0,2,4,8) widens the two-phase region",
          13],
       PlotRange -> {{0, 1}, {0, chiC*3}}, ImageSize -> 640,
       Epilog -> {ink, PointSize[0.018], Point[{Sqrt[nN2/nN1]/(1 + Sqrt[nN2/nN1]), chiC}],
          Text[Style["valence \[UpArrow] \[RightArrow] demix at lower \[Chi]", 11, accent],
             {0.55, chiC*2.2}]}],
     ImageResolution -> 200];
   Print["  saved act_wl_stickers.png"];
];

(* ==================================================================
   4.  1D reaction-Cahn-Hilliard via NDSolve:
       d phi/dt = M d^2/dx^2 (mu) - kReact phi,  mu = -a phi + b phi^3 - kappa phi_xx.
       Passive (kReact=0) coarsens to one domain; reaction (kReact>0) arrests
       into a steady periodic emulsion.
   ================================================================== *)
Module[{a = 1.0, b = 1.0, kappa = 1.0, M = 1.0, L = 100.0, T = 600.0,
        nPts = 256, coeffs, icf, solve, sp0, spK},
   SeedRandom[5];
   coeffs = RandomReal[{-1, 1}, 30];
   icf[x_] := 0.05 Sum[coeffs[[n]] Sin[2 Pi n x/L], {n, 1, 30}];

   solve[kr_] := Quiet @ NDSolveValue[{
        D[ph[x, t], t] ==
           M D[(-a ph[x, t] + b ph[x, t]^3) - kappa D[ph[x, t], x, x], x, x]
           - kr ph[x, t],
        ph[x, 0] == icf[x],
        ph[0, t] == ph[L, t],
        Derivative[1, 0][ph][0, t] == Derivative[1, 0][ph][L, t]},
      ph, {x, 0, L}, {t, 0, T},
      Method -> {"MethodOfLines",
         "SpatialDiscretization" -> {"TensorProductGrid",
            "MinPoints" -> nPts, "MaxPoints" -> nPts}},
      MaxSteps -> 800000, AccuracyGoal -> 5, PrecisionGoal -> 5];

   Print["  Solving 1D reaction-Cahn-Hilliard (NDSolve)..."];
   sp0 = solve[0.0];
   spK = solve[0.05];
   If[Head[sp0] === InterpolatingFunction && Head[spK] === InterpolatingFunction,
      Export[FileNameJoin[{imgDir, "act_wl_reaction_1d.png"}],
        GraphicsRow[{
           DensityPlot[sp0[x, t], {x, 0, L}, {t, 0, T},
              ColorFunction -> "BalancedHue", PlotPoints -> 90,
              PlotRange -> {-1.2, 1.2}, FrameLabel -> {"x", "t"},
              PlotLabel -> Style["(a) passive: coarsens to one domain", 13],
              AspectRatio -> 0.8, ImageSize -> 430],
           DensityPlot[spK[x, t], {x, 0, L}, {t, 0, T},
              ColorFunction -> "BalancedHue", PlotPoints -> 90,
              PlotRange -> {-1.2, 1.2}, FrameLabel -> {"x", "t"},
              PlotLabel -> Style["(b) reaction k=0.05: arrested emulsion", 13],
              AspectRatio -> 0.8, ImageSize -> 430]},
           Spacings -> 20, ImageSize -> 920],
        ImageResolution -> 200];
      Print["  saved act_wl_reaction_1d.png"],
      Print["  WARN: NDSolve failed; skipping act_wl_reaction_1d.png"]];
];

Print["Active-extension figures complete."];
