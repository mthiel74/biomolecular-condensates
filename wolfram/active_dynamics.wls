#!/usr/bin/env wolframscript
(* ::Package:: *)

(*  Active / chemically-driven condensate dynamics -- pure Wolfram Language
    companion to src/active_models.py and src/fh_extensions.py (Phase 2).

    Output (docs/images/):
      nld_wl_active_dispersion.png   reaction-modified dispersion omega(k)
      nld_wl_dispersion_surface.png  3D dispersion surface omega(k, phi0)
      nld_wl_stickers.png            stickers-and-spacers phase diagram
      nld_wl_reaction_1d.png         1D reaction-CH: arrest of coarsening
*)

repoRoot = ParentDirectory @ DirectoryName[$InputFileName];
imgDir   = FileNameJoin[{repoRoot, "docs", "images"}];
If[!DirectoryQ[imgDir], CreateDirectory[imgDir, CreateIntermediateDirectories -> True]];

Print["Active-dynamics figures (Wolfram Language)..."];

accent  = RGBColor[0.78, 0.12, 0.36];
accent2 = RGBColor[0.12, 0.44, 0.70];
ink     = GrayLevel[0.1];

(* ==================================================================
   1.  Reaction-modified dispersion: omega(k) = -M k^2 (f'' + kappa k^2) - k_react
       (double well at phi0=0 => f'' = -a = -1)
   ================================================================== *)
Module[{fpp = -1.0, kappa = 1.0, M = 1.0, krs},
   krs = {0.0, 0.05, 0.15, 0.25};
   Export[FileNameJoin[{imgDir, "nld_wl_active_dispersion.png"}],
     Show[
       Plot[Evaluate@Table[-M k^2 (fpp + kappa k^2) - kr, {kr, krs}], {k, 0, 1.4},
         PlotStyle -> (Directive[Thick, ColorData["SunsetColors"][#]] & /@
            Rescale[Range@Length@krs]),
         PlotLegends -> Placed[
            Map[Row[{"\!\(\*SubscriptBox[\(k\), \(react\)]\) = ", #}] &, krs],
            {0.25, 0.28}],
         Frame -> True, GridLines -> Automatic,
         GridLinesStyle -> Directive[LightGray, Dotted],
         FrameLabel -> {"wavenumber  k", "growth rate  \[Omega](k)"},
         PlotLabel -> Style[
            "Reaction stabilises long-wavelength modes \[LongDash] arrested coarsening", 14],
         PlotRange -> {{0, 1.4}, {-0.35, 0.30}}, ImageSize -> 620],
       Graphics[{Gray, Line[{{0, 0}, {1.4, 0}}],
          Text[Style["finite unstable band \[Rule] selected length scale", 10, accent],
             {0.7, 0.22}]}]],
     ImageResolution -> 200];
   Print["  saved nld_wl_active_dispersion.png"];
];

(* ==================================================================
   2.  3D dispersion surface  omega(k, phi0)  (Flory-Huggins)
   ================================================================== *)
Module[{chi = 1.5, nN1 = 100, nN2 = 1, kappa = 0.5, M = 1.0, fpp},
   fpp[phi_] := 1/(nN1 phi) + 1/(nN2 (1 - phi)) - 2 chi;
   Export[FileNameJoin[{imgDir, "nld_wl_dispersion_surface.png"}],
     Plot3D[
       Clip[-M k^2 (fpp[phi] + kappa k^2), {-1.5, 1.5}],
       {phi, 0.02, 0.66}, {k, 0, 2.2},
       PlotPoints -> 60, MaxRecursion -> 2,
       ColorFunction -> "SunsetColors", PlotRange -> {-1.0, 1.4},
       MeshFunctions -> {#3 &}, Mesh -> {{0}},
       MeshStyle -> Directive[Cyan, Thick],
       AxesLabel -> {"\!\(\*SubscriptBox[\(\[Phi]\), \(0\)]\)", "k", "\[Omega]"},
       PlotLabel -> Style[
          "Dispersion surface  \[Omega](k, \!\(\*SubscriptBox[\(\[Phi]\), \(0\)]\))  (cyan ridge: neutral stability)", 14],
       ViewPoint -> {-2.4, -1.6, 1.4}, Lighting -> "Neutral",
       ImageSize -> 640, BoxRatios -> {1, 1, 0.7}],
     ImageResolution -> 200];
   Print["  saved nld_wl_dispersion_surface.png"];
];

(* ==================================================================
   3.  Stickers-and-spacers phase diagram (symmetric FH baseline)
       Semenov-Rubinstein association: p/(1-p)^2 = Ka phi_s,  phi_s = v phi
   ================================================================== *)
Module[{pbond, fbond, fbpp, chiSpin, cases, curves, pts, h = 10.^-4},
   pbond[phis_, Ka_] := With[{x = Ka Max[phis, 10.^-12]},
      ((2 x + 1) - Sqrt[4 x + 1])/(2 x)];
   fbond[phi_, v_, Ka_] := If[Ka <= 0, 0.,
      Module[{phis = v phi, p}, p = pbond[phis, Ka];
         phis (Log[Max[1 - p, 10.^-12]] + p/2)]];
   fbpp[phi_, v_, Ka_] :=
      (fbond[phi + h, v, Ka] - 2 fbond[phi, v, Ka] + fbond[phi - h, v, Ka])/h^2;
   chiSpin[phi_, v_, Ka_] := 0.5 (1/phi + 1/(1 - phi) + If[Ka > 0, fbpp[phi, v, Ka], 0.]);

   cases = {{1, 0, "no stickers (FH)"}, {2, 2, "v=2, \!\(\*SubscriptBox[\(K\), \(a\)]\)=2"},
            {3, 2, "v=3, \!\(\*SubscriptBox[\(K\), \(a\)]\)=2"},
            {4, 4, "v=4, \!\(\*SubscriptBox[\(K\), \(a\)]\)=4"}};
   curves = Table[
      With[{v = c[[1]], Ka = c[[2]]},
        Table[{phi, chiSpin[phi, v, Ka]}, {phi, 0.02, 0.98, 0.005}]],
      {c, cases}];

   Export[FileNameJoin[{imgDir, "nld_wl_stickers.png"}],
     ListLinePlot[curves,
       PlotStyle -> (Directive[Thick, ColorData["SunsetColors"][#]] & /@
          Rescale[Range@Length@cases]),
       PlotLegends -> Placed[cases[[All, 3]], {0.78, 0.78}],
       Frame -> True, GridLines -> Automatic,
       GridLinesStyle -> Directive[LightGray, Dotted],
       FrameLabel -> {"volume fraction  \[Phi]", "interaction  \[Chi]"},
       PlotLabel -> Style[
          "Stickers-and-spacers: multivalency expands the two-phase region", 14],
       PlotRange -> {{0, 1}, {0, 3.2}}, ImageSize -> 620],
     ImageResolution -> 200];
   Print["  saved nld_wl_stickers.png"];
];

(* ==================================================================
   4.  1D reaction-Cahn-Hilliard: arrest of coarsening (NDSolve)
       d phi/dt = M d^2/dx^2 [ f'(phi) - kappa phi_xx ] - kreact (phi - phi_ss)
       double well f'(phi) = -a phi + b phi^3
   ================================================================== *)
Module[{a = 1.0, b = 1.0, kappa = 1.0, M = 1.0, L = 100.0, T = 600.0,
        coeffs, icf, mu, solve, panels},
   mu[ph_] := -a ph + b ph^3;
   SeedRandom[11];
   coeffs = RandomReal[{-1, 1}, 30];
   icf[x_] := 0.03 Sum[coeffs[[n]] Sin[2 Pi n x/L], {n, 1, 30}];

   solve[kr_] := Quiet @ NDSolveValue[{
        D[ph[x, t], t] ==
           M D[mu[ph[x, t]] - kappa D[ph[x, t], x, x], x, x]
           - kr (ph[x, t] - 0),
        ph[x, 0] == icf[x],
        ph[0, t] == ph[L, t],
        Derivative[1, 0][ph][0, t] == Derivative[1, 0][ph][L, t]},
      ph, {x, 0, L}, {t, 0, T},
      Method -> {"MethodOfLines",
         "SpatialDiscretization" -> {"TensorProductGrid",
            "MinPoints" -> 200, "MaxPoints" -> 200}},
      MaxSteps -> 500000, AccuracyGoal -> 5, PrecisionGoal -> 5];

   Print["  Solving 1D reaction-CH (passive + active)..."];
   panels = Table[
      With[{sol = solve[kr]},
        If[Head[sol] === InterpolatingFunction,
          DensityPlot[sol[x, t], {x, 0, L}, {t, 0, T},
            ColorFunction -> "BalancedHue", PlotPoints -> 90, PlotRange -> {-1.2, 1.2},
            FrameLabel -> {"position  x", "time  t"},
            PlotLabel -> Style[
               If[kr == 0, "(a) passive: coarsens to one domain",
                  "(b) active (" <> ToString[kr] <> "): arrested pattern"], 13],
            AspectRatio -> 1, ImageSize -> 430],
          Graphics[Text["NDSolve failed"], ImageSize -> 430]]],
      {kr, {0.0, 0.03}}];

   Export[FileNameJoin[{imgDir, "nld_wl_reaction_1d.png"}],
     GraphicsRow[panels, Spacings -> 20, ImageSize -> 920],
     ImageResolution -> 200];
   Print["  saved nld_wl_reaction_1d.png"];
];

Print["Active-dynamics figures complete."];
