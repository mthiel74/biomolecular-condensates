#!/usr/bin/env wolframscript
(* ::Package:: *)

(*  Phase-3 extensions of the condensate model.

    Pure Wolfram Language companion to the Phase-3 Python modules
    (src/multicomponent.py, src/darwinian_selection.py,
    src/reaction_diffusion_cargo.py, src/wet_dry_cycling.py,
    src/parameter_space.py).  WL is used where it is strongest: the symbolic
    critical point of the N-component mixture, the chemostat selection ODE
    solved natively with NDSolve, the Turing-threshold linear stability, a 1D
    wet-dry Cahn-Hilliard solved with NDSolve, and the analytic regime planes.

    Output (docs/images/):
      p3_wl_multicomponent.png   spinodal overlay + critical-point migration
      p3_wl_selection.png        chemostat competition: mean fitness rises
      p3_wl_turing.png           Schnakenberg Turing growth rate: inside vs out
      p3_wl_wetdry.png           1D CH with cyclic a(t): dissolve / reform
      p3_wl_paramspace.png       (chi,N) two-phase + (K,p) rate enhancement
*)

repoRoot = ParentDirectory @ DirectoryName[$InputFileName];
imgDir   = FileNameJoin[{repoRoot, "docs", "images"}];
If[!DirectoryQ[imgDir], CreateDirectory[imgDir, CreateIntermediateDirectories -> True]];

Print["Phase-3 extension figures (Wolfram Language)..."];

accent  = RGBColor[0.78, 0.12, 0.36];
accent2 = RGBColor[0.12, 0.44, 0.70];
gold    = RGBColor[0.91, 0.64, 0.24];
ink     = GrayLevel[0.1];
seq     = "SunsetColors";

(* ==================================================================
   1.  Multi-component condensates (Qian & Knowles 2022).
       Effective co-condensation free energy of N equal components gives
       spinodal  chiS(phi) = 2 (1/(N phi) + 1/(1-phi))   [QK convention],
       with critical point  phiC = 1/(1+Sqrt[N]),  chiC = 2 (1+Sqrt[N])^2/N.
   ================================================================== *)
Module[{Ns = {2, 3, 5, 10}, chiS, phiC, chiC, cols, curves, pts, migr},
   chiS[phi_, n_] := 2 (1/(n phi) + 1/(1 - phi));
   phiC[n_] := 1/(1 + Sqrt[n]);
   chiC[n_] := 2 (1 + Sqrt[n])^2/n;
   (* symbolic verification: critical point minimises the spinodal *)
   Do[
     If[Simplify[(D[chiS[phi, n], phi] /. phi -> phiC[n])] =!= 0,
        Print["  [warn] critical-point check failed for N=", n]],
     {n, Ns}];
   cols = Table[ColorData[seq][0.12 + 0.62 (i - 1)/(Length[Ns] - 1)],
                {i, Length[Ns]}];
   curves = Table[chiS[phi, Ns[[i]]], {i, Length[Ns]}];
   pts = Table[{phiC[n], chiC[n]}, {n, Ns}];
   migr = Plot[Evaluate@curves, {phi, 0.01, 0.99},
      PlotStyle -> (Directive[#, Thick] & /@ cols),
      Frame -> True, GridLines -> Automatic,
      GridLinesStyle -> Directive[LightGray, Dotted],
      FrameLabel -> {"total condensed fraction  \[Phi]",
                     "interaction  \[Chi]  (Qian\[Dash]Knowles)"},
      PlotRange -> {{0, 1}, {0, 18}},
      PlotLegends -> Placed[
        (Style["N=" <> ToString[#], 11] & /@ Ns), {0.82, 0.7}],
      PlotLabel -> Style[
        "Multi-component spinodals: critical point migrates as 1/(1+\!\(\*SqrtBox[\(N\)]\))",
        13], ImageSize -> 560];
   Export[FileNameJoin[{imgDir, "p3_wl_multicomponent.png"}],
     Show[migr, Graphics[{Black, PointSize[0.022], Point[pts],
        White, PointSize[0.012], Point[pts]}]],
     ImageResolution -> 200];
   Print["  saved p3_wl_multicomponent.png"];
];

(* ==================================================================
   2.  Darwinian selection -- chemostat competition ODE (NDSolve).
       K genotypes g_i compete for a shared nutrient with dilution:
         c_i' = r0 g_i c_i Nutr - death c_i
         Nutr' = supply (1-Nutr) - use Nutr Sum[r0 g_i c_i]
       The fittest genotype takes over; the cargo-weighted mean fitness rises.
   ================================================================== *)
Module[{r0 = 0.6, death = 0.16, supply = 0.22, use = 0.05, gs, nC, vars,
        eqs, ic, sol, tmax = 120, meanFit, ctot, fitPlot, cPlot},
   gs = {0.6, 0.85, 1.1, 1.35, 1.6, 1.8};
   nC = Length[gs];
   vars = Table[Subscript[c, i][t], {i, nC}];
   eqs = Join[
     Table[Subscript[c, i]'[t] ==
        r0 gs[[i]] Subscript[c, i][t] nut[t] - death Subscript[c, i][t],
        {i, nC}],
     {nut'[t] == supply (1 - nut[t]) -
        use nut[t] Sum[r0 gs[[j]] Subscript[c, j][t], {j, nC}]}];
   ic = Join[Table[Subscript[c, i][0] == 0.1, {i, nC}], {nut[0] == 1.0}];
   sol = NDSolve[Join[eqs, ic], Append[Table[Subscript[c, i], {i, nC}], nut],
                 {t, 0, tmax}][[1]];
   meanFit[tt_] := (Sum[gs[[i]] (Subscript[c, i][tt] /. sol), {i, nC}])/
                   (Sum[(Subscript[c, i][tt] /. sol), {i, nC}]);
   fitPlot = Plot[meanFit[t], {t, 0, tmax},
      PlotStyle -> Directive[accent, Thickness[0.01]],
      Frame -> True, GridLines -> Automatic,
      GridLinesStyle -> Directive[LightGray, Dotted],
      FrameLabel -> {"time", "mean fitness  \[LeftAngleBracket]g\[RightAngleBracket]"},
      PlotRange -> {{0, tmax}, {0.55, 1.85}},
      Epilog -> {gold, Dashed, Line[{{0, Max[gs]}, {tmax, Max[gs]}}],
                 Text[Style["g_max", 10, gold], {tmax 0.9, Max[gs] - 0.05}]},
      PlotLabel -> Style["Chemostat selection: mean fitness rises", 13],
      ImageSize -> 430];
   cPlot = Plot[Evaluate[Table[Subscript[c, i][t] /. sol, {i, nC}]],
      {t, 0, tmax},
      PlotStyle -> Table[Directive[
         ColorData[seq][0.1 + 0.8 (i - 1)/(nC - 1)], Thick], {i, nC}],
      Frame -> True, GridLines -> Automatic,
      GridLinesStyle -> Directive[LightGray, Dotted],
      FrameLabel -> {"time", "genotype abundance  \!\(\*SubscriptBox[\(c\), \(i\)]\)"},
      PlotLegends -> Placed[
        (Style["g=" <> ToString[#], 9] & /@ gs), Right],
      PlotLabel -> Style["The fittest replicator wins", 13],
      ImageSize -> 470];
   Export[FileNameJoin[{imgDir, "p3_wl_selection.png"}],
     GraphicsRow[{fitPlot, cPlot}, Spacings -> 20, ImageSize -> 940],
     ImageResolution -> 200];
   Print["  saved p3_wl_selection.png"];
];

(* ==================================================================
   3.  Intra-droplet Turing patterns -- Schnakenberg linear stability.
       f = a - u + u^2 v,  g = b - u^2 v.  Growth rate of mode k:
         lambda(k) = larger eigenvalue of  J - k^2 diag(Du,Dv).
       Inside the droplet Dv/Du is above the Turing threshold (a band of k is
       unstable); outside it is below (all modes stable).
   ================================================================== *)
Module[{a = 0.1, b = 1.05, us, vs, fu, fv, gu, gv, growth, DuIn = 0.30,
        DvIn = 9.0, DuOut = 2.4, DvOut = 14.0, gIn = 3.4, gOut = 1.4, kmax},
   us = a + b; vs = b/us^2;
   fu = -1 + 2 us vs; fv = us^2; gu = -2 us vs; gv = -us^2;
   growth[k_, Du_, Dv_, gam_] := Module[{J, tr, det},
     J = gam {{fu, fv}, {gu, gv}} - k^2 {{Du, 0}, {0, Dv}};
     tr = Tr[J]; det = Det[J];
     (tr + Sqrt[tr^2 - 4 det])/2];
   kmax = 6;
   Export[FileNameJoin[{imgDir, "p3_wl_turing.png"}],
     Plot[{Re@growth[k, DuIn, DvIn, gIn], Re@growth[k, DuOut, DvOut, gOut]},
       {k, 0, kmax},
       PlotStyle -> {Directive[accent, Thickness[0.01]],
                     Directive[accent2, Thickness[0.008], Dashed]},
       Frame -> True, GridLines -> Automatic,
       GridLinesStyle -> Directive[LightGray, Dotted],
       FrameLabel -> {"wavenumber  k", "growth rate  \[Lambda](k)"},
       PlotRange -> {{0, kmax}, {-2.5, 1.2}},
       PlotLegends -> Placed[{"inside droplet (dense): unstable band",
          "dilute phase: all modes stable"}, {0.62, 0.25}],
       Epilog -> {Gray, Line[{{0, 0}, {kmax, 0}}]},
       PlotLabel -> Style[
         "Schnakenberg Turing instability: patterns only inside the condensate",
         12.5], ImageSize -> 620],
     ImageResolution -> 200];
   Print["  saved p3_wl_turing.png"];
];

(* ==================================================================
   4.  Wet-dry cycling -- 1D Cahn-Hilliard with time-dependent depth a(t).
       phi_t = M d^2/dx^2 [ -a(t) phi + b phi^3 - kappa phi_xx ],  periodic.
       Cyclic a(t) dissolves (a<0) and re-forms (a>0) the domains; the
       space-time plot shows the population repeatedly reset.
   ================================================================== *)
Module[{L = 60.0, M = 1.0, b = 1.0, kappa = 1.0, a0 = 0.6, dA = 1.4,
        period = 60.0, tmax = 180.0, aF, sol, dens, nx = 160, seed},
   aF[t_] := a0 + dA Sin[2 Pi t/period];
   SeedRandom[7];
   Quiet@Check[
     sol = NDSolveValue[{
        D[phi[t, x], t] == M D[
           -aF[t] phi[t, x] + b phi[t, x]^3 - kappa D[phi[t, x], {x, 2}],
           {x, 2}],
        phi[0, x] == 0.05 RandomReal[{-1, 1}] +
           0.04 Sum[Sin[2 Pi n x/L + n], {n, 3, 11}],
        phi[t, 0] == phi[t, L], Derivative[0, 1][phi][t, 0] == Derivative[0, 1][phi][t, L],
        Derivative[0, 2][phi][t, 0] == Derivative[0, 2][phi][t, L],
        Derivative[0, 3][phi][t, 0] == Derivative[0, 3][phi][t, L]},
        phi, {t, 0, tmax}, {x, 0, L},
        Method -> {"MethodOfLines",
           "SpatialDiscretization" -> {"TensorProductGrid",
              "MinPoints" -> nx, "MaxPoints" -> nx, "DifferenceOrder" -> 4}}],
     $Failed];
   If[sol === $Failed,
     Print["  [skip p3_wl_wetdry.png: NDSolve failed]"],
     dens = DensityPlot[sol[t, x], {t, 0, tmax}, {x, 0, L},
        PlotPoints -> 90, ColorFunction -> "RedBlueTones", AspectRatio -> 0.6,
        Frame -> True, FrameLabel -> {"time  t", "position  x"},
        PlotLabel -> Style[
          "1D wet\[Dash]dry Cahn\[Dash]Hilliard: domains dissolve and re-form each cycle",
          12.5], ImageSize -> 680, PlotLegends -> BarLegend[Automatic]];
     Export[FileNameJoin[{imgDir, "p3_wl_wetdry.png"}], dens,
        ImageResolution -> 200];
     Print["  saved p3_wl_wetdry.png"]];
];

(* ==================================================================
   5.  Parameter-space planes (analytic).
       (chi, N): two-phase iff chi > chiC(N) = (1/2)(1 + 1/Sqrt[N])^2.
       (K, p):  rate enhancement  R = (p K^2 + (1-p))/(1 + p(K-1))^2,
                optimal ridge p* ~ 1/K.
   ================================================================== *)
Module[{chiC, twoPhase, ridge, ratePlot, K},
   chiC[n_] := (1/2)(1 + 1/Sqrt[n])^2;
   twoPhase = RegionPlot[chi > chiC[n], {chi, 0, 3.5}, {n, 1, 60},
      PlotPoints -> 60, BoundaryStyle -> Directive[Cyan, Thick],
      PlotStyle -> accent, Frame -> True,
      FrameLabel -> {"interaction  \[Chi]", "chain length  N"},
      PlotLabel -> Style["(\[Chi],N): one-phase / two-phase", 12.5],
      Epilog -> {Text[Style["two-phase", 12, White], {2.6, 45}],
                 Text[Style["one-phase", 12, ink], {0.6, 45}]},
      ImageSize -> 430];
   ratePlot = DensityPlot[
      (p K^2 + (1 - p))/(1 + p (K - 1))^2, {K, 1, 300}, {p, 0.001, 0.999},
      ScalingFunctions -> {"Log", None}, PlotPoints -> 70,
      ColorFunction -> "SunsetColors", Frame -> True,
      FrameLabel -> {"partition coefficient  K", "dense fraction  p"},
      PlotLabel -> Style["(K,p): rate enhancement \!\(\*SubscriptBox[\(k\), \(eff\)]\)/\!\(\*SubscriptBox[\(k\), \(out\)]\)", 12.5],
      Epilog -> {White, Dashed, Thick,
         Line[Table[{K, 1/K}, {K, 1, 300, 1}]]},
      PlotLegends -> BarLegend[Automatic], ImageSize -> 470];
   Export[FileNameJoin[{imgDir, "p3_wl_paramspace.png"}],
     GraphicsRow[{twoPhase, ratePlot}, Spacings -> 24, ImageSize -> 940],
     ImageResolution -> 200];
   Print["  saved p3_wl_paramspace.png"];
];

Print["Done (Phase-3 Wolfram figures)."];
