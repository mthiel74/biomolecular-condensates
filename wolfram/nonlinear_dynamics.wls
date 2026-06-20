#!/usr/bin/env wolframscript
(* ::Package:: *)

(*  Nonlinear-dynamics analysis of the Cahn-Hilliard + reaction model.

    Pure Wolfram Language companion to src/nonlinear_analysis.py.  WL is
    used where it is strongest: the analytical linear-stability, bifurcation
    and reaction-coupling results, plus a 1D double-well Cahn-Hilliard solved
    natively with NDSolve to exhibit the monotone decay of the Lyapunov
    functional.

    Output (docs/images/):
      nld_wl_dispersion.png      dispersion relation omega(k)
      nld_wl_mode_selection.png  fastest-growing wavelength over (phi0, kappa)
      nld_wl_phase_diagram.png   spinodal + binodal (common tangent) + pitchfork
      nld_wl_reaction.png        mass-conserving rate enhancement k_eff/k_out
      nld_wl_energy.png          1D double-well CH: F[phi] decreases monotonically
*)

repoRoot = ParentDirectory @ DirectoryName[$InputFileName];
imgDir   = FileNameJoin[{repoRoot, "docs", "images"}];
If[!DirectoryQ[imgDir], CreateDirectory[imgDir, CreateIntermediateDirectories -> True]];

Print["Nonlinear-dynamics figures (Wolfram Language)..."];

(* ---- shared Flory-Huggins thermodynamics ---- *)
nN1 = 100; nN2 = 1;
phiC = Sqrt[nN2/nN1]/(1 + Sqrt[nN2/nN1]);
chiC = 0.5 (1/Sqrt[nN1] + 1/Sqrt[nN2])^2;

fpp[phi_, chi_] := 1/(nN1 phi) + 1/(nN2 (1 - phi)) - 2 chi;          (* d2f/dphi2 *)
omega[k_, phi_, chi_, kappa_, M_] := -M k^2 (fpp[phi, chi] + kappa k^2);
kStar[phi_, chi_, kappa_] := Sqrt[-fpp[phi, chi]/(2 kappa)];

accent  = RGBColor[0.78, 0.12, 0.36];
accent2 = RGBColor[0.12, 0.44, 0.70];
ink     = GrayLevel[0.1];

(* ==================================================================
   1.  Dispersion relation  omega(k)  (Flory-Huggins, chi = 1.5)
   ================================================================== *)
Module[{chi = 1.5, kappa = 0.5, M = 1.0, phis, curves, dots},
   phis = {0.06, 0.10, 0.20, 0.30, 0.45, 0.60};
   curves = Table[
      omega[k, p, chi, kappa, M], {p, phis}];
   dots = Table[
      If[fpp[p, chi] < 0,
        {kStar[p, chi, kappa], omega[kStar[p, chi, kappa], p, chi, kappa, M]},
        Nothing], {p, phis}];
   Export[FileNameJoin[{imgDir, "nld_wl_dispersion.png"}],
     Show[
       Plot[Evaluate@curves, {k, 0, 3},
         PlotStyle -> (Directive[Thick, ColorData["SunsetColors"][#]] & /@
            Rescale[Range@Length@phis]),
         PlotLegends -> Placed[
            Map[Row[{"\!\(\*SubscriptBox[\(\[Phi]\), \(0\)]\) = ", #}] &, phis],
            {0.18, 0.30}],
         Frame -> True, GridLines -> Automatic,
         GridLinesStyle -> Directive[LightGray, Dotted],
         FrameLabel -> {"wavenumber  k", "growth rate  \[Omega](k)"},
         PlotLabel -> Style[
            "Cahn\[Dash]Hilliard dispersion  \[Omega](k) = \[Minus]M \!\(\*SuperscriptBox[\(k\), \(2\)]\)(f''(\!\(\*SubscriptBox[\(\[Phi]\), \(0\)]\)) + \[Kappa] \!\(\*SuperscriptBox[\(k\), \(2\)]\))",
            14],
         PlotRange -> {{0, 3}, {-1.6, 1.6}}, ImageSize -> 620],
       Graphics[{ink, PointSize[0.012], Point /@ dots,
          Text[Style["fastest mode k*", 10, ink], {1.7, 1.25}]}],
       Graphics[{Gray, Line[{{0, 0}, {3, 0}}]}]],
     ImageResolution -> 200];
   Print["  saved nld_wl_dispersion.png"];
];

(* ==================================================================
   2.  Mode selection: lambda* = 2 pi / k* over (phi0, kappa)
   ================================================================== *)
Module[{chi = 1.5},
   Export[FileNameJoin[{imgDir, "nld_wl_mode_selection.png"}],
     Show[
       DensityPlot[
         If[fpp[phi, chi] < 0, Min[2 Pi/kStar[phi, chi, kappa], 25.], Indeterminate],
         {phi, 0.005, 0.72}, {kappa, 0.05, 2.0},
         ColorFunction -> "SunsetColors", PlotPoints -> 90,
         PlotRange -> All, PlotLegends -> Automatic,
         FrameLabel -> {"mean composition  \!\(\*SubscriptBox[\(\[Phi]\), \(0\)]\)",
                        "gradient coefficient  \[Kappa]"},
         PlotLabel -> Style[
            "selected wavelength  \[Lambda]* = 2\[Pi]/k*  (white = stable, f'' > 0)", 14],
         ImageSize -> 620],
       ContourPlot[fpp[phi, chi] == 0, {phi, 0.005, 0.72}, {kappa, 0.05, 2.0},
         ContourStyle -> Directive[White, Thick]]],
     ImageResolution -> 200];
   Print["  saved nld_wl_mode_selection.png"];
];

(* ==================================================================
   3.  Phase diagram: spinodal (analytic) + binodal (common tangent)
       and the symmetric double-well pitchfork.
   ================================================================== *)
Module[{spin, gg, hh, mu, Pizero, dense, dilute, binL, binR, chis,
        ddwA, ddwB, panelFH, panelDW},

   (* spinodal: chi = 1/2 (1/(N1 phi) + 1/(N2 (1-phi))) *)
   spin[phi_] := 0.5 (1/(nN1 phi) + 1/(nN2 (1 - phi)));

   (* common tangent via osmotic pressure Pi=0 (dense) + mu match (dilute) *)
   gg[p_] := (Log[p] + 1)/nN1 - (Log[1 - p] + 1)/nN2;
   hh[p_] := p Log[p]/nN1 + (1 - p) Log[1 - p]/nN2;
   mu[p_, chi_] := gg[p] + chi (1 - 2 p);
   Pizero[p_, chi_] := hh[p] - p gg[p] + chi p^2;

   (* Parametrise the binodal by the dense composition phiR (monotonic).
      Pi(p,chi)=0 is LINEAR in chi  ->  chiOf[pR] = (pR g(pR) - h(pR))/pR^2.
      Dilute branch phiL from mu-matching in log phi (hugs phi=0). *)
   chiOf[pR_] := (pR gg[pR] - hh[pR])/pR^2;
   With[{tris = Table[Module[{chi = chiOf[pR], muc, u, pL},
        muc = mu[pR, chi];
        u = u /. Quiet@FindRoot[mu[Exp[u], chi] == muc,
              {u, -8, -700, Log[phiC] - 0.01}];
        pL = Exp[u];
        If[NumberQ[pL] && chiC < chi <= chiC*3 && 0 < pL < phiC,
          {chi, pL, pR}, Nothing]],
        {pR, Subdivide[phiC + 0.003, 0.93, 120]}]},
      binL = SortBy[{#[[2]], #[[1]]} & /@ tris, Last];   (* {phiL, chi} *)
      binR = SortBy[{#[[3]], #[[1]]} & /@ tris, Last]];   (* {phiR, chi} *)

   panelFH = Show[
      Plot[spin[phi], {phi, 0.004, 0.75},
        PlotStyle -> Directive[accent, Thick, Dashed],
        Frame -> True, GridLines -> Automatic,
        GridLinesStyle -> Directive[LightGray, Dotted],
        FrameLabel -> {"volume fraction  \[Phi]", "\[Chi]  (\[Proportional] 1/T)"},
        PlotRange -> {{0, 0.8}, {0, chiC*3}}, ImageSize -> 430],
      ListLinePlot[{binL, binR}, PlotStyle -> Directive[accent2, Thick]],
      Graphics[{
         ink, PointSize[0.02], Point[{phiC, chiC}],
         Text[Style["critical point", 10], {phiC + 0.20, chiC + 0.12}],
         Directive[RGBColor[0.85, 0.6, 0.1]], PointSize[0.028],
         Point[{0.30, 1.5}],
         Text[Style["solver operating point", 9], {0.45, 1.62}],
         Text[Style["two-phase", 11, accent2], {0.45, 1.1}],
         Text[Style["one phase", 11, Darker@Green], {0.42, 0.4}]}],
      PlotLabel -> Style["(a) Flory\[Dash]Huggins phase diagram", 13]];

   (* double-well pitchfork: binodal +-sqrt(a/b), spinodal +-sqrt(a/3b) *)
   ddwA = Plot[{Sqrt[a], -Sqrt[a], Sqrt[a/3], -Sqrt[a/3]}, {a, 0, 1.5},
      PlotStyle -> {Directive[accent2, Thick], Directive[accent2, Thick],
                    Directive[accent, Thick, Dashed], Directive[accent, Thick, Dashed]},
      Frame -> True, GridLines -> Automatic,
      GridLinesStyle -> Directive[LightGray, Dotted],
      FrameLabel -> {"control parameter  a  (\[Proportional] \[Chi] \[Minus] \!\(\*SubscriptBox[\(\[Chi]\), \(c\)]\))",
                     "order parameter  \[Phi] \[Minus] \!\(\*SubscriptBox[\(\[Phi]\), \(c\)]\)"},
      PlotRange -> {{0, 1.5}, {-1.35, 1.35}}, ImageSize -> 430,
      PlotLabel -> Style["(b) double-well pitchfork", 13],
      Filling -> {3 -> {4}}, FillingStyle -> Directive[Opacity[0.12, accent]]];
   ddwB = Graphics[{ink, PointSize[0.02], Point[{0, 0}],
      Text[Style["binodal", 10, accent2], {1.15, 1.18}],
      Text[Style["spinodal", 10, accent], {1.2, 0.45}]}];
   panelDW = Show[ddwA, ddwB];

   Export[FileNameJoin[{imgDir, "nld_wl_phase_diagram.png"}],
     GraphicsRow[{panelFH, panelDW}, Spacings -> 20, ImageSize -> 920],
     ImageResolution -> 200];
   Print["  saved nld_wl_phase_diagram.png"];
];

(* ==================================================================
   4.  Reaction coupling: mass-conserving rate enhancement
       E(KA,KB,p) = (p KA KB + (1-p)) / ((1+p(KA-1))(1+p(KB-1)))
   ================================================================== *)
Module[{E, p0 = 0.05, heat, curves, Kvals},
   E[KA_, KB_, p_] := (p KA KB + (1 - p))/((1 + p (KA - 1)) (1 + p (KB - 1)));

   heat = DensityPlot[Log10[E[10^lKA, 10^lKB, p0]], {lKA, 0, 3}, {lKB, 0, 3},
      ColorFunction -> "SunsetColors", PlotPoints -> 70, PlotRange -> All,
      FrameLabel -> {"\!\(\*SubscriptBox[\(log\), \(10\)]\) \!\(\*SubscriptBox[\(K\), \(A\)]\)",
                     "\!\(\*SubscriptBox[\(log\), \(10\)]\) \!\(\*SubscriptBox[\(K\), \(B\)]\)"},
      PlotLegends -> Placed[BarLegend[Automatic,
         LegendLabel -> "\!\(\*SubscriptBox[\(log\), \(10\)]\)(\!\(\*SubscriptBox[\(k\), \(eff\)]\)/\!\(\*SubscriptBox[\(k\), \(out\)]\))"], Right],
      PlotLabel -> Style["(a) enhancement over (\!\(\*SubscriptBox[\(K\), \(A\)]\), \!\(\*SubscriptBox[\(K\), \(B\)]\)),  p = 0.05", 13],
      ImageSize -> 430];

   Kvals = {10, 30, 100, 300, 1000};
   curves = LogPlot[Evaluate@Table[E[Kk, Kk, p], {Kk, Kvals}], {p, 0.0005, 1},
      PlotStyle -> (Directive[Thick, ColorData["TemperatureMap"][#]] & /@
         Rescale[Range@Length@Kvals]),
      PlotLegends -> Placed[
         Map[Row[{"K = ", #}] &, Kvals], {0.82, 0.72}],
      Frame -> True, GridLines -> Automatic,
      GridLinesStyle -> Directive[LightGray, Dotted],
      FrameLabel -> {"dense-phase area fraction  p",
                     "\!\(\*SubscriptBox[\(k\), \(eff\)]\)/\!\(\*SubscriptBox[\(k\), \(out\)]\)"},
      PlotLabel -> Style["(b) optimum at p* \[TildeTilde] 1/K  (gain \[TildeTilde] K/4, not \!\(\*SuperscriptBox[\(K\), \(2\)]\))", 13],
      PlotRange -> {{0, 1}, {1, 300}}, ImageSize -> 430];

   Export[FileNameJoin[{imgDir, "nld_wl_reaction.png"}],
     GraphicsRow[{heat, curves}, Spacings -> 20, ImageSize -> 920],
     ImageResolution -> 200];
   Print["  saved nld_wl_reaction.png"];
];

(* ==================================================================
   5.  1D double-well Cahn-Hilliard via NDSolve:
       Lyapunov functional F[phi] decays monotonically.
   ================================================================== *)
Module[{a = 1.0, b = 1.0, kappa = 1.0, M = 1.0, L = 80.0, T = 400.0,
        mu, sol, coeffs, icf, nPts = 256, Ftot, ts, Fs, spacetime},
   mu[ph_] := -a ph + b ph^3;
   SeedRandom[7];
   coeffs = RandomReal[{-1, 1}, 24];
   (* smooth, exactly periodic initial condition (sines vanish at x=0,L) *)
   icf[x_] := 0.02 Sum[coeffs[[n]] Sin[2 Pi n x/L], {n, 1, 24}];

   Print["  Solving 1D double-well Cahn-Hilliard (NDSolve)..."];
   sol = Quiet @ NDSolveValue[{
        D[ph[x, t], t] == M D[mu[ph[x, t]] - kappa D[ph[x, t], x, x], x, x],
        ph[x, 0] == icf[x],
        ph[0, t] == ph[L, t],
        Derivative[1, 0][ph][0, t] == Derivative[1, 0][ph][L, t]},
      ph, {x, 0, L}, {t, 0, T},
      Method -> {"MethodOfLines",
         "SpatialDiscretization" -> {"TensorProductGrid",
            "MinPoints" -> nPts, "MaxPoints" -> nPts}},
      MaxSteps -> 500000, AccuracyGoal -> 5, PrecisionGoal -> 5];

   If[Head[sol] === InterpolatingFunction,
      Print["  NDSolve OK."];
      Ftot[t_?NumericQ] := NIntegrate[
         -a/2 sol[x, t]^2 + b/4 sol[x, t]^4
            + kappa/2 Derivative[1, 0][sol][x, t]^2,
         {x, 0, L}, AccuracyGoal -> 4];
      ts = Subdivide[0.5, T, 60];
      Fs = Ftot /@ ts;
      spacetime = DensityPlot[sol[x, t], {x, 0, L}, {t, 0, T},
         ColorFunction -> "BalancedHue", PlotPoints -> 90, PlotRange -> {-1, 1},
         FrameLabel -> {"position  x", "time  t"},
         PlotLabel -> Style["(a) 1D spinodal decomposition + coarsening", 13],
         PlotLegends -> Automatic, AspectRatio -> 0.7, ImageSize -> 430];
      Export[FileNameJoin[{imgDir, "nld_wl_energy.png"}],
         GraphicsRow[{
            spacetime,
            ListLinePlot[Transpose[{ts, Fs}],
               PlotStyle -> Directive[accent, Thick], Filling -> Bottom,
               FillingStyle -> Directive[Opacity[0.12, accent]],
               Frame -> True, GridLines -> Automatic,
               GridLinesStyle -> Directive[LightGray, Dotted],
               FrameLabel -> {"time  t", "free energy  F[\[Phi]]"},
               PlotLabel -> Style["(b) Lyapunov functional F[\[Phi]] decreases", 13],
               ImageSize -> 430]},
            Spacings -> 20, ImageSize -> 920],
         ImageResolution -> 200];
      Print["  saved nld_wl_energy.png"],
      Print["  WARN: NDSolve failed; skipping nld_wl_energy.png"]];
];

Print["Nonlinear-dynamics figures complete."];
