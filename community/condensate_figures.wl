(* ::Package:: *)

(*  condensate_figures.wl  ---------------------------------------------

    Model functions and figure builders for the "Biomolecular Condensates
    and the Origin of Life" Wolfram Community notebook.

    Loading this file defines functions only (no side effects, no Export,
    no Print). The notebook does

        Get["condensate_figures.wl"]

    once, then calls a figXxx[...] builder just before each figure so that
    every output is generated *inside the notebook*. The same builders are
    the single source of truth for the docs/images/*.png exports.

    Pure Wolfram Language; figures are returned, never written to disk.
    --------------------------------------------------------------------- *)

BeginPackage["CondensateFigures`"];

(* ---- Flory-Huggins thermodynamics (\[Section]2) ---- *)
fFH::usage          = "fFH[\[Phi],\[Chi],N1,N2] Flory-Huggins free energy density.";
muFH::usage         = "muFH[\[Phi],\[Chi],N1,N2] exchange chemical potential.";
d2fFH::usage        = "d2fFH[\[Phi],\[Chi],N1,N2] second derivative of f.";
criticalPoint::usage= "criticalPoint[N1,N2] -> <|\[Phi]c,\[Chi]c|>.";
spinodalCurve::usage= "spinodalCurve[N1,N2,npts] spinodal {\[Phi],\[Chi]} points.";
binodalCurve::usage = "binodalCurve[N1,N2,npts] binodal {\[Chi],\[Phi]L,\[Phi]R} points.";
figFHFreeEnergy::usage    = "figFHFreeEnergy[N1,N2] free-energy curves.";
figFHPhaseDiagram::usage  = "figFHPhaseDiagram[N1,N2] \[Chi]-\[Phi] phase diagram.";
figFHChainLength::usage   = "figFHChainLength[Nvals] chain-length family.";

(* ---- Cahn-Hilliard (\[Section]3) ---- *)
solveCahnHilliard1D::usage = "solveCahnHilliard1D[] -> InterpolatingFunction (double-well Model B).";
figCH1DEvolution::usage    = "figCH1DEvolution[] space-time + profiles of spinodal decomposition.";

(* ---- Kinetics (\[Section]4-5) ---- *)
figRateEnhancement::usage = "figRateEnhancement[] K^2 / K^3 rate enhancement.";
figRibozyme::usage        = "figRibozyme[kDeg] replication vs degradation vs K.";
figPrebioticSpace::usage  = "figPrebioticSpace[kDeg] net-replication K-c phase space.";

(* ---- Disease (\[Section]7) ---- *)
agingModel::usage          = "agingModel[kGel,kAgg,tMax] three-state NDSolve solution.";
figDiseaseTransition::usage= "figDiseaseTransition[] liquid->gel->aggregate panels.";
figDiseasePhaseDiagram::usage = "figDiseasePhaseDiagram[] material-state phase diagram.";

(* ---- Protocell dynamics (\[Section]6) ---- *)
figProtoOstwald::usage       = "figProtoOstwald[] Ostwald ripening of 100 droplets.";
figProtoSelection::usage     = "figProtoSelection[] proto-Darwinian selection trends.";
figProtoGrowthDivision::usage= "figProtoGrowthDivision[] protocell life-cycle schematic.";

(* ---- Nonlinear dynamics (\[Section]8) ---- *)
figDispersion::usage         = "figDispersion[] Cahn-Hilliard dispersion relation.";
figModeSelection::usage      = "figModeSelection[] selected wavelength over (\[Phi]0,\[Kappa]).";
figNLPhaseDiagram::usage     = "figNLPhaseDiagram[] FH phase diagram + double-well pitchfork.";
figReactionEnhancement::usage= "figReactionEnhancement[] mass-conserving rate enhancement.";
figEnergyLyapunov::usage     = "figEnergyLyapunov[] 1D CH: free energy decays monotonically.";

(* ---- Active condensates (\[Section]9) ---- *)
figActiveDispersion::usage   = "figActiveDispersion[] reaction-modified dispersion.";
figDispersionSurface::usage  = "figDispersionSurface[] 3D dispersion surface.";
figStickers::usage           = "figStickers[] stickers-and-spacers phase diagram.";
figReaction1D::usage         = "figReaction1D[] 1D reaction-CH: arrested coarsening.";
figThermalFission::usage     = "figThermalFission[] 1D thermal-gradient droplet fission.";

(* ---- Emergent complexity (\[Section]10) ---- *)
figMulticomponent::usage     = "figMulticomponent[] multi-component spinodals (Qian-Knowles).";
figCompartmentalization::usage="figCompartmentalization[] ternary free energy with two minima.";
figSelectionChemostat::usage = "figSelectionChemostat[] chemostat competition: fitness rises.";
figLineage::usage            = "figLineage[] protocell phylogeny coloured by fitness.";
figTuring::usage             = "figTuring[] Schnakenberg Turing growth rate inside vs outside.";
figWetDry::usage             = "figWetDry[] 1D wet-dry Cahn-Hilliard space-time.";
figParamSpace::usage         = "figParamSpace[] (\[Chi],N) two-phase + (K,p) rate enhancement.";

Begin["`Private`"];

(* shared accent colours *)
accent  = RGBColor[0.78, 0.12, 0.36];
accent2 = RGBColor[0.12, 0.44, 0.70];
gold    = RGBColor[0.91, 0.64, 0.24];
ink     = GrayLevel[0.1];
seq     = "SunsetColors";

(* ================================================================
   Flory-Huggins
   ================================================================ *)

fFH[\[Phi]_, \[Chi]_, N1_, N2_] :=
   \[Phi] Log[\[Phi]]/N1 + (1 - \[Phi]) Log[1 - \[Phi]]/N2 + \[Chi] \[Phi] (1 - \[Phi]);

muFH[\[Phi]_, \[Chi]_, N1_, N2_] :=
   (Log[\[Phi]] + 1)/N1 - (Log[1 - \[Phi]] + 1)/N2 + \[Chi] (1 - 2 \[Phi]);

d2fFH[\[Phi]_, \[Chi]_, N1_, N2_] :=
   1/(N1 \[Phi]) + 1/(N2 (1 - \[Phi])) - 2 \[Chi];

piOsmotic[\[Phi]_, \[Chi]_, N1_, N2_] :=
   fFH[\[Phi], \[Chi], N1, N2] - \[Phi] muFH[\[Phi], \[Chi], N1, N2];

criticalPoint[N1_, N2_] := Module[{sqr = Sqrt[N2/N1]},
   <|"\[Phi]c" -> sqr/(1 + sqr),
     "\[Chi]c" -> (1/Sqrt[N1] + 1/Sqrt[N2])^2/2|>];

spinodalCurve[N1_, N2_, npts_: 300] := Module[
   {cp = criticalPoint[N1, N2], \[Phi]vals, \[Chi]vals},
   \[Phi]vals = Subdivide[0.001, 0.999, npts];
   \[Chi]vals = (1/(N1 #) + 1/(N2 (1 - #)))/2 & /@ \[Phi]vals;
   Select[Transpose[{\[Phi]vals, \[Chi]vals}], #[[2]] >= cp["\[Chi]c"] &]
];

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

figFHFreeEnergy[N1_: 100, N2_: 1] := Module[{cp, \[Chi]vals},
   cp = criticalPoint[N1, N2];
   \[Chi]vals = {0.3, 0.5, 1.0, 1.5, 2.5} cp["\[Chi]c"];
   Plot[
     Evaluate[Table[fFH[\[Phi], \[Chi], N1, N2], {\[Chi], \[Chi]vals}]],
     {\[Phi], 0.001, 0.999},
     PlotStyle -> (Directive[Thick, #] & /@ {
        Darker[Blue, 0.3], Blue, Darker[Green], Darker[Orange], Darker[Red]}),
     PlotLegends -> Placed[
        LineLegend[
          (Directive[Thick, #] & /@ {
             Darker[Blue, 0.3], Blue, Darker[Green], Darker[Orange], Darker[Red]}),
          (Row[{"\[Chi] = ", NumberForm[#, {3, 3}]}] & /@ \[Chi]vals),
          LegendMarkerSize -> 22,
          Background -> Directive[GrayLevel[0, 0.55]],
          LabelStyle -> Directive[White, 11]],
        Scaled[{0.52, 0.74}]],
     PlotLabel -> Style[Row[{"Flory\[Dash]Huggins free energy  (", Subscript["N", "1"],
        " = ", N1, ", ", Subscript["N", "2"], " = ", N2, ")"}], 14],
     PlotRange -> {Automatic, {-0.15, 0.02}},
     GridLines -> {{}, {0}}, GridLinesStyle -> Directive[Gray, Dashed],
     ImageSize -> 720, Frame -> True,
     FrameLabel -> {"Volume fraction \[Phi]", "Free energy density  f(\[Phi])/kT"}]
];

figFHPhaseDiagram[N1_: 100, N2_: 1] := Module[{cp, sp, bn, pSpin, pBin, pCrit},
   cp = criticalPoint[N1, N2]; sp = spinodalCurve[N1, N2]; bn = binodalCurve[N1, N2];
   pSpin = ListLinePlot[sp, PlotStyle -> Directive[Red, Dashed, Thick],
      PlotLegends -> {"Spinodal"}];
   pBin = If[Length[bn] > 0,
     ListLinePlot[{bn[[All, {2, 1}]], bn[[All, {3, 1}]]},
       PlotStyle -> Directive[Blue, Thick], PlotLegends -> {"Binodal", None}],
     Graphics[]];
   pCrit = ListPlot[{{cp["\[Phi]c"], cp["\[Chi]c"]}},
      PlotStyle -> {Black, PointSize[Large]}, PlotLegends -> {"Critical point"}];
   Show[pBin, pSpin, pCrit,
     PlotRange -> {{0, 1}, {0, cp["\[Chi]c"] * 3}}, Frame -> True,
     FrameLabel -> {"Volume fraction \[Phi]",
        "Flory\[Dash]Huggins parameter \[Chi]  (\[Proportional] 1/T)"},
     PlotLabel -> Style[Row[{"Phase diagram  (", Subscript["N", "1"],
        " = ", N1, ", ", Subscript["N", "2"], " = ", N2, ")"}], 14],
     ImageSize -> 640, GridLines -> Automatic, GridLinesStyle -> Directive[LightGray]]
];

figFHChainLength[Nvals_: {1, 5, 20, 100, 500}] := Module[{colors, plots},
   colors = ColorData["DarkRainbow"] /@ Subdivide[0, 1, Length[Nvals] - 1];
   plots = MapThread[
     Module[{N1 = #1, sp, T\[Phi]},
       sp = spinodalCurve[N1, 1, 200];
       T\[Phi] = {#[[1]], 1/#[[2]]} & /@ sp;
       ListLinePlot[T\[Phi], PlotStyle -> Directive[Thick, #2],
         PlotLegends -> {Row[{"N = ", N1}]}]] &,
     {Nvals, colors}];
   Show[plots, PlotRange -> {{0, 1}, {0, Automatic}}, Frame -> True,
     FrameLabel -> {"Volume fraction \[Phi]", "Reduced temperature T*"},
     PlotLabel -> Style["Chain length widens the two-phase region", 14],
     ImageSize -> 640, GridLines -> Automatic, GridLinesStyle -> Directive[LightGray]]
];

(* ================================================================
   Cahn-Hilliard (double-well Model B; the stiff FH solver freezes)
   ================================================================ *)

solveCahnHilliard1D[] := Module[
   {a = 1.0, b = 1.0, kappa = 1.0, M = 1.0, L = 80.0, T = 400.0, nPts = 256,
    mu, coeffs, icf},
   mu[\[Psi]_] := -a \[Psi] + b \[Psi]^3;
   SeedRandom[17];
   coeffs = RandomReal[{-1, 1}, 24];
   icf[x_] := 0.02 Sum[coeffs[[n]] Sin[2 Pi n x/L], {n, 1, 24}];
   {L, T, Quiet @ NDSolveValue[{
      D[\[Psi][x, t], t] == M D[mu[\[Psi][x, t]] - kappa D[\[Psi][x, t], x, x], x, x],
      \[Psi][x, 0] == icf[x],
      \[Psi][0, t] == \[Psi][L, t],
      Derivative[1, 0][\[Psi]][0, t] == Derivative[1, 0][\[Psi]][L, t]},
      \[Psi], {x, 0, L}, {t, 0, T},
      Method -> {"MethodOfLines",
        "SpatialDiscretization" -> {"TensorProductGrid",
          "MinPoints" -> nPts, "MaxPoints" -> nPts}},
      MaxSteps -> 500000, AccuracyGoal -> 5, PrecisionGoal -> 5]}
];

figCH1DEvolution[] := Module[
   {sol, L, T, phiOf, \[Phi]0 = 0.4, amp = 0.28, tSnap, densPlot, profPlot},
   {L, T, sol} = solveCahnHilliard1D[];
   phiOf[x_, t_] := Clip[\[Phi]0 + amp sol[x, t], {0., 0.8}];
   tSnap = {0, 8, 40, 150, 400};
   densPlot = DensityPlot[phiOf[x, t], {x, 0, L}, {t, 0, T},
      PlotRange -> {0, 0.8}, ColorFunction -> "TemperatureMap", PlotPoints -> 90,
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
      PlotPoints -> 200, MaxRecursion -> 2, Frame -> True,
      FrameLabel -> {"Position x", "Volume fraction \[Phi]"},
      PlotLabel -> Style["Composition profiles at selected times", 13],
      PlotRange -> {0, 0.85}, GridLines -> Automatic,
      GridLinesStyle -> Directive[LightGray], ImageSize -> 620];
   GraphicsColumn[{densPlot, profPlot}, Spacings -> 15]
];

(* ================================================================
   Kinetics
   ================================================================ *)

figRateEnhancement[] := Module[{Krange},
   Krange = Subdivide[1., 200., 300];
   LogPlot[{#^2, #^3} &[K] // Evaluate, {K, 1, 200},
     PlotStyle -> {Directive[Thick, Blue], Directive[Thick, Dashed, Red]},
     PlotLegends -> Placed[LineLegend[{Blue, Red},
        {"Bimolecular  (K\:00b2)", "Termolecular  (K\:00b3, illustrative)"}], Scaled[{0.3, 0.8}]],
     Frame -> True, FrameLabel -> {"Partition coefficient K", "Rate enhancement"},
     PlotLabel -> Style["Reaction-rate enhancement inside condensates (open-system limit)", 13],
     GridLines -> Automatic, GridLinesStyle -> Directive[LightGray], ImageSize -> 640]
];

figRibozyme[kDeg_: 0.5] := Module[{Vmax = 1.0, Km = 10.0, Sfixed = 1.0, Krange, vRep, vDeg},
   Krange = Subdivide[1., 200., 300];
   vRep = Table[Module[{Seff = K Sfixed, Veff = K Vmax}, Veff Seff/(Km + Seff)], {K, Krange}];
   vDeg = kDeg # Sfixed & /@ Krange;
   Show[
     ListLinePlot[{Transpose[{Krange, vRep}], Transpose[{Krange, vDeg}]},
       PlotStyle -> {Directive[Thick, Blue], Directive[Thick, Dashed, Red]},
       PlotLegends -> {"Replication rate", "Degradation rate"},
       Filling -> {1 -> {{2}, Directive[Opacity[0.15], Green]}}],
     Frame -> True, FrameLabel -> {"Partition coefficient K", "Rate (\[Mu]M/min)"},
     PlotLabel -> Style["Replication vs degradation: the origin-of-life window", 14],
     GridLines -> Automatic, GridLinesStyle -> Directive[LightGray], ImageSize -> 640]
];

figPrebioticSpace[kDeg_: 1.0] := Module[{Vmax = 1.0, Km = 10.0, Krange, Crange, netRate, pts},
   Krange = Subdivide[1., 200., 150];
   Crange = 10^Subdivide[-3., 1., 150];
   netRate = Table[Module[{Seff = K c, Veff = K Vmax, rep, deg},
       rep = Veff Seff/(Km + Seff); deg = kDeg K c; rep - deg], {K, Krange}, {c, Crange}];
   pts = Flatten[Table[{Krange[[i]], Log10[Crange[[j]]], netRate[[i, j]]},
       {i, Length[Krange]}, {j, Length[Crange]}], 1];
   Show[
     ListContourPlot[pts, Contours -> 20, ColorFunction -> "RedBlueTones",
       ContourLabels -> None,
       PlotLegends -> Placed[BarLegend[{"RedBlueTones", Automatic},
         LegendLabel -> "Net rate"], Right],
       FrameLabel -> {"Partition coefficient K",
          Row[{"log", Subscript["", "10"], "([Substrate]", Subscript["", "bulk"], " / \[Mu]M)"}]},
       PlotLabel -> Style["Net replication: replication \[Minus] degradation", 14]],
     ListContourPlot[pts, Contours -> {0}, ContourStyle -> Directive[Black, Thick],
       ContourLabels -> None],
     Frame -> True, ImageSize -> 640]
];

(* ================================================================
   Disease
   ================================================================ *)

agingModel[kGel_, kAgg_, tMax_: 100.] :=
   NDSolve[{l'[t] == -kGel l[t], g'[t] == kGel l[t] - kAgg g[t], a'[t] == kAgg g[t],
     l[0] == 1, g[0] == 0, a[0] == 0}, {l, g, a}, {t, 0, tMax}][[1]];

figDiseaseTransition[] := Module[{makePanel},
   makePanel[sol_, title_] := Plot[Evaluate[{l[t], g[t], a[t]} /. sol], {t, 0, 100},
      PlotStyle -> {Directive[Thick, Blue], Directive[Thick, Orange], Directive[Thick, Red]},
      Filling -> {1 -> Axis, 2 -> Axis, 3 -> Axis},
      FillingStyle -> {Directive[Opacity[0.15], LightBlue],
        Directive[Opacity[0.15], LightOrange], Directive[Opacity[0.15], LightRed]},
      PlotLegends -> Placed[{"Liquid", "Gel", "Aggregate"}, {0.8, 0.5}],
      Frame -> True, FrameLabel -> {"Time (a.u.)", "Fraction"},
      PlotLabel -> Style[title, 12], PlotRange -> {{0, 100}, {0, 1.05}},
      GridLines -> Automatic, GridLinesStyle -> Directive[LightGray], ImageSize -> 350];
   Labeled[
     GraphicsGrid[{{
       makePanel[agingModel[0.02, 0.005], "Slow maturation\n(wild type)"],
       makePanel[agingModel[0.05, 0.01], "Moderate\n(risk variant)"],
       makePanel[agingModel[0.1, 0.02], "Fast (disease mutant)\n(FUS/TDP-43)"]}},
       ImageSize -> 1050, Spacings -> 15],
     Style["Condensate ageing: liquid \[Rule] gel \[Rule] aggregate (neurodegeneration)",
       14, Bold], Top]
];

figDiseasePhaseDiagram[] := Module[{c, boundLL, boundGel, boundAgg},
   c = Subdivide[0., 10., 200];
   boundLL = 3.0 + 2.0 Exp[-0.5 #] & /@ c;
   boundGel = 2.0 + 1.0 Exp[-0.3 #] & /@ c;
   boundAgg = 1.0 + 0.5 Exp[-0.2 #] & /@ c;
   Show[
     ListLinePlot[{Transpose[{c, boundLL}], Transpose[{c, boundGel}], Transpose[{c, boundAgg}]},
       PlotStyle -> {Directive[Thick, Blue], Directive[Thick, Orange], Directive[Thick, Red]},
       Filling -> {1 -> {10, Directive[Opacity[0.1], LightBlue]},
         2 -> {{1}, Directive[Opacity[0.1], LightOrange]},
         3 -> {{2}, Directive[Opacity[0.1], LightRed]}},
       PlotLegends -> {"Liquid condensate", "Gel transition", "Solid aggregate"}],
     Graphics[{
       Text[Style["Mixed\n(one phase)", 13, Bold, Blue], {5, 8}],
       Text[Style["Liquid condensate", 12, Bold, Darker[Orange]], {7, 4.2}],
       Text[Style["Gel", 12, Bold, Darker[Red, 0.3]], {7, 2.5}],
       Text[Style["Solid aggregate", 11, Bold, Darker[Red]], {7, 0.7}],
       {Red, Thick, Arrowheads[0.04], Arrow[{{2, 4.5}, {4, 3}}]},
       Text[Style["Disease mutations\nshift boundaries", 10, Italic, Red], {1.5, 5.2}]}],
     Frame -> True,
     FrameLabel -> {"Protein concentration (a.u.)", "Temperature / Interaction strength (a.u.)"},
     PlotLabel -> Style["Condensate material state phase diagram", 14],
     PlotRange -> {{0, 10}, {0, 10}}, GridLines -> Automatic,
     GridLinesStyle -> Directive[LightGray], ImageSize -> 640]
];

(* ================================================================
   Protocell dynamics (\[Section]6)
   ================================================================ *)

ostwaldRipening[nDroplets_, tMax_, dt_, seed_] := Module[
   {R, Rstar, dRdt, alive, t = 0., lc = 1.0, D0 = 1.0, history = {}},
   SeedRandom[seed];
   R = Exp[RandomVariate[NormalDistribution[0, 0.5], nDroplets]];
   R = Clip[R, {0.01, \[Infinity]}];
   AppendTo[history, {0., R}];
   While[t < tMax && Count[R, x_ /; x > 0.01] > 2,
     alive = UnitStep[R - 0.01];
     Rstar = Total[R alive]/Max[1, Total[alive]];
     dRdt = alive D0 lc (1/Rstar - 1/Clip[R, {0.01, \[Infinity]}])/Clip[R, {0.01, \[Infinity]}];
     R += dRdt dt; R = Clip[R, {0, \[Infinity]}]; R = R UnitStep[R - 0.01];
     t += dt;
     If[Mod[Round[t/dt], Round[1/dt]] == 0, AppendTo[history, {t, R}]]];
   history
];

figProtoOstwald[] := Module[{history, rd},
   history = ostwaldRipening[100, 100., 0.01, 42];
   rd = {#[[1]], Mean[Select[#[[2]], # > 0.01 &]], Max[#[[2]]],
      Count[#[[2]], x_ /; x > 0.01]} & /@ history;
   GraphicsGrid[{{
     ListLinePlot[{rd[[All, {1, 2}]], rd[[All, {1, 3}]]},
       PlotStyle -> {Directive[Thick, Blue], Directive[Thick, Dashed, Red]},
       PlotLegends -> {"Mean R", "Max R"}, Frame -> True,
       FrameLabel -> {"Time", "Radius"}, PlotLabel -> Style["Droplet coarsening", 13],
       GridLines -> Automatic, GridLinesStyle -> Directive[LightGray], ImageSize -> 400],
     ListLinePlot[rd[[All, {1, 4}]], PlotStyle -> Directive[Thick, Darker[Green]],
       Frame -> True, FrameLabel -> {"Time", "Surviving droplets"},
       PlotLabel -> Style["Ostwald ripening", 13], GridLines -> Automatic,
       GridLinesStyle -> Directive[LightGray], ImageSize -> 400]}}, ImageSize -> 800]
];

protocellSim[nInit_, tMax_, dt_, cBulk_, seed_] := Module[
   {cells, t = 0., history = {}, maxPop = 300, newCells, fitness},
   SeedRandom[seed];
   cells = Table[{RandomReal[{0.5, 2.0}], RandomReal[{10, 200}],
      RandomReal[{0.01, 0.2}], 0., 0}, {nInit}];
   While[t < tMax && Length[cells] > 0,
     cells = Map[Module[
         {r = #[[1]], K = #[[2]], kR = #[[3]], age = #[[4]], gen = #[[5]], cEff, netRate, vol},
         cEff = K cBulk; netRate = kR cEff - 0.05 cEff; vol = 4/3 Pi r^3;
         If[netRate > 0, vol += netRate vol dt 0.01]; r = (3 vol/(4 Pi))^(1/3);
         {r, K, kR, age + dt, gen}] &, cells];
     newCells = {};
     Do[If[cell[[1]] > 3.0,
         Module[{rNew = cell[[1]]/2^(1/3), K = cell[[2]], kR = cell[[3]], gen = cell[[5]], Kn, kn},
           Kn = 1 + 0.05 RandomReal[NormalDistribution[]];
           kn = 1 + 0.05 RandomReal[NormalDistribution[]];
           AppendTo[newCells, {rNew, K Kn, kR kn, 0., gen + 1}];
           AppendTo[newCells, {rNew, K/Kn, kR/kn, 0., gen + 1}]],
         AppendTo[newCells, cell]], {cell, cells}];
     cells = Select[newCells, #[[1]] > 0.1 &];
     If[Length[cells] > maxPop, fitness = #[[2]] #[[3]] & /@ cells;
       cells = cells[[Ordering[fitness, -maxPop]]]];
     If[Mod[Round[t/dt], 10] == 0 && Length[cells] > 0,
       AppendTo[history, {t, Length[cells], Mean[cells[[All, 2]]],
         Mean[cells[[All, 3]]], Max[cells[[All, 5]]]}]];
     t += dt];
   history
];

figProtoSelection[] := Module[{h},
   h = protocellSim[50, 500., 0.5, 0.1, 42];
   GraphicsGrid[{{
     ListLinePlot[h[[All, {1, 3}]], PlotStyle -> Directive[Thick, Red], Frame -> True,
       FrameLabel -> {"Time", "Mean K"}, PlotLabel -> Style["Selection for higher K", 13],
       GridLines -> Automatic, GridLinesStyle -> Directive[LightGray], ImageSize -> 350],
     ListLinePlot[h[[All, {1, 4}]], PlotStyle -> Directive[Thick, Darker[Green]], Frame -> True,
       FrameLabel -> {"Time", "Mean replication rate"},
       PlotLabel -> Style["Selection for faster replication", 13],
       GridLines -> Automatic, GridLinesStyle -> Directive[LightGray], ImageSize -> 350],
     ListLinePlot[h[[All, {1, 5}]], PlotStyle -> Directive[Thick, Purple], Frame -> True,
       FrameLabel -> {"Time", "Max generation"}, PlotLabel -> Style["Lineage depth", 13],
       GridLines -> Automatic, GridLinesStyle -> Directive[LightGray], ImageSize -> 350]}},
     ImageSize -> 1000]
];

figProtoGrowthDivision[] := Graphics[{
   {EdgeForm[Directive[Thick, Blue]], FaceForm[Directive[Opacity[0.3], LightBlue]], Disk[{0, 0}, 0.5]},
   Text[Style["Small\ncondensate", 10], {0, -0.9}],
   {Thick, Darker[Green], Arrowheads[0.03], Arrow[{{0.7, 0}, {2.3, 0}}]},
   Text[Style["Growth\n(accretion)", 9, Italic, Darker[Green]], {1.5, 0.3}],
   {EdgeForm[Directive[Thick, Blue]], FaceForm[Directive[Opacity[0.3], LightBlue]], Disk[{3.5, 0}, 1.0]},
   Text[Style["Large\ncondensate", 10], {3.5, -1.4}],
   {Thick, Red, Arrowheads[0.03], Arrow[{{4.7, 0}, {6.3, 0}}]},
   Text[Style["Division\n(threshold)", 9, Italic, Red], {5.5, 0.35}],
   {EdgeForm[Directive[Thick, Blue]], FaceForm[Directive[Opacity[0.3], LightBlue]], Disk[{7.5, 0.5}, 0.6]},
   {EdgeForm[Directive[Thick, Blue]], FaceForm[Directive[Opacity[0.3], LightBlue]], Disk[{7.5, -0.5}, 0.6]},
   Text[Style["Daughters\n(with variation)", 10], {7.5, -1.4}],
   {Thick, Orange, Arrowheads[0.03], Arrow[{{8.3, 0.5}, {9.8, 0.5}}]},
   {Thick, Gray, Dashed, Arrowheads[0.03], Arrow[{{8.3, -0.5}, {9.0, -1.2}}]},
   Text[Style["Fitter\nsurvives", 9, Italic, Orange], {9.3, 0.8}],
   Text[Style["Weaker\ndissolves", 9, Italic, Gray], {9.5, -1.0}],
   {EdgeForm[Directive[Thick, Darker[Green]]], FaceForm[Directive[Opacity[0.3], LightGreen]], Disk[{10.5, 0.5}, 0.55]}},
   PlotRange -> {{-1.2, 11.5}, {-2, 1.8}}, ImageSize -> 800,
   PlotLabel -> Style["Protocell life cycle: growth \[RightArrow] division \[RightArrow] selection", 14]];

(* ================================================================
   Nonlinear dynamics (\[Section]8)
   ================================================================ *)

cN1 = 100; cN2 = 1;
cphiC = Sqrt[cN2/cN1]/(1 + Sqrt[cN2/cN1]);
cchiC = 0.5 (1/Sqrt[cN1] + 1/Sqrt[cN2])^2;
cfpp[phi_, chi_] := 1/(cN1 phi) + 1/(cN2 (1 - phi)) - 2 chi;
comega[k_, phi_, chi_, kappa_, M_] := -M k^2 (cfpp[phi, chi] + kappa k^2);
ckStar[phi_, chi_, kappa_] := Sqrt[-cfpp[phi, chi]/(2 kappa)];

figDispersion[] := Module[{chi = 1.5, kappa = 0.5, M = 1.0, phis, curves, dots},
   phis = {0.06, 0.10, 0.20, 0.30, 0.45, 0.60};
   curves = Table[comega[k, p, chi, kappa, M], {p, phis}];
   dots = Table[If[cfpp[p, chi] < 0,
       {ckStar[p, chi, kappa], comega[ckStar[p, chi, kappa], p, chi, kappa, M]}, Nothing], {p, phis}];
   Show[
     Plot[Evaluate@curves, {k, 0, 3},
       PlotStyle -> (Directive[Thick, ColorData["SunsetColors"][#]] & /@ Rescale[Range@Length@phis]),
       PlotLegends -> Placed[Map[Row[{"\!\(\*SubscriptBox[\(\[Phi]\), \(0\)]\) = ", #}] &, phis], {0.18, 0.30}],
       Frame -> True, GridLines -> Automatic, GridLinesStyle -> Directive[LightGray, Dotted],
       FrameLabel -> {"wavenumber  k", "growth rate  \[Omega](k)"},
       PlotLabel -> Style["Cahn\[Dash]Hilliard dispersion  \[Omega](k) = \[Minus]M k\:00b2(f''(\[Phi]\:2080) + \[Kappa] k\:00b2)", 14],
       PlotRange -> {{0, 3}, {-1.6, 1.6}}, ImageSize -> 620],
     Graphics[{ink, PointSize[0.012], Point /@ dots, Text[Style["fastest mode k*", 10, ink], {1.7, 1.25}]}],
     Graphics[{Gray, Line[{{0, 0}, {3, 0}}]}]]
];

figModeSelection[] := Module[{chi = 1.5},
   Show[
     DensityPlot[If[cfpp[phi, chi] < 0, Min[2 Pi/ckStar[phi, chi, kappa], 25.], Indeterminate],
       {phi, 0.005, 0.72}, {kappa, 0.05, 2.0}, ColorFunction -> "SunsetColors", PlotPoints -> 90,
       PlotRange -> All, PlotLegends -> Automatic,
       FrameLabel -> {"mean composition  \[Phi]\:2080", "gradient coefficient  \[Kappa]"},
       PlotLabel -> Style["selected wavelength  \[Lambda]* = 2\[Pi]/k*  (white = stable, f'' > 0)", 14],
       ImageSize -> 620],
     ContourPlot[cfpp[phi, chi] == 0, {phi, 0.005, 0.72}, {kappa, 0.05, 2.0},
       ContourStyle -> Directive[White, Thick]]]
];

figNLPhaseDiagram[] := Module[{spin, gg, hh, mu, Pizero, chiOf, binL, binR, panelFH, ddwA, ddwB, panelDW},
   spin[phi_] := 0.5 (1/(cN1 phi) + 1/(cN2 (1 - phi)));
   gg[p_] := (Log[p] + 1)/cN1 - (Log[1 - p] + 1)/cN2;
   hh[p_] := p Log[p]/cN1 + (1 - p) Log[1 - p]/cN2;
   mu[p_, chi_] := gg[p] + chi (1 - 2 p);
   Pizero[p_, chi_] := hh[p] - p gg[p] + chi p^2;
   chiOf[pR_] := (pR gg[pR] - hh[pR])/pR^2;
   With[{tris = Table[Module[{chi = chiOf[pR], muc, u, pL},
        muc = mu[pR, chi];
        u = u /. Quiet@FindRoot[mu[Exp[u], chi] == muc, {u, -8, -700, Log[cphiC] - 0.01}];
        pL = Exp[u];
        If[NumberQ[pL] && cchiC < chi <= cchiC*3 && 0 < pL < cphiC, {chi, pL, pR}, Nothing]],
        {pR, Subdivide[cphiC + 0.003, 0.93, 120]}]},
     binL = SortBy[{#[[2]], #[[1]]} & /@ tris, Last];
     binR = SortBy[{#[[3]], #[[1]]} & /@ tris, Last]];
   panelFH = Show[
     Plot[spin[phi], {phi, 0.004, 0.75}, PlotStyle -> Directive[accent, Thick, Dashed],
       Frame -> True, GridLines -> Automatic, GridLinesStyle -> Directive[LightGray, Dotted],
       FrameLabel -> {"volume fraction  \[Phi]", "\[Chi]  (\[Proportional] 1/T)"},
       PlotRange -> {{0, 0.8}, {0, cchiC*3}}, ImageSize -> 430],
     ListLinePlot[{binL, binR}, PlotStyle -> Directive[accent2, Thick]],
     Graphics[{ink, PointSize[0.02], Point[{cphiC, cchiC}],
        Text[Style["critical point", 10], {cphiC + 0.20, cchiC + 0.12}],
        Directive[RGBColor[0.85, 0.6, 0.1]], PointSize[0.028], Point[{0.30, 1.5}],
        Text[Style["solver operating point", 9], {0.45, 1.62}],
        Text[Style["two-phase", 11, accent2], {0.45, 1.1}],
        Text[Style["one phase", 11, Darker@Green], {0.42, 0.4}]}],
     PlotLabel -> Style["(a) Flory\[Dash]Huggins phase diagram", 13]];
   ddwA = Plot[{Sqrt[a], -Sqrt[a], Sqrt[a/3], -Sqrt[a/3]}, {a, 0, 1.5},
     PlotStyle -> {Directive[accent2, Thick], Directive[accent2, Thick],
       Directive[accent, Thick, Dashed], Directive[accent, Thick, Dashed]},
     Frame -> True, GridLines -> Automatic, GridLinesStyle -> Directive[LightGray, Dotted],
     FrameLabel -> {"control parameter  a  (\[Proportional] \[Chi] \[Minus] \[Chi]\:1D9C)",
        "order parameter  \[Phi] \[Minus] \[Phi]\:1D9C"},
     PlotRange -> {{0, 1.5}, {-1.35, 1.35}}, ImageSize -> 430,
     PlotLabel -> Style["(b) double-well pitchfork", 13],
     Filling -> {3 -> {4}}, FillingStyle -> Directive[Opacity[0.12, accent]]];
   ddwB = Graphics[{ink, PointSize[0.02], Point[{0, 0}],
     Text[Style["binodal", 10, accent2], {1.15, 1.18}], Text[Style["spinodal", 10, accent], {1.2, 0.45}]}];
   panelDW = Show[ddwA, ddwB];
   GraphicsRow[{panelFH, panelDW}, Spacings -> 20, ImageSize -> 920]
];

figReactionEnhancement[] := Module[{enh, p0 = 0.05, heat, curves, Kvals},
   enh[KA_, KB_, p_] := (p KA KB + (1 - p))/((1 + p (KA - 1)) (1 + p (KB - 1)));
   heat = DensityPlot[Log10[enh[10^lKA, 10^lKB, p0]], {lKA, 0, 3}, {lKB, 0, 3},
     ColorFunction -> "SunsetColors", PlotPoints -> 70, PlotRange -> All,
     FrameLabel -> {"\!\(\*SubscriptBox[\(log\), \(10\)]\) \!\(\*SubscriptBox[\(K\), \(A\)]\)",
        "\!\(\*SubscriptBox[\(log\), \(10\)]\) \!\(\*SubscriptBox[\(K\), \(B\)]\)"},
     PlotLegends -> Placed[BarLegend[Automatic,
        LegendLabel -> "\!\(\*SubscriptBox[\(log\), \(10\)]\)(\!\(\*SubscriptBox[\(k\), \(eff\)]\)/\!\(\*SubscriptBox[\(k\), \(out\)]\))"], Right],
     PlotLabel -> Style["(a) enhancement over (\!\(\*SubscriptBox[\(K\), \(A\)]\), \!\(\*SubscriptBox[\(K\), \(B\)]\)),  p = 0.05", 13],
     ImageSize -> 430];
   Kvals = {10, 30, 100, 300, 1000};
   curves = LogPlot[Evaluate@Table[enh[Kk, Kk, p], {Kk, Kvals}], {p, 0.0005, 1},
     PlotStyle -> (Directive[Thick, ColorData["TemperatureMap"][#]] & /@ Rescale[Range@Length@Kvals]),
     PlotLegends -> Placed[Map[Row[{"K = ", #}] &, Kvals], {0.82, 0.72}],
     Frame -> True, GridLines -> Automatic, GridLinesStyle -> Directive[LightGray, Dotted],
     FrameLabel -> {"dense-phase area fraction  p", "\!\(\*SubscriptBox[\(k\), \(eff\)]\)/\!\(\*SubscriptBox[\(k\), \(out\)]\)"},
     PlotLabel -> Style["(b) optimum at p* \[TildeTilde] 1/K  (gain \[TildeTilde] K/4, not K\:00b2)", 13],
     PlotRange -> {{0, 1}, {1, 300}}, ImageSize -> 430];
   GraphicsRow[{heat, curves}, Spacings -> 20, ImageSize -> 920]
];

figEnergyLyapunov[] := Module[
   {a = 1.0, b = 1.0, kappa = 1.0, M = 1.0, L = 80.0, T = 400.0, mu, sol, coeffs, icf, nPts = 256, Ftot, ts, Fs, spacetime},
   mu[ph_] := -a ph + b ph^3;
   SeedRandom[7]; coeffs = RandomReal[{-1, 1}, 24];
   icf[x_] := 0.02 Sum[coeffs[[n]] Sin[2 Pi n x/L], {n, 1, 24}];
   sol = Quiet@NDSolveValue[{
       D[ph[x, t], t] == M D[mu[ph[x, t]] - kappa D[ph[x, t], x, x], x, x],
       ph[x, 0] == icf[x], ph[0, t] == ph[L, t],
       Derivative[1, 0][ph][0, t] == Derivative[1, 0][ph][L, t]},
     ph, {x, 0, L}, {t, 0, T},
     Method -> {"MethodOfLines", "SpatialDiscretization" -> {"TensorProductGrid", "MinPoints" -> nPts, "MaxPoints" -> nPts}},
     MaxSteps -> 500000, AccuracyGoal -> 5, PrecisionGoal -> 5];
   Ftot[t_?NumericQ] := NIntegrate[-a/2 sol[x, t]^2 + b/4 sol[x, t]^4 + kappa/2 Derivative[1, 0][sol][x, t]^2, {x, 0, L}, AccuracyGoal -> 4];
   ts = Subdivide[0.5, T, 60]; Fs = Ftot /@ ts;
   spacetime = DensityPlot[sol[x, t], {x, 0, L}, {t, 0, T}, ColorFunction -> "BalancedHue",
     PlotPoints -> 90, PlotRange -> {-1, 1}, FrameLabel -> {"position  x", "time  t"},
     PlotLabel -> Style["(a) 1D spinodal decomposition + coarsening", 13],
     PlotLegends -> Automatic, AspectRatio -> 0.7, ImageSize -> 430];
   GraphicsRow[{spacetime,
     ListLinePlot[Transpose[{ts, Fs}], PlotStyle -> Directive[accent, Thick], Filling -> Bottom,
       FillingStyle -> Directive[Opacity[0.12, accent]], Frame -> True, GridLines -> Automatic,
       GridLinesStyle -> Directive[LightGray, Dotted], FrameLabel -> {"time  t", "free energy  F[\[Phi]]"},
       PlotLabel -> Style["(b) Lyapunov functional F[\[Phi]] decreases", 13], ImageSize -> 430]},
     Spacings -> 20, ImageSize -> 920]
];

(* ================================================================
   Active and chemically-maintained condensates (\[Section]9)
   ================================================================ *)

figActiveDispersion[] := Module[{fpp = -1.0, kappa = 1.0, M = 1.0, krs},
   krs = {0.0, 0.05, 0.15, 0.25};
   Show[
     Plot[Evaluate@Table[-M k^2 (fpp + kappa k^2) - kr, {kr, krs}], {k, 0, 1.4},
       PlotStyle -> (Directive[Thick, ColorData["SunsetColors"][#]] & /@ Rescale[Range@Length@krs]),
       PlotLegends -> Placed[Map[Row[{"\!\(\*SubscriptBox[\(k\), \(react\)]\) = ", #}] &, krs], {0.25, 0.28}],
       Frame -> True, GridLines -> Automatic, GridLinesStyle -> Directive[LightGray, Dotted],
       FrameLabel -> {"wavenumber  k", "growth rate  \[Omega](k)"},
       PlotLabel -> Style["Reaction stabilises long-wavelength modes \[LongDash] arrested coarsening", 14],
       PlotRange -> {{0, 1.4}, {-0.35, 0.30}}, ImageSize -> 620],
     Graphics[{Gray, Line[{{0, 0}, {1.4, 0}}],
        Text[Style["finite unstable band \[Rule] selected length scale", 10, accent], {0.7, 0.22}]}]]
];

figDispersionSurface[] := Module[{chi = 1.5, nN1 = 100, nN2 = 1, kappa = 0.5, M = 1.0, fpp},
   fpp[phi_] := 1/(nN1 phi) + 1/(nN2 (1 - phi)) - 2 chi;
   Plot3D[Clip[-M k^2 (fpp[phi] + kappa k^2), {-1.5, 1.5}], {phi, 0.02, 0.66}, {k, 0, 2.2},
     PlotPoints -> 60, MaxRecursion -> 2, ColorFunction -> "SunsetColors", PlotRange -> {-1.0, 1.4},
     MeshFunctions -> {#3 &}, Mesh -> {{0}}, MeshStyle -> Directive[Cyan, Thick],
     AxesLabel -> {"\!\(\*SubscriptBox[\(\[Phi]\), \(0\)]\)", "k", "\[Omega]"},
     PlotLabel -> Style["Dispersion surface  \[Omega](k, \[Phi]\:2080)  (cyan ridge: neutral stability)", 14],
     ViewPoint -> {-2.4, -1.6, 1.4}, Lighting -> "Neutral", ImageSize -> 640, BoxRatios -> {1, 1, 0.7}]
];

figStickers[] := Module[{pbond, fbond, fbpp, chiSpin, cases, curves, h = 10.^-4},
   pbond[phis_, Ka_] := With[{x = Ka Max[phis, 10.^-12]}, ((2 x + 1) - Sqrt[4 x + 1])/(2 x)];
   fbond[phi_, v_, Ka_] := If[Ka <= 0, 0.,
      Module[{phis = v phi, p}, p = pbond[phis, Ka]; phis (Log[Max[1 - p, 10.^-12]] + p/2)]];
   fbpp[phi_, v_, Ka_] := (fbond[phi + h, v, Ka] - 2 fbond[phi, v, Ka] + fbond[phi - h, v, Ka])/h^2;
   chiSpin[phi_, v_, Ka_] := 0.5 (1/phi + 1/(1 - phi) + If[Ka > 0, fbpp[phi, v, Ka], 0.]);
   cases = {{1, 0, "no stickers (FH)"}, {2, 2, "v=2, \!\(\*SubscriptBox[\(K\), \(a\)]\)=2"},
            {3, 2, "v=3, \!\(\*SubscriptBox[\(K\), \(a\)]\)=2"}, {4, 4, "v=4, \!\(\*SubscriptBox[\(K\), \(a\)]\)=4"}};
   curves = Table[With[{v = c[[1]], Ka = c[[2]]}, Table[{phi, chiSpin[phi, v, Ka]}, {phi, 0.02, 0.98, 0.005}]], {c, cases}];
   ListLinePlot[curves, PlotStyle -> (Directive[Thick, ColorData["SunsetColors"][#]] & /@ Rescale[Range@Length@cases]),
     PlotLegends -> Placed[cases[[All, 3]], {0.78, 0.78}], Frame -> True, GridLines -> Automatic,
     GridLinesStyle -> Directive[LightGray, Dotted], FrameLabel -> {"volume fraction  \[Phi]", "interaction  \[Chi]"},
     PlotLabel -> Style["Stickers-and-spacers: multivalency expands the two-phase region", 14],
     PlotRange -> {{0, 1}, {0, 3.2}}, ImageSize -> 620]
];

figReaction1D[] := Module[{a = 1.0, b = 1.0, kappa = 1.0, M = 1.0, L = 100.0, T = 600.0, coeffs, icf, mu, solve, panels},
   mu[ph_] := -a ph + b ph^3; SeedRandom[11]; coeffs = RandomReal[{-1, 1}, 30];
   icf[x_] := 0.03 Sum[coeffs[[n]] Sin[2 Pi n x/L], {n, 1, 30}];
   solve[kr_] := Quiet@NDSolveValue[{
       D[ph[x, t], t] == M D[mu[ph[x, t]] - kappa D[ph[x, t], x, x], x, x] - kr (ph[x, t] - 0),
       ph[x, 0] == icf[x], ph[0, t] == ph[L, t], Derivative[1, 0][ph][0, t] == Derivative[1, 0][ph][L, t]},
     ph, {x, 0, L}, {t, 0, T},
     Method -> {"MethodOfLines", "SpatialDiscretization" -> {"TensorProductGrid", "MinPoints" -> 200, "MaxPoints" -> 200}},
     MaxSteps -> 500000, AccuracyGoal -> 5, PrecisionGoal -> 5];
   panels = Table[With[{sol = solve[kr]},
       If[Head[sol] === InterpolatingFunction,
         DensityPlot[sol[x, t], {x, 0, L}, {t, 0, T}, ColorFunction -> "BalancedHue", PlotPoints -> 90,
           PlotRange -> {-1.2, 1.2}, FrameLabel -> {"position  x", "time  t"},
           PlotLabel -> Style[If[kr == 0, "(a) passive: coarsens to one domain",
              "(b) active (" <> ToString[kr] <> "): arrested pattern"], 13], AspectRatio -> 1, ImageSize -> 430],
         Graphics[Text["NDSolve failed"], ImageSize -> 430]]], {kr, {0.0, 0.03}}];
   GraphicsRow[panels, Spacings -> 20, ImageSize -> 920]
];

(* Thermal-gradient fission (\[Section]9.4): native-Wolfram replacement for the
   former Python figure. A position-dependent well a(x) with a central warm
   band (a<0) dissolves a droplet's neck and splits it in two. *)
figThermalFission[] := Module[{b = 1.0, kappa = 1.0, M = 1.0, L = 60.0, T = 140.0, nx = 160, icf, solve, pairs, panels},
   icf[x_] := 0.85 (Tanh[(x - 16)/2] - Tanh[(x - 44)/2] - 1);
   solve[aF_] := Quiet@NDSolveValue[{
       D[ph[x, t], t] == M D[-aF[x] ph[x, t] + b ph[x, t]^3 - kappa D[ph[x, t], x, x], x, x],
       ph[x, 0] == icf[x], ph[0, t] == ph[L, t], Derivative[1, 0][ph][0, t] == Derivative[1, 0][ph][L, t]},
     ph, {x, 0, L}, {t, 0, T},
     Method -> {"MethodOfLines", "SpatialDiscretization" -> {"TensorProductGrid", "MinPoints" -> nx, "MaxPoints" -> nx}},
     MaxSteps -> 500000, AccuracyGoal -> 5, PrecisionGoal -> 5];
   pairs = {{Function[x, 1.0], "(a) uniform temperature: domain persists"},
            {Function[x, 1.0 - 1.9 Exp[-((x - 30)/4)^2]], "(b) warm band at x=30: neck dissolves \[Rule] two daughters"}};
   panels = Table[Module[{aF = pr[[1]], lbl = pr[[2]], sol},
       sol = solve[aF];
       If[Head[sol] === InterpolatingFunction,
         DensityPlot[sol[x, t], {x, 0, L}, {t, 0, T}, ColorFunction -> "BalancedHue", PlotPoints -> 80,
           PlotRange -> {-1.1, 1.1}, FrameLabel -> {"position  x", "time  t"},
           PlotLabel -> Style[lbl, 12.5], AspectRatio -> 1, ImageSize -> 430],
         Graphics[Text["NDSolve failed"], ImageSize -> 430]]], {pr, pairs}];
   GraphicsRow[panels, Spacings -> 20, ImageSize -> 920]
];

(* ================================================================
   Emergent complexity (\[Section]10)
   ================================================================ *)

figMulticomponent[] := Module[{Ns = {2, 3, 5, 10}, chiS, phiCf, chiCf, cols, curves, pts},
   chiS[phi_, n_] := 2 (1/(n phi) + 1/(1 - phi)); phiCf[n_] := 1/(1 + Sqrt[n]); chiCf[n_] := 2 (1 + Sqrt[n])^2/n;
   cols = Table[ColorData[seq][0.12 + 0.62 (i - 1)/(Length[Ns] - 1)], {i, Length[Ns]}];
   curves = Table[chiS[phi, Ns[[i]]], {i, Length[Ns]}];
   pts = Table[{phiCf[n], chiCf[n]}, {n, Ns}];
   Show[
     Plot[Evaluate@curves, {phi, 0.01, 0.99}, PlotStyle -> (Directive[#, Thick] & /@ cols), Frame -> True,
       GridLines -> Automatic, GridLinesStyle -> Directive[LightGray, Dotted],
       FrameLabel -> {"total condensed fraction  \[Phi]", "interaction  \[Chi]  (Qian\[Dash]Knowles)"},
       PlotRange -> {{0, 1}, {0, 18}}, PlotLegends -> Placed[(Style["N=" <> ToString[#], 11] & /@ Ns), {0.82, 0.7}],
       PlotLabel -> Style["Multi-component spinodals: critical point migrates as 1/(1+\!\(\*SqrtBox[\(N\)]\))", 13], ImageSize -> 560],
     Graphics[{Black, PointSize[0.022], Point[pts], White, PointSize[0.012], Point[pts]}]]
];

(* Ternary compartmentalisation (\[Section]10.1): native-Wolfram replacement.
   Two mutually-repelling solutes (chiAB large) make an A-rich and a B-rich
   condensate -- the physical basis of organelle identity. *)
figCompartmentalization[] := Module[{ff, chiAB = 3.5, chiAS = 0.5, chiBS = 0.5},
   ff[a_, b_] := a Log[a] + b Log[b] + (1 - a - b) Log[1 - a - b] + chiAB a b + chiAS a (1 - a - b) + chiBS b (1 - a - b);
   ContourPlot[ff[a, b], {a, 0.02, 0.94}, {b, 0.02, 0.94},
     RegionFunction -> Function[{a, b}, a + b < 0.97], Contours -> 22, ColorFunction -> "SunsetColors",
     PlotPoints -> 60, FrameLabel -> {"component A fraction", "component B fraction"},
     PlotLabel -> Style["Ternary free energy: two immiscible minima \[Rule] distinct condensates", 12.5],
     Epilog -> {ink, PointSize[0.022], Point[{0.7, 0.08}], Point[{0.08, 0.7}],
        Text[Style["A-rich\ncondensate", 10, ink], {0.66, 0.24}],
        Text[Style["B-rich\ncondensate", 10, ink], {0.24, 0.66}]}, ImageSize -> 560]
];

figSelectionChemostat[] := Module[{r0 = 0.6, death = 0.16, supply = 0.22, use = 0.05, gs, nC, eqs, ic, sol, tmax = 120, meanFit, fitPlot, cPlot},
   gs = {0.6, 0.85, 1.1, 1.35, 1.6, 1.8}; nC = Length[gs];
   eqs = Join[Table[Subscript[c, i]'[t] == r0 gs[[i]] Subscript[c, i][t] nut[t] - death Subscript[c, i][t], {i, nC}],
     {nut'[t] == supply (1 - nut[t]) - use nut[t] Sum[r0 gs[[j]] Subscript[c, j][t], {j, nC}]}];
   ic = Join[Table[Subscript[c, i][0] == 0.1, {i, nC}], {nut[0] == 1.0}];
   sol = NDSolve[Join[eqs, ic], Append[Table[Subscript[c, i], {i, nC}], nut], {t, 0, tmax}][[1]];
   meanFit[tt_] := (Sum[gs[[i]] (Subscript[c, i][tt] /. sol), {i, nC}])/(Sum[(Subscript[c, i][tt] /. sol), {i, nC}]);
   fitPlot = Plot[meanFit[t], {t, 0, tmax}, PlotStyle -> Directive[accent, Thickness[0.01]], Frame -> True,
     GridLines -> Automatic, GridLinesStyle -> Directive[LightGray, Dotted],
     FrameLabel -> {"time", "mean fitness  \[LeftAngleBracket]g\[RightAngleBracket]"}, PlotRange -> {{0, tmax}, {0.55, 1.85}},
     Epilog -> {gold, Dashed, Line[{{0, Max[gs]}, {tmax, Max[gs]}}], Text[Style["g_max", 10, gold], {tmax 0.9, Max[gs] - 0.05}]},
     PlotLabel -> Style["Chemostat selection: mean fitness rises", 13], ImageSize -> 430];
   cPlot = Plot[Evaluate[Table[Subscript[c, i][t] /. sol, {i, nC}]], {t, 0, tmax},
     PlotStyle -> Table[Directive[ColorData[seq][0.1 + 0.8 (i - 1)/(nC - 1)], Thick], {i, nC}], Frame -> True,
     GridLines -> Automatic, GridLinesStyle -> Directive[LightGray, Dotted],
     FrameLabel -> {"time", "genotype abundance  \!\(\*SubscriptBox[\(c\), \(i\)]\)"},
     PlotLegends -> Placed[(Style["g=" <> ToString[#], 9] & /@ gs), Right],
     PlotLabel -> Style["The fittest replicator wins", 13], ImageSize -> 470];
   GraphicsRow[{fitPlot, cPlot}, Spacings -> 20, ImageSize -> 940]
];

(* Protocell phylogeny (\[Section]10.2): native-Wolfram replacement for the
   former Python lineage figure. Fitter lineages branch more; node colour = fitness. *)
figLineage[] := Module[{edges = {}, fitOf = <|0 -> 1.0|>, alive = {0}, nid = 1, newAlive, nKids, cid},
   SeedRandom[5];
   Do[
     newAlive = {};
     Do[
       nKids = RandomChoice[{Max[0.1, 0.65 - 0.45 (fitOf[p] - 1)], 0.5, Max[0.05, 0.45 (fitOf[p] - 0.6)]} -> {0, 1, 2}];
       Do[cid = nid++; AppendTo[fitOf, cid -> Clip[fitOf[p] + RandomReal[{-0.1, 0.16}], {0.5, 2.0}]];
          AppendTo[edges, p -> cid]; AppendTo[newAlive, cid], {nKids}];
       If[nKids == 0 && RandomReal[] < 0.5, AppendTo[newAlive, p]], {p, alive}];
     alive = If[Length[newAlive] > 60, RandomSample[newAlive, 60], newAlive];
     If[alive === {}, Break[]], {7}];
   Graph[edges,
     VertexStyle -> (# -> ColorData["TemperatureMap"][Rescale[fitOf[#], {0.5, 2.0}]] & /@ Keys[fitOf]),
     VertexSize -> Medium, GraphLayout -> "LayeredEmbedding", ImageSize -> 620,
     PlotLabel -> Style["Protocell phylogeny: fit (yellow) lineages branch; low-fitness (blue) ones die out", 12.5]]
];

figTuring[] := Module[{a = 0.1, b = 1.05, us, vs, fu, fv, gu, gv, growth, DuIn = 0.30, DvIn = 9.0, DuOut = 2.4, DvOut = 14.0, gIn = 3.4, gOut = 1.4, kmax = 6},
   us = a + b; vs = b/us^2; fu = -1 + 2 us vs; fv = us^2; gu = -2 us vs; gv = -us^2;
   growth[k_, Du_, Dv_, gam_] := Module[{J, tr, det}, J = gam {{fu, fv}, {gu, gv}} - k^2 {{Du, 0}, {0, Dv}}; tr = Tr[J]; det = Det[J]; (tr + Sqrt[tr^2 - 4 det])/2];
   Plot[{Re@growth[k, DuIn, DvIn, gIn], Re@growth[k, DuOut, DvOut, gOut]}, {k, 0, kmax},
     PlotStyle -> {Directive[accent, Thickness[0.01]], Directive[accent2, Thickness[0.008], Dashed]}, Frame -> True,
     GridLines -> Automatic, GridLinesStyle -> Directive[LightGray, Dotted], FrameLabel -> {"wavenumber  k", "growth rate  \[Lambda](k)"},
     PlotRange -> {{0, kmax}, {-2.5, 1.2}},
     PlotLegends -> Placed[{"inside droplet (dense): unstable band", "dilute phase: all modes stable"}, {0.62, 0.25}],
     Epilog -> {Gray, Line[{{0, 0}, {kmax, 0}}]},
     PlotLabel -> Style["Schnakenberg Turing instability: patterns only inside the condensate", 12.5], ImageSize -> 620]
];

figWetDry[] := Module[{L = 60.0, M = 1.0, b = 1.0, kappa = 1.0, a0 = 0.6, dA = 1.4, period = 60.0, tmax = 180.0, aF, sol, nx = 160},
   aF[t_] := a0 + dA Sin[2 Pi t/period]; SeedRandom[7];
   sol = Quiet@Check[NDSolveValue[{
       D[phi[t, x], t] == M D[-aF[t] phi[t, x] + b phi[t, x]^3 - kappa D[phi[t, x], {x, 2}], {x, 2}],
       phi[0, x] == 0.05 RandomReal[{-1, 1}] + 0.04 Sum[Sin[2 Pi n x/L + n], {n, 3, 11}],
       phi[t, 0] == phi[t, L], Derivative[0, 1][phi][t, 0] == Derivative[0, 1][phi][t, L],
       Derivative[0, 2][phi][t, 0] == Derivative[0, 2][phi][t, L], Derivative[0, 3][phi][t, 0] == Derivative[0, 3][phi][t, L]},
     phi, {t, 0, tmax}, {x, 0, L},
     Method -> {"MethodOfLines", "SpatialDiscretization" -> {"TensorProductGrid", "MinPoints" -> nx, "MaxPoints" -> nx, "DifferenceOrder" -> 4}}], $Failed];
   If[sol === $Failed, Graphics[Text["NDSolve failed"], ImageSize -> 680],
     DensityPlot[sol[t, x], {t, 0, tmax}, {x, 0, L}, PlotPoints -> 90, ColorFunction -> "RedBlueTones", AspectRatio -> 0.6,
       Frame -> True, FrameLabel -> {"time  t", "position  x"},
       PlotLabel -> Style["1D wet\[Dash]dry Cahn\[Dash]Hilliard: domains dissolve and re-form each cycle", 12.5],
       ImageSize -> 680, PlotLegends -> BarLegend[Automatic]]]
];

figParamSpace[] := Module[{chiCf, twoPhase, ratePlot},
   chiCf[n_] := (1/2) (1 + 1/Sqrt[n])^2;
   twoPhase = RegionPlot[chi > chiCf[n], {chi, 0, 3.5}, {n, 1, 60}, PlotPoints -> 60, BoundaryStyle -> Directive[Cyan, Thick],
     PlotStyle -> accent, Frame -> True, FrameLabel -> {"interaction  \[Chi]", "chain length  N"},
     PlotLabel -> Style["(\[Chi],N): one-phase / two-phase", 12.5],
     Epilog -> {Text[Style["two-phase", 12, White], {2.6, 45}], Text[Style["one-phase", 12, ink], {0.6, 45}]}, ImageSize -> 430];
   ratePlot = DensityPlot[(p K^2 + (1 - p))/(1 + p (K - 1))^2, {K, 1, 300}, {p, 0.001, 0.999},
     ScalingFunctions -> {"Log", None}, PlotPoints -> 70, ColorFunction -> "SunsetColors", Frame -> True,
     FrameLabel -> {"partition coefficient  K", "dense fraction  p"},
     PlotLabel -> Style["(K,p): rate enhancement \!\(\*SubscriptBox[\(k\), \(eff\)]\)/\!\(\*SubscriptBox[\(k\), \(out\)]\)", 12.5],
     Epilog -> {White, Dashed, Thick, Line[Table[{K, 1/K}, {K, 1, 300, 1}]]}, PlotLegends -> BarLegend[Automatic], ImageSize -> 470];
   GraphicsRow[{twoPhase, ratePlot}, Spacings -> 24, ImageSize -> 940]
];

End[];
EndPackage[];
