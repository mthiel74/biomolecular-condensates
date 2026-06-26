#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Master runner: execute all Wolfram analysis scripts in order, then
   copy the Python-generated figures from figures/ into docs/images/ so
   the notebook builder finds every figure it embeds.

   NOTE: the Python track (src/*.py) must be run first — it writes the
   darwin_*, multicomp_*, active_*, rd_*, wetdry_* and param_space_grand
   figures into figures/. See the Reproducibility section of the notebook.

   Usage:  wolframscript -file wolfram/run_all.wls
*)

repoRoot = ParentDirectory @ DirectoryName[$InputFileName];
scripts  = {
   "phase_diagram.wls",
   "cahn_hilliard.wls",
   "kinetics.wls",
   "protocell.wls",
   "disease.wls",
   "nonlinear_dynamics.wls",
   "active_dynamics.wls",
   "phase3_extensions.wls"
};

(* Python-generated figures the notebook embeds (written to figures/). *)
pythonFigures = {
   "darwin_fitness.png", "darwin_generations.png", "darwin_lineage.png",
   "multicomp_phase_diagram.png", "multicomp_compartmentalization.png",
   "active_morphology_distribution.png", "active_thermal_fission.png",
   "rd_interface_hotspots.png", "rd_turing_patterns.png",
   "wetdry_population.png", "wetdry_morphology.png", "param_space_grand.png"
};

Print["=" <> StringRepeat["=", 59]];
Print["  Biomolecular condensates — full Wolfram pipeline"];
Print["=" <> StringRepeat["=", 59]];
Print[""];

startAll = AbsoluteTime[];

Do[
   Module[{path, t0, exitCode},
     path = FileNameJoin[{repoRoot, "wolfram", script}];
     Print["\n>>> Running ", script, " ..."];
     t0 = AbsoluteTime[];
     exitCode = Run["wolframscript -file " <> path];
     Print["<<< ", script, " done (",
       NumberForm[AbsoluteTime[] - t0, {5, 1}], " s, exit ", exitCode, ")"];
   ],
   {script, scripts}];

(* Copy the Python-generated figures into docs/images/ *)
Print["\n>>> Copying Python-generated figures into docs/images/ ..."];
Module[{src, dst, figDir, imgDir},
   figDir = FileNameJoin[{repoRoot, "figures"}];
   imgDir = FileNameJoin[{repoRoot, "docs", "images"}];
   Do[
     src = FileNameJoin[{figDir, f}];
     dst = FileNameJoin[{imgDir, f}];
     If[FileExistsQ[src],
       CopyFile[src, dst, OverwriteTarget -> True];
       Print["    copied ", f],
       Print["    WARNING: missing ", FileNameJoin[{"figures", f}],
         " \[LongDash] run the Python track first (see Reproducibility)"]];
     ,
     {f, pythonFigures}];
];

Print["\n" <> StringRepeat["=", 60]];
Print["  All scripts complete in ",
   NumberForm[AbsoluteTime[] - startAll, {6, 1}], " s"];
Print[StringRepeat["=", 60]];
