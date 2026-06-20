#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Master runner: execute all analysis scripts in order.
   Writes all figures to docs/images/.

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
   "active_dynamics.wls"
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

Print["\n" <> StringRepeat["=", 60]];
Print["  All scripts complete in ",
   NumberForm[AbsoluteTime[] - startAll, {6, 1}], " s"];
Print[StringRepeat["=", 60]];
