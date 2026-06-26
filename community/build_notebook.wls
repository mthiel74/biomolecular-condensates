#!/usr/bin/env wolframscript
(* ::Package:: *)

(*  Build community/biomolecular_condensates.nb — the educational
    long-form notebook for the Wolfram Community submission.

    Architecture: identical to the ENSO-emergence and
    scalinglawsincities builders. Accumulates Cell expressions,
    writes the Notebook expression to disk. Images embedded as
    compressed RasterBox data.

    Sections
    --------
      Abstract
      1. The membrane-free organelle
      2. Thermodynamics of phase separation
         2.1 The Flory\[Dash]Huggins free energy
         2.2 Spinodal and binodal: the phase diagram
         2.3 Chain length and multivalency
      3. Cahn\[Dash]Hilliard dynamics
         3.1 The equation
         3.2 Spinodal decomposition
         3.3 Coarsening: Lifshitz-Slyozov
      4. Concentration enhancement and catalysis
         4.1 The partition coefficient
         4.2 Bimolecular rate enhancement
         4.3 Prebiotic chemistry phase space
      5. RNA world inside condensates
         5.1 Ribozyme kinetics
         5.2 Replication vs degradation
         5.3 The origin-of-life window
      6. Protocell dynamics
         6.1 Growth and division
         6.2 Proto-Darwinian selection
      7. Condensates and disease
         7.1 Liquid-to-gel-to-aggregate
         7.2 Phase diagram of material states
      8. Nonlinear dynamics: stability, bifurcation, coarsening
      9. Active and chemically-maintained condensates
     10. Emergent complexity: selection, patterns and cycles
     11. Conclusions
     12. References
        Further reading (New Scientist)
     13. Reproducibility
*)

repoRoot = ParentDirectory @ DirectoryName[$InputFileName];
imgDir   = FileNameJoin[{repoRoot, "docs", "images"}];
outNb    = FileNameJoin[{repoRoot, "community", "biomolecular_condensates.nb"}];
outPdf   = FileNameJoin[{repoRoot, "community", "biomolecular_condensates.pdf"}];


(* ==================================================================
   Helpers — identical to the ENSO/scaling notebook builders
   ================================================================== *)

title[txt_]    := Cell[txt, "Title"];
subtitle[txt_] := Cell[txt, "Subtitle"];
hd1[txt_]      := Cell[txt, "Section"];
hd2[txt_]      := Cell[txt, "Subsection"];
hd3[txt_]      := Cell[txt, "Subsubsection"];

(* Body text and abstract inherit the RiverNetworkStatistics stylesheet
   set on the Notebook (narrow margins, gray Abstract). No per-cell
   margin/background overrides, so the stylesheet governs the look. *)
para[stuff_]   := Cell[TextData[stuff], "Text"];

abstractCell[stuff_] := Cell[TextData[stuff], "Abstract"];

(* Opening hero caption: small gray text, matching the RiverNetwork
   post's lead-in caption (FontSize 12, GrayLevel 0.4). *)
heroCaption[stuff_] := Cell[TextData[stuff], "Text",
   FontSize -> 12, FontColor -> GrayLevel[0.4]];

bold[s_]   := StyleBox[s, FontWeight -> Bold];
ital[s_]   := StyleBox[s, FontSlant -> Italic];
math[s_]   := StyleBox[s, FontFamily -> "Times", FontSlant -> Italic];
mono[s_]   := StyleBox[s, FontFamily -> "Courier"];

imath[boxes_] := Cell[BoxData @ FormBox[boxes, TraditionalForm],
   "InlineFormula"];
sub[base_, k_] := imath @ SubscriptBox[base, k];

displayMath[expr_] := Cell[
   BoxData @ FormBox[expr, TraditionalForm],
   "DisplayFormula",
   CellMargins -> {{80, 50}, {8, 8}}
];

(* Image embedding *)
imgCacheDir = FileNameJoin[{repoRoot, "community", ".img-cache"}];
If[!DirectoryQ[imgCacheDir], CreateDirectory[imgCacheDir]];

(* Embedded raster width in pixels. High-resolution source PNGs
   (docs/images/*.png are typically 1700-4200 px wide) are preserved
   up to this cap; only genuinely huge panoramas are downsampled. At a
   ~560 pt display size this gives ~180 DPI, a big improvement over the
   old 520 px cap (which rendered at ~58 DPI). *)
targetWidth[path_] := 1400;

resizedCopy[path_] := Module[{cached, w, h, im, target},
   cached = FileNameJoin[{imgCacheDir, FileNameTake[path]}];
   If[FileExistsQ[cached] &&
      FileDate[cached, "Modification"] > FileDate[path, "Modification"],
     Return[cached]];
   im = Import[path];
   {w, h} = ImageDimensions[im];
   target = targetWidth[path];
   If[w > target, im = ImageResize[im, target]];
   Export[cached, im, "PNG"];
   cached
];

$rasterSubstitutions = <||>;
$rasterCounter = 0;

imageBoxes[path_String] := Module[
   {small, im, na, h, w, compressed, key, boxStr},
   small = resizedCopy[path];
   im = Image[Import[small], "Byte"];
   na = Reverse @ ImageData[im, "Byte"];
   {h, w} = Take[Dimensions[na], 2];
   compressed = Compress[na];
   $rasterCounter += 1;
   key = "RasterPh" <> ToString[$rasterCounter] <> "X";
   boxStr = StringJoin[
      "GraphicsBox[",
        "TagBox[",
          "RasterBox[",
            "CompressedData[\"", compressed, "\"], ",
            "{{0, 0}, {", ToString[w], ", ", ToString[h], "}}, ",
            "{0, 255}, ",
            "ColorFunction -> RGBColor",
          "],",
          "BoxForm`ImageTag[\"Byte\", ColorSpace -> \"RGB\", Interleaving -> True],",
          " Selectable -> False",
        "],",
        " DefaultBaseStyle -> \"ImageGraphics\",",
        " ImageSize -> Full,",
        " ImageSizeRaw -> {", ToString[w], ", ", ToString[h], "},",
        " PlotRange -> {{0, ", ToString[w], "}, {0, ", ToString[h], "}}",
      "]"
   ];
   $rasterSubstitutions[key] = boxStr;
   key
];

imgCell[file_, opts___] := Module[{path = FileNameJoin[{imgDir, file}], key},
   If[FileExistsQ[path],
     key = imageBoxes[path];
     Cell[BoxData @ Symbol[key],
        "Output",
        TextAlignment -> Center,
        CellMargins -> {{Automatic, Automatic}, {10, 6}}],
     Cell[
       TextData[{bold["[ missing figure: " <> file <> " ]"]}],
       "Text", FontColor -> RGBColor[0.7, 0.1, 0.1],
       TextAlignment -> Center]
   ]
];

captionCell[args__] := Cell[
   TextData[{args}], "Text",
   FontSlant -> Italic, FontSize -> 11,
   FontColor -> GrayLevel[0.35],
   CellMargins -> {{70, 70}, {2, 14}}
];

codeIn[code_String] := Cell[code, "Program",
   CellMargins -> {{60, 30}, {6, 6}},
   FontFamily -> "Courier", FontSize -> 11,
   ShowStringCharacters -> False,
   Background -> RGBColor[0.95, 0.96, 0.98]];

(* Live, evaluatable Wolfram Language input cell. Unlike codeIn (a
   static "Program" listing), these are real "Input" cells the reader
   can evaluate in place. *)
wlInput[code_String] := Cell[BoxData[code], "Input"];

tryNote[stuff_] := Cell[TextData[stuff], "Text",
   FontSlant -> Italic, FontColor -> RGBColor[0.153, 0.51, 0.64]];

link[txt_String, url_String] := Sequence[
   StyleBox[txt, FontWeight -> "SemiBold"],
   " (", StyleBox[url, FontFamily -> "Courier", FontSize -> 10], ")"
];

(* Cell accumulator *)
$allCells = {};
writeAll[cs_List] := Module[{},
   $allCells = Join[$allCells, cs];
   Print["  writeAll: +", Length[cs], " cells (total = ", Length[$allCells], ")"];
];


(* ============================================================
   TITLE & ABSTRACT
   ============================================================ *)

writeAll[{

title["Biomolecular Condensates and the Origin of Life"],

subtitle["From Flory\[Dash]Huggins thermodynamics to protocellular \
compartments \[LongDash] first-principles mathematical modelling of \
liquid\[Dash]liquid phase separation in biology"],

Cell["Marco Thiel", "Subsubtitle", FontSlant -> Italic],

Cell["Wolfram Community post \[LongDash] June 2026",
   "Subsubtitle", FontSlant -> Italic, FontColor -> GrayLevel[0.4]],

(* ABSTRACT *)

hd1["Abstract"],

abstractCell[{"Biomolecular condensates are ", bold["membrane-free \
liquid droplets"], " inside cells, formed by liquid\[Dash]liquid \
phase separation (LLPS) of proteins and nucleic acids. Discovered as \
P granules in ", ital["C. elegans"], " germ cells (Brangwynne et al., \
2009), they have since been found to organise virtually every aspect of \
cellular function: gene regulation, signal transduction, stress \
response, and the formation of structures from the nucleolus to \
synaptic densities."}],

abstractCell[{"This notebook develops the mathematical physics of \
condensates from first principles. We implement the ",
   bold["Flory\[Dash]Huggins"], " free energy to derive the phase \
diagram, solve the ", bold["Cahn\[Dash]Hilliard"], " partial \
differential equation for spinodal decomposition, compute the ",
   bold["concentration enhancement"], " that makes condensates natural \
chemical reactors, and model the ", bold["protocell dynamics"],
   " (growth, division, competition) that connect condensates to the \
origin of life. We conclude with the ", bold["liquid \[Rule] gel \
\[Rule] aggregate"], " transition that links condensate dysfunction \
to neurodegenerative disease."}],

abstractCell[{ital["Audience:"], " the prose is written for two \
overlapping audiences. Each section opens with an intuitive explanation \
that requires no background, then deepens into a quantitative treatment \
with equations, numerical integration, and visualisations. Biophysicists \
will recognise the canonical references; readers from other backgrounds \
should be able to follow the physical narrative."}],

abstractCell[{"All code and figures are original and open-source: ",
   link["github.com/mthiel74/biomolecular-condensates",
     "https://github.com/mthiel74/biomolecular-condensates"], "."}],

para[{bold["Roadmap. "], "\[Section]2\[Dash]3 build the equilibrium and \
dynamical theory (Flory\[Dash]Huggins phase diagrams, Cahn\[Dash]Hilliard \
spinodal decomposition and coarsening). \[Section]4\[Dash]6 turn condensates \
into chemistry: concentration enhancement, an RNA-world replication window, \
and protocell growth, division and proto-Darwinian selection. \[Section]7 \
covers the liquid\[Rule]gel\[Rule]aggregate disease transition. \
\[Section]8\[Dash]10 add the nonlinear-dynamics toolkit (stability, \
bifurcation, Lyapunov functionals), active/chemically-maintained \
condensates, and the emergent layer where selection, Turing patterns and \
wet\[Dash]dry cycles arise. \[Section]11 concludes."}],

(* Initialisation cell, just after the abstract. It must precede the hero
   figure so the figure package is loaded before figParamSpace is called. *)
Cell[TextData[{"Initialisation \[LongDash] evaluate this cell first. Every \
figure below is generated by the code cell shown just above it, using the \
model functions in ", StyleBox["wolfram/condensate_figures.wl",
   FontFamily -> "Courier"], " (which sits beside this notebook in the \
repository)."}], "Text", FontSize -> 12, FontColor -> GrayLevel[0.4]],

wlInput["Module[{dir = Quiet@NotebookDirectory[]},\n  Get@FileNameJoin[{If[StringQ[dir], ParentDirectory[dir], " <> ToString[repoRoot, InputForm] <> "], \"wolfram\", \"condensate_figures.wl\"}]]"],

wlInput["figParamSpace[]"],

heroCaption[{"The phase behaviour of biomolecular condensates, generated live \
here: (left) the interaction\[Dash]chain-length plane that decides whether a \
mixture demixes at all, and (right) how strongly a condensate accelerates the \
reactions inside it. Beneath it all runs a single idea \[LongDash] the same \
Flory\[Dash]Huggins physics that unmixes oil from water also organises the \
living cell, and may have lit the first spark of life."}]
}];


(* ============================================================
   \[Section]1 The membrane-free organelle
   ============================================================ *)

writeAll[{

hd1["1.  The membrane-free organelle"],

para[{"For over a century, cell biology was built on the assumption that \
the cell's internal organisation depends on ", bold["membranes"],
   ": lipid bilayers that enclose the nucleus, the mitochondria, the \
endoplasmic reticulum. Compartmentalisation by membranes was the \
organising principle of eukaryotic life."}],

para[{"Then, in 2009, ", bold["Brangwynne, Eckmann, Courson, Rybarska, \
Hoege, Gharakhani, J\[UDoubleDot]licher, and Hyman"], " showed that \
P granules \[LongDash] RNA\[Dash]protein bodies in the germ cells of \
the nematode ", ital["C. elegans"], " \[LongDash] are not membrane-bound \
organelles at all. They are ", bold["liquid droplets"], ". They flow, \
they fuse, they drip, they dissolve and re-form. They are phase-separated \
condensates of RNA and proteins, held together not by a lipid membrane \
but by the same physics that separates oil from water."}],

para[{"This discovery launched a revolution. Within a decade, dozens of \
cellular structures were reclassified as condensates: stress granules, \
nucleoli, Cajal bodies, PML bodies, centrosomes, synaptic densities, \
heterochromatin domains, transcriptional hubs. Protein science has been \
transformed in parallel \[LongDash] the 2024 Nobel Prize in Chemistry went \
to David Baker for computational protein design and jointly to Demis \
Hassabis and John Jumper for protein-structure prediction \[LongDash] and \
the condensate revolution adds a complementary lesson: that proteins are \
not only rigid machines but can phase-separate into functional liquids."}],

para[{ital["The physics."], " The formation of a condensate is a \
thermodynamic phase transition. Below a critical temperature (or above \
a critical interaction strength), a homogeneous solution of \
macromolecules becomes thermodynamically unstable and spontaneously \
separates into a dilute phase and a dense phase. The dense phase is the \
condensate. The transition is governed by a competition between the \
entropy of mixing (which favours homogeneity) and the enthalpic \
interactions between molecules (which favour like-with-like \
aggregation). The mathematical framework for this competition is \
", bold["Flory\[Dash]Huggins theory"], ", which we derive in the \
next section."}],

para[{ital["The origin-of-life connection."], " Perhaps the most \
provocative implication of the condensate revolution is its connection \
to the ", bold["origin of life"], ". Aleksandr Oparin proposed (in 1924, \
and in molecular detail in his 1936 book) that life began inside ",
   bold["coacervate droplets"], " \[LongDash] \
spontaneously forming colloidal structures that concentrate organic \
molecules from the dilute primordial ocean. Oparin's hypothesis was \
dismissed for decades because the known coacervates were too disordered \
to sustain biochemistry. But modern condensates, with their exquisite \
sequence-dependent phase behaviour and their ability to concentrate \
specific clients by factors of up to a thousand, revive Oparin's vision \
in molecular detail. We show in \[Section]4 that condensate-assisted \
chemistry solves the ", bold["concentration problem"], " that has \
bedevilled prebiotic chemistry: how do you get reactions going in a \
picomolar ocean?"}]
}];


(* ============================================================
   \[Section]2 Thermodynamics of phase separation
   ============================================================ *)

writeAll[{

hd1["2.  Thermodynamics of phase separation"],

hd2["2.1  The Flory\[Dash]Huggins free energy"],

para[{ital["The lattice model."], " Consider a lattice of ", math["N\:209C"],
   " sites, each occupied by either a monomer of species 1 (the polymer, \
chain length ", math["N\:2081"], ") or a monomer of species 2 (the \
solvent, with ", math["N\:2082 = 1"], "). Let ",
   math["\[Phi]"], " be the volume fraction of species 1. The \
Flory\[Dash]Huggins free energy of mixing per lattice site is:"}],

displayMath[
   RowBox[{
     FractionBox["f", "kT"], "=",
     FractionBox[RowBox[{"\[Phi]", " ", RowBox[{"ln", " ", "\[Phi]"}]}],
        SubscriptBox["N", "1"]],
     "+",
     FractionBox[RowBox[{"(", RowBox[{"1", "-", "\[Phi]"}], ")",
        " ", RowBox[{"ln", "(", RowBox[{"1", "-", "\[Phi]"}], ")"}]}],
        SubscriptBox["N", "2"]],
     "+",
     RowBox[{"\[Chi]", " ", "\[Phi]", "(", RowBox[{"1", "-", "\[Phi]"}], ")"}]
   }]
],

para[{"The first two terms are the ", bold["entropy of mixing"],
   " \[LongDash] always negative (favour mixing) and scaled by the \
inverse chain lengths ", math["1/N\:2081"], " and ", math["1/N\:2082"],
   ". Longer chains have less translational entropy per monomer, so \
polymers phase-separate more readily than small molecules. The third \
term is the ", bold["enthalpy of mixing"], ", parameterised by the \
Flory interaction parameter ", math["\[Chi]"],
   ". For many condensate-forming proteins, \[Chi] > 0 (unfavourable \
mixing) and \[Chi] decreases with temperature (\[Chi] \[Proportional] \
1/T), giving an upper critical solution temperature (UCST) phase \
diagram \[LongDash] though a substantial class of intrinsically \
disordered proteins (e.g. elastin-like sequences) instead show the \
opposite, lower critical solution temperature (LCST), behaviour."}],

wlInput["figFHFreeEnergy[100, 1]"],
captionCell["Figure 2.1  Flory\[Dash]Huggins free energy density for a \
polymer\[Dash]solvent system (N\:2081 = 100, N\:2082 = 1) at several \
values of \[Chi]. Below the critical \[Chi], the free energy has a \
single minimum (stable mixed state). Above it, a double well appears: \
the system phase-separates into dilute (\[Phi]\:2097) and dense \
(\[Phi]\:1D48) phases at the two minima."],

tryNote[{"Evaluate the cell below to reproduce Figure 2.1 \[LongDash] the \
whole of Flory\[Dash]Huggins is two lines of Wolfram Language:"}],

wlInput["fFH[\[Phi]_, n1_, n2_, \[Chi]_] := \[Phi] Log[\[Phi]]/n1 + (1 - \[Phi]) Log[1 - \[Phi]]/n2 + \[Chi] \[Phi] (1 - \[Phi]);\nPlot[Evaluate@Table[fFH[\[Phi], 100, 1, \[Chi]], {\[Chi], {0.4, 0.605, 0.9, 1.5}}], {\[Phi], 0.001, 0.999},\n  PlotRange -> {-0.04, 0.04}, Frame -> True, FrameLabel -> {\"\[Phi]\", \"f/kT\"},\n  PlotLegends -> {\"\[Chi]=0.4\", \"\[Chi]=0.605 (critical)\", \"\[Chi]=0.9\", \"\[Chi]=1.5\"}]"],

tryNote[{"And the double well opening up as \[Chi] crosses the critical \
value is best felt interactively (drag \[Chi] through \[TildeTilde] 0.6):"}],

wlInput["Manipulate[\n Plot[\[Phi] Log[\[Phi]]/n1 + (1 - \[Phi]) Log[1 - \[Phi]] + \[Chi] \[Phi] (1 - \[Phi]),\n  {\[Phi], 0.001, 0.999}, PlotRange -> All, Frame -> True,\n  FrameLabel -> {\"\[Phi]\", \"f/kT\"}, ImageSize -> 480],\n {{\[Chi], 0.605, \"\[Chi]\"}, 0.3, 1.6, 0.005},\n {{n1, 100, \"chain length N\"}, 1, 300, 1}]"],

hd2["2.2  Spinodal and binodal: the phase diagram"],

para[{"The thermodynamic stability of the mixed state is determined by \
the curvature of the free energy. The ", bold["spinodal"], " is the \
locus of inflection points:"}],

displayMath[
   RowBox[{
     FractionBox[RowBox[{SuperscriptBox["\[PartialD]", "2"], "f"}],
        RowBox[{"\[PartialD]", SuperscriptBox["\[Phi]", "2"]}]],
     "=",
     FractionBox["1", RowBox[{SubscriptBox["N", "1"], " ", "\[Phi]"}]],
     "+",
     FractionBox["1", RowBox[{SubscriptBox["N", "2"], " ", "(", RowBox[{"1", "-", "\[Phi]"}], ")"}]],
     "-",
     RowBox[{"2", "\[Chi]"}],
     "=", "0"
   }]
],

para[{"Inside the spinodal, the free energy has negative curvature (",
   math["\[PartialD]\:00b2f/\[PartialD]\[Phi]\:00b2 < 0"],
   ") and the system is ", bold["unstable"], " to infinitesimal \
composition fluctuations \[LongDash] this is ", bold["spinodal \
decomposition"], ". Between the spinodal and the ", bold["binodal"],
   " (the coexistence curve), the system is ", bold["metastable"],
   " \[LongDash] it requires a finite-amplitude fluctuation (nucleation) \
to phase-separate."}],

para[{"The binodal is found by the ", bold["common tangent construction"],
   ": the two coexisting compositions \[Phi]\:2097 and \[Phi]\:1D48 \
satisfy equal chemical potential and equal osmotic pressure:"}],

displayMath[
   RowBox[{
     RowBox[{"\[Mu]", "(", SubscriptBox["\[Phi]", "L"], ")"}],
     "=",
     RowBox[{"\[Mu]", "(", SubscriptBox["\[Phi]", "R"], ")"}],
     ",", "     ",
     RowBox[{"\[CapitalPi]", "(", SubscriptBox["\[Phi]", "L"], ")"}],
     "=",
     RowBox[{"\[CapitalPi]", "(", SubscriptBox["\[Phi]", "R"], ")"}]
   }]
],

para[{"The ", bold["critical point"], " is where the spinodal and \
binodal meet:"}],

displayMath[
   RowBox[{
     SubscriptBox["\[Phi]", "c"], "=",
     FractionBox[SqrtBox[SubscriptBox["N", "2"]],
        RowBox[{SqrtBox[SubscriptBox["N", "1"]], "+",
           SqrtBox[SubscriptBox["N", "2"]]}]],
     ",", "     ",
     SubscriptBox["\[Chi]", "c"], "=",
     FractionBox[SuperscriptBox[RowBox[{"(",
        FractionBox["1", SqrtBox[SubscriptBox["N", "1"]]], "+",
        FractionBox["1", SqrtBox[SubscriptBox["N", "2"]]], ")"}], "2"], "2"]
   }]
],

para[{"For a long polymer (", math["N\:2081 = 100"], ", ",
   math["N\:2082 = 1"], "): ",
   math["\[Phi]\:1D9C = 0.091"], ", ",
   math["\[Chi]\:1D9C = 0.605"], ". The critical concentration is low \
(~9%), meaning that even dilute polymer solutions phase-separate when \
\[Chi] exceeds a modest threshold."}],

wlInput["figFHPhaseDiagram[100, 1]"],
captionCell["Figure 2.2  Phase diagram in \[Chi]\[Dash]\[Phi] \
coordinates. The blue curve is the binodal (coexistence); the red \
dashed curve is the spinodal (stability limit). Between them: \
metastable (nucleation and growth). Inside the spinodal: unstable \
(spinodal decomposition). The black dot is the critical point."],

hd2["2.3  Chain length and multivalency"],

para[{"A key prediction of Flory\[Dash]Huggins theory is that ",
   bold["longer polymers phase-separate more readily"], ". The critical \
interaction falls toward ", math["1/2"], " as chain length grows, with \
the excess ", math["\[Chi]\:1D9C - 1/2"], " scaling as ", math["N^(-1/2)"],
   ": a 1000-mer demixes at ", math["\[Chi]\:1D9C \[TildeTilde] 0.53"],
   " versus ", math["\[TildeTilde] 0.87"], " for a 10-mer, so long chains \
phase-separate at appreciably weaker interactions. This is biologically \
significant: intrinsically disordered \
regions (IDRs) of condensate-forming proteins are typically 50\[Dash]500 \
residues, and their effective ", math["N"], " is further increased by \
multivalent interactions (multiple binding sites per molecule, \
RNA scaffolds, post-translational modifications)."}],

wlInput["figFHChainLength[{1, 5, 20, 100, 500}]"],
captionCell["Figure 2.3  Spinodal curves for different chain lengths. \
Longer chains produce wider two-phase regions at higher temperatures \
(lower \[Chi] threshold). This is why proteins with long disordered \
regions are prone to phase separation."]
}];


(* ============================================================
   \[Section]3 Cahn\[Dash]Hilliard dynamics
   ============================================================ *)

writeAll[{

hd1["3.  Cahn\[Dash]Hilliard dynamics of droplet formation"],

hd2["3.1  The equation"],

para[{"The Flory\[Dash]Huggins theory tells us ", ital["whether"],
   " phase separation occurs, but not ", ital["how"], ". The dynamics \
of a conserved order parameter (composition must be conserved \[LongDash] \
matter doesn't appear or disappear) are governed by the ",
   bold["Cahn\[Dash]Hilliard equation"], ":"}],

displayMath[
   RowBox[{
     FractionBox[RowBox[{"\[PartialD]", "\[Phi]"}],
        RowBox[{"\[PartialD]", "t"}]],
     "=",
     RowBox[{"M", " ",
       SuperscriptBox["\[Del]", "2"],
       RowBox[{"(", RowBox[{
         FractionBox[RowBox[{"\[Delta]", "F"}],
            RowBox[{"\[Delta]", "\[Phi]"}]],
         "-",
         RowBox[{"\[Kappa]", " ", SuperscriptBox["\[Del]", "2"], "\[Phi]"}]
       }], ")"}]
     }]
   }]
],

para[{"where ", math["M"], " is the mobility (diffusion coefficient), \
", math["\[Delta]F/\[Delta]\[Phi]"], " is the chemical potential from \
the Flory\[Dash]Huggins free energy, and ", math["\[Kappa]"],
   " is the gradient energy coefficient that penalises sharp interfaces \
(it sets the surface tension \[Gamma] of the condensate)."}],

para[{ital["Physical content."], " The Cahn\[Dash]Hilliard equation \
says: material flows down chemical potential gradients (as in ordinary \
diffusion), but the chemical potential now includes a penalty for \
composition gradients (the \[Kappa] term). Inside the spinodal region, \
the bulk chemical potential drives ", bold["uphill diffusion"],
   " \[LongDash] material flows from dilute to concentrated regions, \
amplifying any fluctuation. The gradient penalty prevents infinitely \
sharp interfaces and selects a characteristic wavelength for the \
initial instability."}],

hd2["3.2  Spinodal decomposition"],

para[{"We solve the 1D Cahn\[Dash]Hilliard equation using Wolfram's \
built-in ", mono["NDSolve"], " with the Method of Lines. The Flory\
\[Dash]Huggins chemical potential is logarithmically stiff near the pure \
phases and the strongly-stabilised 1D solver freezes; we therefore use \
the standard symmetric double-well (Model B) order parameter ",
   math["\[Psi]"], ", started from a uniform state plus a small \
multi-mode periodic perturbation, and map it to a physical volume \
fraction ", math["\[Phi] = \[Phi]\:2080 + A\[ThinSpace]\[Psi]"],
   " for display."}],

codeIn["mu[\[Psi]_] := -a \[Psi] + b \[Psi]^3;   (* double-well Model B *)
sol = NDSolveValue[{
   D[\[Psi][x, t], t] == M D[
      mu[\[Psi][x, t]] - \[Kappa] D[\[Psi][x, t], x, x], x, x],
   \[Psi][x, 0] == icf[x],            (* smooth periodic perturbation *)
   \[Psi][0, t] == \[Psi][L, t],      (* periodic BCs *)
   \[Psi]'[0, t] == \[Psi]'[L, t]},
   \[Psi], {x, 0, L}, {t, 0, T},
   Method -> {\"MethodOfLines\",
     \"SpatialDiscretization\" -> {\"TensorProductGrid\",
       \"MinPoints\" -> 256}}];"],

wlInput["figCH1DEvolution[]"],
captionCell["Figure 3.1  Top: space\[Dash]time density plot of the 1D \
spinodal decomposition. The initially uniform composition rapidly \
develops concentration fluctuations that coarsen into well-separated \
dense (red) and dilute (blue) domains. Bottom: composition profiles at \
selected times showing the progression from small perturbation to \
fully phase-separated condensates."],

hd2["3.3  Coarsening: Lifshitz\[Dash]Slyozov\[Dash]Wagner"],

para[{"Once phase separation is complete, the system continues to \
evolve: large droplets grow at the expense of small ones. This is ",
   bold["Ostwald ripening"], ", described by the Lifshitz\[Dash]Slyozov\
\[Dash]Wagner (LSW) theory. The critical radius ",
   math["R*"], " is set by the Gibbs\[Dash]Thomson equation:"}],

displayMath[
   RowBox[{
     SubscriptBox["c", "eq"], "(", "R", ")", "=",
     SubscriptBox["c", "\[Infinity]"],
     RowBox[{"exp", "(", FractionBox[
       RowBox[{"2", "\[Gamma]", " ", SubscriptBox["V", "m"]}],
       RowBox[{"R", " ", "k", "T"}]], ")"}],
     "\[TildeTilde]",
     SubscriptBox["c", "\[Infinity]"],
     RowBox[{"(", RowBox[{"1", "+",
       FractionBox[SubscriptBox["l", "c"], "R"]}], ")"}]
   }]
],

para[{"Droplets smaller than the critical radius (", math["R < R*"],
   ") have an equilibrium concentration above the surrounding matrix \
supersaturation \[LongDash] they lose material and dissolve. Droplets with \
", math["R > R*"], " sit below it \[LongDash] they grow. The mean domain size grows as ", math["\[LeftAngleBracket]R\
\[RightAngleBracket] \[Proportional] t^(1/3)"],
   ", the classic LSW scaling."}],

wlInput["figProtoOstwald[]"],
captionCell["Figure 3.2  Ostwald ripening of a polydisperse condensate \
population (100 droplets). Left: the mean and maximum droplet radius \
grow with time as large droplets feed on small ones. Right: the number \
of surviving droplets falls as the smallest dissolve \[LongDash] the \
hallmark of LSW coarsening (\[LeftAngleBracket]R\[RightAngleBracket] \
\[Proportional] t^(1/3))."]
}];


(* ============================================================
   \[Section]4 Concentration enhancement and catalysis
   ============================================================ *)

writeAll[{

hd1["4.  Concentration enhancement and catalysis"],

hd2["4.1  The partition coefficient"],

para[{"When a condensate forms, molecules partition between the dense \
and dilute phases according to their chemical affinity. The ",
   bold["partition coefficient"], " ", math["K"], " is defined as:"}],

displayMath[
   RowBox[{
     "K", "=",
     FractionBox[SubscriptBox["c", "in"], SubscriptBox["c", "out"]]
   }]
],

para[{"For condensate-forming proteins, K is typically 10\[Dash]100. \
For RNA molecules that interact strongly with the condensate scaffold, \
K can reach 100\[Dash]1000. A molecule at 1 nM bulk concentration \
with K = 100 is at 100 nM inside the condensate. A molecule at 1 nM \
with K = 1000 is at 1 \[Mu]M. This is the concentration miracle."}],

hd2["4.2  Bimolecular rate enhancement"],

para[{"For a bimolecular reaction A + B \[Rule] C with second-order \
rate constant ", math["k"], ":"}],

displayMath[
   RowBox[{
     SubscriptBox["rate", "dilute"], "=",
     RowBox[{"k", " ",
       SubscriptBox["[A]", "out"], " ",
       SubscriptBox["[B]", "out"]}]
   }]
],

displayMath[
   RowBox[{
     SubscriptBox["rate", "condensate"], "=",
     RowBox[{"k", " ",
       RowBox[{SubscriptBox["K", "A"], " ",
         SubscriptBox["[A]", "out"]}], " ",
       RowBox[{SubscriptBox["K", "B"], " ",
         SubscriptBox["[B]", "out"]}]}],
     "=",
     RowBox[{SubscriptBox["K", "A"], " ",
       SubscriptBox["K", "B"], " ",
       SubscriptBox["rate", "dilute"]}]
   }]
],

para[{"If both reactants partition equally (",
   math["K_A \[TildeTilde] K_B \[TildeTilde] K"], "), the \
reaction rate is enhanced by ", bold["K\:00b2"], ". With \
K = 100, this is a ", bold["10\:2074-fold enhancement"],
   ". This is the ", ital["open-system"], " (buffered-reservoir) limit, \
and it is an upper bound: it assumes the surrounding dilute phase is never \
depleted. As \[Section]8.5 shows, imposing mass conservation caps the \
realised gain at order ", math["K/4"], ". Even so, the mechanism reframes \
the concentration problem in prebiotic chemistry: you do not need a \
concentrated primordial soup if you have phase separation."}],

wlInput["figRateEnhancement[]"],
captionCell["Figure 4.1  Reaction rate enhancement inside condensates \
as a function of the partition coefficient K, in the open-system limit \
(mass conservation lowers the realised gain \[LongDash] see \[Section]8.5). \
Bimolecular reactions (blue) scale as K\:00b2; the termolecular curve (red, \
K\:00b3) is shown only to illustrate the scaling, as three-body reactions \
are kinetically negligible in solution. At the typical condensate \
K \[TildeTilde] 100, the bimolecular bound is 10\:2074."],

hd2["4.3  Prebiotic chemistry phase space"],

para[{"Which combinations of partition coefficient and bulk \
concentration allow net chemical production? We map the net replication \
rate (replication minus degradation) across the K\[Dash]c parameter \
space:"}],

wlInput["figPrebioticSpace[1.0]"],
captionCell["Figure 4.2  Net replication rate as a function of partition \
coefficient K and bulk substrate concentration. The black contour \
marks the boundary between net growth (blue) and net decay (red). \
At low bulk concentrations (the likely prebiotic regime), only \
condensates with K > 50\[Dash]100 support net replication."]
}];


(* ============================================================
   \[Section]5 RNA world inside condensates
   ============================================================ *)

writeAll[{

hd1["5.  RNA world inside condensates"],

hd2["5.1  Ribozyme kinetics"],

para[{"The RNA world hypothesis proposes that the first self-replicating \
molecules were RNA ribozymes \[LongDash] RNA molecules that catalyse \
their own replication. The kinetics of ribozyme-catalysed replication \
follow ", bold["Michaelis\[Dash]Menten"], " saturation:"}],

displayMath[
   RowBox[{
     "v", "=",
     FractionBox[RowBox[{SubscriptBox["V", "max"], " ", "[S]"}],
        RowBox[{SubscriptBox["K", "m"], "+", "[S]"}]]
   }]
],

para[{"Inside a condensate, both the enzyme (ribozyme) and the substrate \
(template RNA) are concentrated by the partition coefficient K. Writing ",
   math["V\:2098\:2090\:2093"], " for the dilute-phase maximal rate (itself \
proportional to enzyme concentration), the in-condensate ",
   math["V\:2098\:2090\:2093"], " is raised by the factor K (more enzyme) \
and the substrate concentration [S] is likewise raised by K (more \
substrate):"}],

displayMath[
   RowBox[{
     SubscriptBox["v", "condensate"], "=",
     FractionBox[
       RowBox[{"K", " ", SubscriptBox["V", "max"], " ",
         "K", " ", "[S]"}],
       RowBox[{SubscriptBox["K", "m"], "+", "K", " ", "[S]"}]]
   }]
],

para[{"At low substrate concentrations (",
   math["K [S] \[LessLess] K\:2098"], "), the rate scales as ",
   math["K\:00b2"], " \[LongDash] the full bimolecular \
enhancement. At saturation (",
   math["K [S] \[GreaterGreater] K\:2098"], "), the rate scales as K \
(only the enzyme enhancement matters)."}],

hd2["5.2  Replication vs degradation"],

para[{"RNA is chemically unstable. In the prebiotic ocean, hydrolysis \
degrades RNA on timescales of hours to days. The ", bold["origin-of-life \
window"], " exists only when the replication rate exceeds the \
degradation rate:"}],

displayMath[
   RowBox[{
     SubscriptBox["v", "rep"], ">",
     SubscriptBox["v", "deg"],
     "   \[DoubleLongLeftRightArrow]   ",
     "K", ">", SuperscriptBox["K", "*"]
   }]
],

para[{"where ", math["K*"], " is the critical partition coefficient. \
Below K*, degradation wins and no self-replicating system can persist."}],

wlInput["figRibozyme[0.5]"],
captionCell["Figure 5.1  Replication (blue) versus degradation (red) rate \
as a function of the partition coefficient K, at fixed bulk substrate \
[S] = 1 \[Mu]M. Replication is enhanced as K\:00b2 and degradation as K, \
so above a critical K* the green-shaded origin-of-life window opens: \
replication outpaces hydrolysis and a self-replicating system can persist."],

hd2["5.3  The origin-of-life window"],

para[{"The critical partition coefficient K* depends on the ribozyme \
efficiency (V\:2098\:2090\:2093/K\:2098), the degradation rate constant, \
and the substrate level. For plausible prebiotic parameters \
(V\:2098\:2090\:2093 \[TildeTilde] 1 \[Mu]M/min, K\:2098 \[TildeTilde] \
10 \[Mu]M, k\:1D48\:2091\:2092 \[TildeTilde] 0.5/min), we find K* \
\[TildeTilde] 10 at [S] = 1 \[Mu]M, rising toward 50\[Dash]100 at the \
lower substrate levels of a dilute prebiotic pool. Modern coacervate \
droplets routinely achieve K = 50\[Dash]500, comfortably above this \
threshold \[LongDash] the concentration-enhancement mechanism is not \
marginal."}],

para[{"This is a quantitative version of Oparin's coacervate \
hypothesis: phase separation concentrates reactants enough to bootstrap \
biochemistry from a dilute ocean \[LongDash] though, as \[Section]8.5 \
shows, mass conservation bounds how large that boost can be."}]
}];


(* ============================================================
   \[Section]6 Protocell dynamics
   ============================================================ *)

writeAll[{

hd1["6.  Protocell dynamics"],

hd2["6.1  Growth and division"],

para[{"A condensate droplet that sustains net internal replication grows: \
material produced inside increases its volume. A ", ital["passive"],
   " liquid droplet, however, does not divide on its own \[LongDash] \
surface tension drives it to minimise area, so it stays spherical or \
coalesces with its neighbours. Genuine membrane-free division requires a \
", bold["non-equilibrium drive"], ": continuous chemical turnover that \
destabilises large droplets (Zwicker et al. 2017), or an external thermal \
gradient that pinches a droplet in two (", math["\[Section]9.4"],
   "; Ianeselli et al. 2022). In the toy protocell model of this section \
we therefore ", ital["impose"], " division as a simple size threshold \
\[LongDash] a stand-in for those mechanisms \[LongDash] and ask what \
population dynamics follow."}],

wlInput["figProtoGrowthDivision[]"],
captionCell["Figure 6.1  Schematic of the protocell life cycle. Small \
condensates grow by accretion of material from the bulk and by internal \
replication. When they exceed a critical size, an imposed threshold \
\[LongDash] standing in for the chemically- or thermally-driven fission \
of \[Section]9 \[LongDash] splits them into two daughters. Daughters \
inherit (with variation) \
the parent's composition, including its partition coefficient K and \
replication efficiency. Weaker daughters dissolve; fitter ones survive \
and grow."],

hd2["6.2  Proto-Darwinian selection"],

para[{"If daughter droplets inherit their parent's chemical composition \
(approximately), and if that composition affects their fitness \
(replication rate, stability), then we have the three ingredients of \
natural selection: ", bold["variation"], " (daughters differ slightly \
from each other), ", bold["inheritance"], " (daughters resemble parents \
more than random droplets), and ", bold["differential fitness"],
   " (some compositions grow faster than others). This is ",
   bold["proto-Darwinian selection"], " \[LongDash] natural selection \
operating on condensate populations before the invention of genetic \
material."}],

para[{"We simulate this by initialising a population of 50 protocells \
with random partition coefficients K and replication rates, then \
letting them grow, divide, compete for resources, and either survive \
or dissolve:"}],

wlInput["figProtoSelection[]"],
captionCell["Figure 6.2  Proto-Darwinian selection in a simulated \
protocell population. Left: mean partition coefficient K increases \
over time as protocells with higher K out-compete those with lower K. \
Centre: mean replication rate increases. Right: lineage depth (max \
generation number) grows as successful lineages proliferate. These \
trends are the signature of adaptive evolution operating on \
membrane-free liquid droplets."]
}];


(* ============================================================
   \[Section]7 Condensates and disease
   ============================================================ *)

writeAll[{

hd1["7.  Condensates and disease"],

hd2["7.1  Liquid to gel to aggregate"],

para[{"Functional condensates are ", bold["liquid"], ": they flow, \
they exchange material rapidly with the surrounding cytoplasm, and \
they dissolve when the cell no longer needs them. But this liquidity \
is not guaranteed. Over time, or in the presence of disease-associated \
mutations, condensates can undergo a ", bold["maturation"],
   " transition: liquid \[Rule] gel \[Rule] solid aggregate. The \
gel state has reduced dynamics; the solid aggregate is essentially \
irreversible."}],

para[{"This maturation is implicated in several ",
   bold["neurodegenerative diseases"], ", though the pathologies differ. \
FUS and TDP-43 form functional liquid condensates whose mutated or \
mislocalised forms harden into the cytoplasmic ", bold["inclusions"],
   " of ALS and frontotemporal dementia (FTLD); tau matures into the \
neurofibrillary ", bold["tangles"], " of Alzheimer's and other \
tauopathies. (The amyloid ", bold["plaques"], " of Alzheimer's, by \
contrast, are built from a different peptide, amyloid-\[Beta].) In each \
case disease mutations or conditions shift the phase boundaries, \
accelerating the liquid\[Rule]gel\[Rule]solid transition toward these \
pathological end states."}],

para[{"We model this as a three-state kinetic system:"}],

displayMath[
   RowBox[{
     "[liquid]",
     OverscriptBox["\[LongRightArrow]", SubscriptBox["k", "gel"]],
     "[gel]",
     OverscriptBox["\[LongRightArrow]", SubscriptBox["k", "agg"]],
     "[aggregate]"
   }]
],

wlInput["figDiseaseTransition[]"],
captionCell["Figure 7.1  Three kinetic scenarios for condensate ageing. \
Left: wild-type protein with slow maturation (healthy neuron). Centre: \
risk variant with moderate rates. Right: disease mutant (FUS/TDP-43 in \
ALS) with fast maturation \[LongDash] the liquid phase is rapidly \
depleted, and toxic aggregates accumulate."],

hd2["7.2  Phase diagram of material states"],

para[{"The condensate material state depends on protein concentration \
and temperature (or interaction strength). We construct a schematic \
phase diagram showing the four regimes: mixed (one phase), liquid \
condensate, gel, and solid aggregate. Disease mutations effectively \
shift the boundaries, expanding the gel and aggregate regions at the \
expense of the functional liquid state."}],

wlInput["figDiseasePhaseDiagram[]"],
captionCell["Figure 7.2  Condensate material state phase diagram. \
Higher concentration and lower temperature (stronger interactions) \
drive the progression from mixed solution to liquid condensate to gel \
to solid aggregate. Disease mutations shift the boundaries toward \
lower concentrations and higher temperatures, making pathological \
aggregation more likely."]
}];


(* ============================================================
   \[Section]8 Nonlinear dynamics
   ============================================================ *)

writeAll[{

hd1["8.  Nonlinear dynamics: stability, bifurcation, and coarsening"],

para[{"The previous sections built the model; here we analyse it with the \
toolkit of ", bold["nonlinear dynamics"], " \[LongDash] linear stability, \
bifurcation theory, Lyapunov functionals, and scaling laws. Every figure in \
this section is generated by ", mono["wolfram/nonlinear_dynamics.wls"],
   " (with a parallel Python implementation in ",
   mono["src/nonlinear_analysis.py"], ")."}],

hd2["8.1  Linear stability and the dispersion relation"],

para[{"Linearise the Cahn\[Dash]Hilliard equation about a uniform state ",
   math["\[Phi] = \[Phi]\:2080 + \[Delta]\[ThinSpace]exp(i\[ThinSpace]k\
\[CenterDot]x + \[Omega]t)"], ". Each Fourier mode grows or decays \
independently at the rate"}],

displayMath[
   RowBox[{"\[Omega]", "(", "k", ")", "=", RowBox[{"-", "M", " ",
     SuperscriptBox["k", "2"], " ",
     RowBox[{"(", RowBox[{
        SuperscriptBox["f", "\[Prime]\[Prime]"], "(",
        SubscriptBox["\[Phi]", "0"], ")", "+", "\[Kappa]", " ",
        SuperscriptBox["k", "2"]}], ")"}]}]}]
],

para[{"A uniform state is ", bold["unstable"], " precisely when ",
   math["f\[Prime]\[Prime](\[Phi]\:2080) < 0"], " \[LongDash] the ",
   bold["spinodal"], " condition \[LongDash] in which case modes with ",
   math["0 < k < k_c = sqrt(-f\[Prime]\[Prime]/\[Kappa])"],
   " grow. The fastest-growing mode and its rate are ",
   math["k* = sqrt(-f\[Prime]\[Prime]/2\[Kappa])"], " and ",
   math["\[Omega]* = M (f\[Prime]\[Prime])\:00b2/4\[Kappa]"],
   ", setting the initial droplet spacing ", math["\[Lambda]* = 2\[Pi]/k*"],
   ". This single mode is the seed of all condensate morphology."}],

tryNote[{"The dispersion relation, live \[LongDash] the positive lobe (and \
hence the fastest-growing mode) appears exactly when ", math["f\[Prime]\
\[Prime] < 0"], ":"}],

wlInput["\[Omega][k_, fpp_, M_, \[Kappa]_] := -M k^2 (fpp + \[Kappa] k^2);\nPlot[Evaluate@Table[\[Omega][k, fpp, 1, 1], {fpp, {-2., -1., -0.5, 0.5}}],\n  {k, 0, 1.6}, Frame -> True, FrameLabel -> {\"k\", \"\[Omega](k)\"},\n  PlotLegends -> {\"f''=-2\", \"f''=-1\", \"f''=-0.5\", \"f''=+0.5 (stable)\"}]"],

wlInput["figDispersion[]"],
captionCell["Figure 8.1  Cahn\[Dash]Hilliard dispersion relation for the \
Flory\[Dash]Huggins free energy (\[Chi] = 1.5). Curves with a positive lobe \
are unstable; dots mark the fastest-growing mode k*. As the mean composition \
\[Phi]\:2080 approaches the spinodal edge, the growth rate collapses to zero \
and the selected wavelength diverges."],

hd2["8.2  Mode selection across the parameter plane"],

para[{"The selected wavelength ", math["\[Lambda]* = 2\[Pi] sqrt(2\[Kappa]/\
(-f\[Prime]\[Prime]))"], " grows with the surface-tension parameter ",
   math["\[Kappa]"], " and diverges at the spinodal boundary, where ",
   math["f\[Prime]\[Prime] \[Rule] 0"], ". The map below shows how the \
length scale is selected jointly by composition and ", math["\[Kappa]"], "."}],

wlInput["figModeSelection[]"],
captionCell["Figure 8.2  Fastest-growing wavelength \[Lambda]* over the \
(\[Phi]\:2080, \[Kappa]) plane (clipped at 25 for contrast). The white curve \
is the spinodal f\[Prime]\[Prime] = 0; outside it the homogeneous state is \
linearly stable and no length scale is selected."],

hd2["8.3  Bifurcation: common-tangent coexistence vs. the spinodal"],

para[{"Sweeping the control parameter traces out the ", bold["bifurcation \
diagram"], ". For the symmetric double well ", math["f = -a/2\[ThinSpace]\
\[Phi]\:00b2 + b/4\[ThinSpace]\[Phi]\:2074"], " the coexistence amplitude \
follows a supercritical pitchfork: the binodal ", math["\[Phi] = \
\[PlusMinus]sqrt(a/b)"], " encloses the spinodal ",
   math["\[Phi] = \[PlusMinus]sqrt(a/3b)"],
   ". For Flory\[Dash]Huggins the coexistence curve is the ",
   bold["common-tangent"], " (Maxwell) construction. Because the osmotic \
pressure is linear in \[Chi], the dense branch follows the closed form ",
   math["\[Chi](\[Phi]_R) = (\[Phi]_R g(\[Phi]_R) - h(\[Phi]_R))/\[Phi]_R\:00b2"],
   ", which is numerically robust where the repository's earlier \
two-variable solver failed."}],

wlInput["figNLPhaseDiagram[]"],
captionCell["Figure 8.3  (a) Flory\[Dash]Huggins phase diagram: binodal \
(coexistence, blue) outside the spinodal (instability, red dashed), meeting \
at the critical point; the orange dot marks the solver operating point \
(\[Phi] = 0.3, \[Chi] = 1.5). (b) The symmetric double-well pitchfork, with \
the metastable region between binodal and spinodal."],

hd2["8.4  The free energy as a Lyapunov functional"],

para[{"Cahn\[Dash]Hilliard dynamics are a ", bold["gradient flow"],
   " of the free-energy functional ", math["F[\[Phi]] = \[Integral] [ f(\[Phi]\
) + (\[Kappa]/2) |\[Del]\[Phi]|\:00b2 ] dx"], ". Conservation of \[Phi] and \
the structure of the equation guarantee ", math["dF/dt \[LessEqual] 0"],
   ": F is a ", bold["Lyapunov functional"],
   " and the phase-separated state is its minimiser. Each coarsening event \
\[LongDash] two domains merging \[LongDash] appears as a step down in F."}],

wlInput["figEnergyLyapunov[]"],
captionCell["Figure 8.4  (a) Space\[Dash]time evolution of the 1D \
double-well Cahn\[Dash]Hilliard equation (NDSolve): domains form, then \
coarsen by merging. (b) The free energy F[\[Phi]] decreases monotonically, \
falling in discrete steps as domains coalesce \[LongDash] the signature of a \
Lyapunov functional. (Note: the strongly stabilised Flory\[Dash]Huggins \
solver instead freezes; the canonical t\[ThinSpace]^(1/3) coarsening, with ",
   math["L\:00b3 = L\:2080\:00b3 + K\[ThinSpace]t"], " linear in time, is \
reproduced by the double-well Model B in the Python implementation.)"],

hd2["8.5  Mass-conserving reaction enhancement"],

para[{"Section 4 estimated the bimolecular rate boost inside a condensate as \
", math["K_A K_B"], ". That is the ", ital["open-system"], " limit \
(buffered reservoir). Imposing ", bold["mass conservation"], " \[LongDash] \
the condensate concentrates reactants but ", ital["depletes"],
   " the surrounding dilute phase \[LongDash] gives a very different answer. \
With dense area fraction ", math["p"], " and rapid partitioning ",
   math["c_in = K\[ThinSpace]c_out"], " in each phase, the volume-averaged \
rate relative to a well-mixed system is"}],

displayMath[
   RowBox[{FractionBox[SubscriptBox["k", "eff"], SubscriptBox["k", "out"]],
     "=", FractionBox[
       RowBox[{"p", " ", SubscriptBox["K", "A"], SubscriptBox["K", "B"],
         "+", RowBox[{"(", RowBox[{"1", "-", "p"}], ")"}]}],
       RowBox[{RowBox[{"(", RowBox[{"1", "+",
            RowBox[{"p", RowBox[{"(", RowBox[{SubscriptBox["K", "A"], "-",
               "1"}], ")"}]}]}], ")"}],
         RowBox[{"(", RowBox[{"1", "+",
            RowBox[{"p", RowBox[{"(", RowBox[{SubscriptBox["K", "B"], "-",
               "1"}], ")"}]}]}], ")"}]}]]}]
],

para[{"The enhancement vanishes at both ", math["p \[Rule] 0"], " (no \
condensate) and ", math["p \[Rule] 1"], " (everything uniform); it is \
maximal at an ", bold["intermediate dense fraction"], " ", math["p* \
\[TildeTilde] 1/K"], ", where the gain is only ", math["\[TildeTilde] K/4"],
   " \[LongDash] not the naive ", math["K\:00b2"], ". Bulk depletion is the \
hard ceiling on condensate catalysis, and it predicts an optimal droplet \
volume fraction \[LongDash] a falsifiable, quantitative claim."}],

wlInput["figReactionEnhancement[]"],
captionCell["Figure 8.5  Mass-conserving rate enhancement ",
   math["k_eff/k_out"], ". (a) Dependence on the partition coefficients at \
fixed dense fraction p = 0.05. (b) Dependence on dense fraction for several \
K: each curve peaks at p* \[TildeTilde] 1/K with a maximum gain of order \
K/4, far below the open-system K\:00b2."]
}];


(* ============================================================
   \[Section]9 Active and chemically-maintained condensates
   ============================================================ *)

writeAll[{

hd1["9.  Active and chemically-maintained condensates"],

para[{"Living condensates are not at equilibrium: ATP-driven reactions \
continually build up and break down their components. This section extends \
the model with chemistry. Code: ", mono["src/active_models.py"], ", ",
   mono["src/fh_extensions.py"], ", and ", mono["wolfram/active_dynamics.wls"], "."}],

hd2["9.1  Reaction-arrested coarsening"],

para[{"Add a linear chemical reaction that interconverts droplet material and \
a soluble form, ", math["s(\[Phi]) = k_+(1-\[Phi]) - k_-\[Phi] = -k_react\
(\[Phi]-\[Phi]_ss)"], ", to the conserved dynamics:"}],

displayMath[
   RowBox[{FractionBox["\[PartialD]\[Phi]", "\[PartialD]t"], "=",
     "M", " ", SuperscriptBox["\[Del]", "2"], "\[Mu]", "-",
     SubscriptBox["k", "react"],
     RowBox[{"(", RowBox[{"\[Phi]", "-", SubscriptBox["\[Phi]", "ss"]}], ")"}]}]
],

para[{"The reaction adds a uniform ", math["-k_react"], " to the growth \
rate, ", math["\[Omega](k) = -Mk\:00b2(f\[Prime]\[Prime]+\[Kappa]k\:00b2) - k_react"],
   ", which ", bold["stabilises the long-wavelength modes"], ": the unstable \
band becomes finite and a length scale is selected. Ostwald ripening is \
arrested and the system settles into a steady-state emulsion \[LongDash] the \
mechanism behind chemically-active droplets (Zwicker et al. 2015; \
Bauermann et al. 2025)."}],

wlInput["figActiveDispersion[]"],
captionCell["Figure 9.1  Reaction-modified dispersion relation. Increasing \
the reaction rate pushes the small-k (long-wavelength) modes below zero, \
leaving a finite unstable band \[LongDash] coarsening can no longer proceed \
to arbitrarily large domains."],

wlInput["figReaction1D[]"],
captionCell["Figure 9.2  1D reaction\[Dash]Cahn\[Dash]Hilliard (NDSolve). \
(a) Without reaction the domains coarsen toward a single domain. (b) With \
reaction the pattern freezes into a steady-state periodic structure of \
selected wavelength."],

para[{"The full 2D problem (solved with a spectral method in ",
   mono["src/active_models.py"], ") sharpens the same picture: chemical \
turnover turns the coarse, ever-growing morphology into a fine ",
   bold["steady-state emulsion"], " whose mean droplet size shrinks with \
reaction rate, and whose size distribution is narrower than the passive LSW \
one (Bauermann et al. 2025)."}],

hd2["9.2  Active Model B: broken detailed balance"],

para[{"A genuinely non-equilibrium term ", math["\[Lambda]|\[Del]\[Phi]|\:00b2"],
   " added to the chemical potential breaks detailed balance (Active Model B; \
Tjhung et al. 2018). Measured honestly, this term shifts the coexistence densities \
", bold["asymmetrically"], " (the dense phase becomes denser than the dilute \
phase is dilute) and modestly slows coarsening, but does ", ital["not"], " by \
itself arrest it \[LongDash] the robust route to a steady-state emulsion is \
the chemical reaction of \[Section]9.1."}],

hd2["9.3  Stickers and spacers: multivalency drives phase separation"],

para[{"Real biomolecules interact through specific, saturable contacts \
(\"stickers\" joined by flexible \"spacers\"; Harmon & Pappu 2017). Adding a \
saturable association free energy ", math["f_bond = \[Phi]_s [ln(1-p)+\
p/2]"], " (bond fraction ", math["p"], " from the law of mass action ",
   math["p/(1-p)\:00b2 = K_a \[Phi]_s"], ", ", math["\[Phi]_s = v\[Phi]"],
   ") to the Flory\[Dash]Huggins free energy expands the two-phase region and lowers \
the critical interaction: higher sticker valence ", math["v"], " makes phase \
separation far easier."}],

wlInput["figStickers[]"],
captionCell["Figure 9.3  Stickers-and-spacers spinodal curves. Each curve \
bounds the two-phase region; increasing sticker valence pushes the critical \
point down and outward \[LongDash] multivalency drives liquid\[Dash]liquid \
phase separation \
even when simple mixing is unfavourable."],

hd2["9.4  Thermal-gradient fission: membrane-free division"],

para[{"A temperature gradient makes the interaction ", math["\[Chi] ~ 1/T"],
   " position-dependent. A warm band (weakened cohesion) drawn through an \
elongated droplet dissolves its neck and splits it into two daughters \
\[LongDash] division without a membrane, a candidate protocell mechanism \
(Ianeselli 2022). The split is genuinely band-driven: two domains form with \
the band, only one without it."}],

wlInput["figThermalFission[]"],
captionCell["Figure 9.4  Thermal-gradient fission (1D Cahn\[Dash]Hilliard, \
NDSolve). (a) At uniform temperature a single domain persists. (b) A central \
warm band, where the well depth a(x) dips below zero, dissolves the neck and \
splits the domain into two daughters \[LongDash] membrane-free division."],

hd2["9.5  The dispersion surface"],

para[{"Finally, the linear-stability result of \[Section]8 viewed as a \
surface: the growth rate ", math["\[Omega](k,\[Phi]_0)"], " over the whole \
plane of wavenumber and composition, with the neutral-stability ridge ",
   math["\[Omega]=0"], " marking the boundary of the unstable region."}],

wlInput["figDispersionSurface[]"],
captionCell["Figure 9.5  The Cahn\[Dash]Hilliard dispersion surface \[Omega](k, \
\[Phi]\:2080) for the Flory\[Dash]Huggins free energy. The unstable dome shrinks to \
nothing as the composition approaches the spinodal edge; the cyan curve is \
the neutral-stability ridge."]
}];


(* ============================================================
   \[Section]10 Phase 3: selection, patterns and cycles
   ============================================================ *)

writeAll[{

hd1["10.  Emergent complexity: selection, patterns and cycles"],

para[{"The final layer of the model lets genuinely life-like behaviour \
emerge from the same physics. Code: ", mono["src/multicomponent.py"], ", ",
   mono["src/darwinian_selection.py"], ", ",
   mono["src/reaction_diffusion_cargo.py"], ", ",
   mono["src/wet_dry_cycling.py"], ", ", mono["src/parameter_space.py"],
   ", and the Wolfram companion ", mono["wolfram/phase3_extensions.wls"], "."}],

hd2["10.1  Multi-component condensates"],

para[{"Real condensates contain dozens of macromolecules. For ", math["N"],
   " mutually attracting components that co-condense in near-stoichiometric \
proportions, the demixing entropy is reduced as if they were one effective \
polymer of length ", math["N"], ", and the symmetric Flory\[Dash]Huggins critical \
point migrates (Qian & Knowles 2022):"}],

displayMath[
   RowBox[{SubscriptBox["\[Phi]", "c"], "=",
     FractionBox["1", RowBox[{"1", "+", SqrtBox["N"]}]], ",",
     "\[InvisibleSpace]", "\[InvisibleSpace]",
     SubscriptBox["\[Chi]", "c"], "=",
     FractionBox[RowBox[{"2", " ",
        SuperscriptBox[RowBox[{"(", "1", "+", SqrtBox["N"], ")"}], "2"]}],
       "N"]}]
],

para[{"As ", math["N"], " grows the critical concentration falls toward zero \
\[LongDash] multivalent mixtures phase-separate at vanishing concentration \
\[LongDash] and the critical interaction drops, so they demix far more easily. \
Distinct scaffolds give distinct condensates, and each client partitions \
selectively into the one it binds best (the physical basis of organelle \
identity)."}],

wlInput["figMulticomponent[]"],
captionCell["Figure 10.1  Multi-component spinodals for ", math["N=2,3,5,10"],
   " in the Qian\[Dash]Knowles convention; the critical point (white dots) \
migrates as ", math["\[Phi]_c=1/(1+\[Sqrt]N)"], " toward low concentration \
and low interaction as the number of components grows."],

wlInput["figCompartmentalization[]"],
captionCell["Figure 10.2  Selective compartmentalisation: a ternary \
(A, B, solvent) free energy in which two mutually-repelling solutes form an \
A-rich and a B-rich condensate \[LongDash] each client partitions into the \
droplet it binds best, the physical basis of organelle identity."],

hd2["10.2  Darwinian selection among protocells"],

para[{"The intellectual climax. Replicators partition into condensate \
droplets and replicate at a rate set by a heritable genotype ", math["g"],
   ", competing for a ", bold["shared global nutrient"],
   " under a uniform chemostat dilution. Replicators whose growth rate ",
   math["r_0 g"], " falls below the dilution rate are out-competed \
(competitive exclusion), and the population's mean fitness climbs. A full 2D \
spatial version (", mono["src/darwinian_selection.py"], "), in which \
fuel-driven droplets grow, divide under periodic thermal pulses, and pass \
their genotype to daughters, reaches the same outcome with explicit lineages."}],

para[{bold["The result is natural selection from pure physics."],
   " The cargo-weighted mean fitness rises over generations and saturates at \
a mutation-selection balance \[LongDash] no biology is hand-coded; selection \
is an emergent consequence of the dynamics."}],

para[{ital["What the selection requires."], " The effect is not automatic: it \
demands ", bold["Malthusian (per-capita) replication competing for a shared, \
global nutrient"], " together with a uniform chemostat dilution. Replacing the \
shared resource with a ", ital["local"], " carrying capacity, or raising the \
mutation rate too high, erases the inheritance signal \[LongDash] mean fitness \
then drifts rather than climbs and selection collapses. The trajectory below is \
therefore a genuine result of the chosen ecology, not a foregone conclusion."}],

wlInput["figSelectionChemostat[]"],
captionCell["Figure 10.3  Chemostat competition (", mono["NDSolve"], "). \
Left: the population's mean fitness ", math["\[LeftAngleBracket]g\
\[RightAngleBracket]"], " rises toward ", math["g_max"], " as fitter \
replicators take over. Right: the fittest genotype out-competes the rest \
\[LongDash] competitive exclusion, with no biology hand-coded."],

wlInput["figLineage[]"],
captionCell["Figure 10.4  A protocell phylogeny generated by the same \
selection dynamics: fit (yellow) lineages branch and proliferate while \
low-fitness (blue) lineages die out (node colour = fitness)."],

hd2["10.3  Reaction-diffusion inside condensates"],

para[{"A condensate is also a chemical reactor. When crowding slows the \
activator far more than the inhibitor, the diffusion ratio ", math["D_v/D_u"],
   " crosses the ", bold["Turing threshold"], " ", ital["only inside"],
   " the dense phase: a band of wavenumbers becomes unstable there while every \
mode stays stable in the dilute surroundings, so spots self-organise inside \
the condensate and nowhere else. (The full 2D patterns, and the reaction \
hotspot that forms in a shell at the droplet interface, are computed in ",
   mono["src/reaction_diffusion_cargo.py"], ".)"}],

wlInput["figTuring[]"],
captionCell["Figure 10.5  Schnakenberg linear-stability growth rate ",
   math["\[Lambda](k)"], ": inside the dense phase a band of wavenumbers is \
unstable (patterns form), while in the dilute phase every mode is stable. \
The Turing instability is switched on by the condensate."],

hd2["10.4  Wet-dry cycling"],

para[{"On the early Earth, rock pools cycled between flooded (wet) and \
evaporated (dry) states. Cycling the interaction ", math["\[Chi](t)=\[Chi]_0+\
\[CapitalDelta]\[Chi] sin(2\[Pi]t/T)"], " swings the system between a \
one-phase well (dissolve, remix) and a spinodally unstable double well \
(re-nucleate many droplets). Each dissolution erases the size differences that \
drive Ostwald ripening, so a ", bold["population"],
   " of droplets is maintained, while a matched constant-", math["\[Chi]"],
   " run coarsens toward a single winner (Ianeselli 2022)."}],

wlInput["figWetDry[]"],
captionCell["Figure 10.6  1D wet\[Dash]dry Cahn\[Dash]Hilliard space\[Dash]time \
(", mono["NDSolve"], "): with a cyclic well depth ", math["a(t)"],
   ", domains dissolve at each wet half-cycle and re-nucleate at each dry \
half-cycle. The repeated dissolution resets the size differences that would \
otherwise drive Ostwald ripening, so a population of droplets is sustained \
rather than coarsening to a single winner (Ianeselli et al. 2022)."],

hd2["10.5  The accessible parameter space"],

para[{"Finally, two analytic maps of the model's regimes: the ",
   math["(\[Chi],N)"], " plane separating one-phase from two-phase by the \
critical line ", math["\[Chi]_c(N)"], ", and the ", math["(K,p)"],
   " mass-conserving rate-enhancement landscape with its optimal ridge ",
   math["p*\[TildeTilde]1/K"], ". (The full four-panel regime atlas \[LongDash] \
which also includes the thermodynamic ", math["(\[Phi]_0,\[Kappa])"],
   " plane and the active dissolution / microphase / macrophase plane from 2D \
simulations \[LongDash] is built in ", mono["src/parameter_space.py"], ".)"}],

wlInput["figParamSpace[]"],
captionCell["Figure 10.7  The two analytic regime planes. Left, ",
   math["(\[Chi],N)"], ": the one-phase / two-phase boundary along the critical \
line. Right, ", math["(K,p)"], ": the mass-conserving rate enhancement, with \
the optimal ridge ", math["p*\[TildeTilde]1/K"], " (dashed)."]
}];


(* ============================================================
   \[Section]11 Conclusions
   ============================================================ *)

writeAll[{

hd1["11.  Conclusions"],

para[{ital["What we showed."], " (1) The Flory\[Dash]Huggins lattice \
theory, implemented in a few lines of Wolfram Language, correctly \
predicts the qualitative phase behaviour of biomolecular condensates: \
a UCST phase diagram with the critical concentration shifting to lower \
values for longer/more multivalent molecules. (2) The Cahn\[Dash]Hilliard \
equation, solved natively by ", mono["NDSolve"], ", captures spinodal \
decomposition; the subsequent Lifshitz\[Dash]Slyozov t\[ThinSpace]^(1/3) \
coarsening is reproduced with the double-well Model B form, since the \
stiff Flory\[Dash]Huggins solver freezes near the pure phases. \
(3) The concentration-enhancement argument gives a K\:00b2 boost for \
bimolecular reactions in the open-system limit \[LongDash] a partial \
resolution of the dilute-ocean problem, since mass conservation caps the \
realised gain at order K/4 (\[Section]8.5). (4) Protocell simulations \
show that condensate populations can undergo proto-Darwinian selection \
without genetic material. (5) The liquid\[Dash]to\[Dash]gel\[Dash]to\
\[Dash]aggregate transition connects condensate physics to \
neurodegenerative disease. (6) Layering the model \[LongDash] multi-component \
mixtures, a replicating cargo with heredity and a shared nutrient, \
reaction-diffusion, and wet-dry cycling \[LongDash] makes genuinely life-like \
behaviour emerge: ", bold["natural selection, Turing morphogenesis, and \
ripening-resistant droplet populations all fall out of the same physics."]}],

para[{ital["The bigger picture."], " Biomolecular condensates sit at \
the intersection of soft-matter physics, molecular biology, the origin \
of life, and neurodegeneration. The mathematical framework \[LongDash] \
Flory\[Dash]Huggins thermodynamics, Cahn\[Dash]Hilliard dynamics, \
Michaelis\[Dash]Menten kinetics, Lifshitz\[Dash]Slyozov coarsening \
\[LongDash] is classical and well-understood. What is new is the \
recognition that this physics operates inside living cells, that it \
provides the thermodynamic foundation for protocellular compartments, \
and that its failure modes are disease. The condensate revolution is \
not a discovery of new physics; it is the discovery that old physics \
was hiding in plain sight inside biology."}],

para[{bold["Key takeaways."]}],
Cell[TextData[{bold["Phase separation is thermodynamics, not biology. "],
   "A two-parameter Flory\[Dash]Huggins free energy already predicts when \
proteins and RNA demix into condensates."}], "Item"],
Cell[TextData[{bold["Concentration is the payoff \[LongDash] but it is \
bounded. "], "Partitioning boosts bimolecular rates by up to K\:00b2 in the \
open-system limit, while mass conservation caps the realised gain at order \
K/4."}], "Item"],
Cell[TextData[{bold["Life-like behaviour is emergent. "], "Add chemistry, \
heredity and a shared resource and natural selection, Turing patterns and \
ripening-resistant populations fall out of the same equations \[LongDash] \
no biology hand-coded."}], "Item"],
Cell[TextData[{bold["The same physics, run backwards, is disease. "],
   "The liquid\[Rule]gel\[Rule]solid maturation that endangers neurons is \
the flip side of the transition that organises healthy cells."}], "Item"]
}];


(* ============================================================
   \[Section]12 References
   ============================================================ *)

writeAll[{

hd1["12.  References"],

para[{bold["Qian, D., Michaels, T. C. T., & Knowles, T. P. J."], " (2022). \
Analytical solution to the Flory\[Dash]Huggins model for multicomponent phase \
separation. ", ital["Journal of Physical Chemistry Letters"], " 13, \
7853\[Dash]7860. [Multi-component critical point.]"}],

para[{bold["Turing, A. M."], " (1952). The chemical basis of morphogenesis. ",
   ital["Philosophical Transactions of the Royal Society B"], " 237, \
37\[Dash]72. [Reaction-diffusion patterns.]"}],

para[{bold["Eigen, M."], " (1971). Selforganization of matter and the \
evolution of biological macromolecules. ", ital["Naturwissenschaften"],
   " 58, 465\[Dash]523. [Competitive exclusion of replicators.]"}],

para[{bold["Ianeselli, A., Atienza, M., Kudella, P. W., et al."], " (2022). \
Water cycles dissolve and divide protocell-like droplets. ",
   ital["Nature Chemistry"], " 14, 32\[Dash]39. [Wet-dry cycling and \
thermal fission.]"}],

para[{bold["Flory, P. J."], " (1942). Thermodynamics of high polymer \
solutions. ", ital["Journal of Chemical Physics"], " 10, 51\[Dash]61."}],

para[{bold["Huggins, M. L."], " (1942). Some properties of solutions \
of long-chain compounds. ", ital["Journal of Physical Chemistry"],
   " 46, 151\[Dash]158."}],

para[{bold["Cahn, J. W., & Hilliard, J. E."], " (1958). Free energy \
of a nonuniform system. I. Interfacial free energy. ", ital["Journal \
of Chemical Physics"], " 28, 258\[Dash]267."}],

para[{bold["Lifshitz, I. M., & Slyozov, V. V."], " (1961). The \
kinetics of precipitation from supersaturated solid solutions. ",
   ital["Journal of Physics and Chemistry of Solids"], " 19, \
35\[Dash]50."}],

para[{bold["Oparin, A. I."], " (1924). ", ital["Proiskhozhdenie Zhizni"],
   " (The Origin of Life). Moscow: Moskovskii Rabochii."}],

para[{bold["Brangwynne, C. P., Eckmann, C. R., Courson, D. S., \
Rybarska, A., Hoege, C., Gharakhani, J., J\[UDoubleDot]licher, F., \
& Hyman, A. A."], " (2009). Germline P granules are liquid droplets \
that localize by controlled dissolution/condensation. ", ital["Science"],
   " 324, 1729\[Dash]1732. [The founding observation.]"}],

para[{bold["Hyman, A. A., Weber, C. A., & J\[UDoubleDot]licher, F."],
   " (2014). Liquid\[Dash]liquid phase separation in biology. ",
   ital["Annual Review of Cell and Developmental Biology"], " 30, \
39\[Dash]58."}],

para[{bold["Banani, S. F., Lee, H. O., Hyman, A. A., & Rosen, M. K."],
   " (2017). Biomolecular condensates: organisers of cellular biochemistry. ",
   ital["Nature Reviews Molecular Cell Biology"], " 18, 285\[Dash]298. \
[Comprehensive review.]"}],

para[{bold["Shin, Y., & Brangwynne, C. P."], " (2017). Liquid phase \
condensation in cell physiology and disease. ", ital["Science"], " 357, \
eaaf4382."}],

para[{bold["Patel, A., Lee, H. O., Jawerth, L., et al."], " (2015). \
A liquid-to-solid phase transition of the ALS protein FUS accelerated \
by disease mutation. ", ital["Cell"], " 162, 1066\[Dash]1077. [Disease \
connection.]"}],

para[{bold["Zwicker, D., Seyboldt, R., Weber, C. A., Hyman, A. A., \
& J\[UDoubleDot]licher, F."], " (2017). Growth and division of active \
droplets provides a model for protocells. ", ital["Nature Physics"],
   " 13, 408\[Dash]413. [Protocell dynamics.]"}],

para[{bold["Dignon, G. L., Zheng, W., Best, R. B., Kim, Y. C., \
& Mittal, J."], " (2018). Relation between single-molecule properties \
and phase behavior of intrinsically disordered proteins. ",
   ital["Proceedings of the National Academy of Sciences"], " 115, \
9929\[Dash]9934."}],

para[{bold["Wagner, C."], " (1961). Theorie der Alterung von \
Niederschl\[ADoubleDot]gen durch Uml\[ODoubleDot]sen (Ostwald-Reifung). ",
   ital["Zeitschrift f\[UDoubleDot]r Elektrochemie"], " 65, \
581\[Dash]591. [The \"W\" in LSW coarsening.]"}],

para[{bold["Zwicker, D., Hyman, A. A., & J\[UDoubleDot]licher, F."],
   " (2015). Suppression of Ostwald ripening in active emulsions. ",
   ital["Physical Review E"], " 92, 012317. [Chemically active droplets.]"}],

para[{bold["Tjhung, E., Nardini, C., & Cates, M. E."], " (2018). Cluster \
phases and bubbly phase separation in active fluids: reversal of the \
Ostwald process. ", ital["Physical Review X"], " 8, 031080. \
[Active Model B.]"}],

para[{bold["Harmon, T. S., Holehouse, A. S., Rosen, M. K., & Pappu, R. V."],
   " (2017). Intrinsically disordered linkers determine the interplay \
between phase separation and gelation in multivalent proteins. ",
   ital["eLife"], " 6, e30294. [Stickers and spacers.]"}],

para[{bold["Schnakenberg, J."], " (1979). Simple chemical reaction systems \
with limit cycle behaviour. ", ital["Journal of Theoretical Biology"],
   " 81, 389\[Dash]400. [Reaction\[Dash]diffusion model.]"}],

para[{bold["Bauermann, J., Bartolucci, G., Weber, C. A., & \
J\[UDoubleDot]licher, F."], " (2025). The droplet size distribution and \
its dynamics in chemically active emulsions. ",
   ital["Physical Review Letters"], " 135, 148201. [Active-emulsion size \
distributions.]"}]
}];


(* ============================================================
   Further reading (unnumbered, RiverNetworkStatistics style)
   ============================================================ *)

writeAll[{

hd1["Further reading"],

para[{"For an accessible, non-technical tour of the same story \[Dash] \
condensates, coacervates and the idea that the origin of life may have been \
hiding inside our own cells all along \[Dash] see the ", ital["New Scientist"],
   " cover feature that prompted this notebook:"}],

para[{bold["Marshall, M."], " (2026). Hiding in plain sight. ",
   ital["New Scientist"], ", 20 June 2026, pp. 30\[Dash]33 (Features / Cover \
story). Online edition titled ", ital["\[OpenCurlyDoubleQuote]Has the answer \
to life\[CloseCurlyQuote]s origins been hiding in our cells all along?\
\[CloseCurlyDoubleQuote]"], " ",
   link["newscientist.com/article/2529162",
     "https://www.newscientist.com/article/2529162-has-the-answer-to-lifes-origins-been-hiding-in-our-cells-all-along/"],
   "."}]
}];


(* ============================================================
   \[Section]13 Reproducibility
   ============================================================ *)

writeAll[{

hd1["13.  Reproducibility"],

para[{"All source code lives in the repository ",
   link["github.com/mthiel74/biomolecular-condensates",
     "https://github.com/mthiel74/biomolecular-condensates"],
   ". The pipeline has two parallel tracks: Wolfram Language (analysis \
scripts in ", mono["wolfram/"], " and this notebook builder in ",
   mono["community/"], ") and Python (", mono["src/"], ", including \
a full 2D Cahn\[Dash]Hilliard spectral solver). The figures come from both \
tracks: the Wolfram scripts write directly to ", mono["docs/images/"],
   ", while the Python scripts write to ", mono["figures/"], " and are \
copied across. Wolfram Language 13+ is required. To regenerate this \
notebook from scratch, run the Python track first, then the Wolfram \
runner (which now also executes the Phase-3 companion and copies the \
Python figures into place), then the builder:"}],

codeIn["# 1. Python track \[LongRightArrow] writes figures to figures/
pip install -r src/requirements.txt    # numpy, scipy, matplotlib (+ ffmpeg for animations)
python3 src/flory_huggins.py
python3 src/cahn_hilliard.py
python3 src/condensate_kinetics.py
python3 src/protocell.py
python3 src/nonlinear_analysis.py        # nonlinear-dynamics analysis (\[Section]8)
python3 src/fh_extensions.py             # stickers-and-spacers (\[Section]9.3)
python3 src/active_models.py             # active condensates (\[Section]9)
python3 src/multicomponent.py            # multi-component condensates (\[Section]10.1)
python3 src/darwinian_selection.py       # Darwinian selection (\[Section]10.2)
python3 src/reaction_diffusion_cargo.py  # reaction-diffusion (\[Section]10.3)
python3 src/wet_dry_cycling.py           # wet-dry cycling (\[Section]10.4)
python3 src/parameter_space.py           # parameter space (\[Section]10.5)
python3 src/phase3_heroes.py             # hero visualisations

# 2. Wolfram analysis + copy of Python figures \[LongRightArrow] docs/images/
wolframscript -file wolfram/run_all.wls

# 3. Rebuild this notebook
wolframscript -file community/build_notebook.wls"],

para[{bold["Licence and provenance. "], "All code and every figure in this \
notebook are original and generated by the accompanying scripts; the \
notebook and figures are released under the MIT licence (see ",
   mono["LICENSE"], " in the repository). Suggested Wolfram Community tags: ",
   ital["Biology, Physics, Differential Equations, NDSolve, Phase \
Transitions, Origin of Life, Soft Matter"], "."}],

para[{Style["Generated on " <> DateString[Now, {"MonthName", " ",
    "DayShort", ", ", "Year"}] <> " with Wolfram Language " <>
    ToString[$VersionNumber] <> ".",
    FontSlant -> Italic, FontColor -> GrayLevel[0.4]]}]
}];


(* ================================================================
   SAVE
   ================================================================ *)

Print["assembled ", Length[$allCells], " cells"];

(* RiverNetworkStatistics stylesheet: teal "Roboto Condensed" headings,
   narrow text margins, gray Abstract. Mirrors the StyleDefinitions of
   community.wolfram.com post 125916155 so this notebook reads in the
   same house style. *)
riverStyles = Notebook[{
   Cell[StyleData[StyleDefinitions -> "Default.nb"]],
   Cell[StyleData["Text"], CellMargins -> {{0, 10}, {7, 8}}],
   Cell[StyleData["Title"], CellMargins -> {{0, 3}, {15, 0}},
      FontFamily -> "Roboto Condensed",
      FontColor -> RGBColor[0.153, 0.51, 0.64]],
   Cell[StyleData["Subtitle"], CellMargins -> {{0, 3}, {20, -10}},
      FontSize -> 17, FontColor -> RGBColor[0.153, 0.51, 0.64]],
   Cell[StyleData["Abstract"], CellMargins -> {{0, 10}, {0, 0}},
      FontColor -> GrayLevel[0.25]],
   Cell[StyleData["Section"], CellMargins -> {{0, 3}, {10, 18}},
      FontFamily -> "Roboto Condensed", FontSize -> 30,
      FontColor -> RGBColor[0.153, 0.51, 0.64]],
   Cell[StyleData["Subsection"], CellMargins -> {{0, 3}, {8, 15}},
      FontFamily -> "Roboto Condensed",
      FontColor -> RGBColor[0.153, 0.51, 0.64]],
   Cell[StyleData["Subsubsection"], CellMargins -> {{0, 3}, {8, 12}},
      FontFamily -> "Roboto Condensed",
      FontColor -> RGBColor[0.153, 0.51, 0.64]],
   Cell[StyleData["Subsubsubsection"], CellMargins -> {{0, 3}, {8, 10}},
      FontFamily -> "Roboto Condensed",
      FontColor -> RGBColor[0.153, 0.51, 0.64]],
   Cell[StyleData["Item"], CellMargins -> {{22, 10.}, {4., 8.}}],
   Cell[StyleData["ItemNumbered"], CellMargins -> {{53, 10.}, {4., 8.}}]},
   Visible -> False];

nbExpr = Notebook[$allCells,
   WindowSize -> {980, 800},
   StyleDefinitions -> riverStyles];

nbBody = ToString[nbExpr, InputForm, PageWidth -> Infinity,
   CharacterEncoding -> "ASCII"];

Print["substituting ", Length[$rasterSubstitutions],
   " raster placeholders..."];
nbBody = Fold[
   StringReplace[#1, #2[[1]] -> #2[[2]]] &,
   nbBody,
   Normal[$rasterSubstitutions]
];

nbHeader = StringJoin[
   "(* Content-type: application/vnd.wolfram.mathematica *)\n\n",
   "(*** Wolfram Notebook File ***)\n",
   "(* http://www.wolfram.com/nb *)\n\n",
   "(* CreatedBy='biomolecular-condensates build pipeline' *)\n\n",
   "(* Beginning of Notebook Content *)\n"
];
nbFooter = "\n(* End of Notebook Content *)\n";

Export[outNb, nbHeader <> nbBody <> nbFooter, "Text",
   CharacterEncoding -> "ASCII"];
Print["wrote ", outNb, " (", FileByteCount[outNb], " B)"];

(* Evaluate the notebook so every figure output is generated *in place*
   (the house rule: outputs must be produced within the notebook), then
   save the evaluated notebook and export the PDF from it. *)
TimeConstrained[
   UsingFrontEnd[
      (* Backstop: load the figure package into the evaluation kernel so the
         figXxx[] cells render even if the notebook's own init cell cannot
         resolve its path headlessly. *)
      Quiet@Get[FileNameJoin[{repoRoot, "wolfram", "condensate_figures.wl"}]];
      nbObj = NotebookOpen[outNb, Visible -> False];
      Print["evaluating notebook (generating figure outputs in place)..."];
      NotebookEvaluate[nbObj, InsertResults -> True];
      NotebookSave[nbObj];
      Print["wrote evaluated ", outNb, " (", FileByteCount[outNb], " B)"];
      Quiet @ Check[
         Export[outPdf, nbObj, ImageResolution -> 300];
         Print["wrote ", outPdf, " (", FileByteCount[outPdf], " B)"],
         Print["warn: PDF export failed"]];
      NotebookClose[nbObj];
   ],
   1800,
   Print["warn: evaluation/PDF timed out"]
];
