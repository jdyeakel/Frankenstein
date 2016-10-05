(* Content-type: application/vnd.wolfram.cdf.text *)

(*** Wolfram CDF File ***)
(* http://www.wolfram.com/cdf *)

(* CreatedBy='Mathematica 10.0' *)

(*************************************************************************)
(*                                                                       *)
(*  The Mathematica License under which this file was created prohibits  *)
(*  restricting third parties in receipt of this file from republishing  *)
(*  or redistributing it by any means, including but not limited to      *)
(*  rights management or terms of use, without the express consent of    *)
(*  Wolfram Research, Inc. For additional information concerning CDF     *)
(*  licensing and redistribution see:                                    *)
(*                                                                       *)
(*        www.wolfram.com/cdf/adopting-cdf/licensing-options.html        *)
(*                                                                       *)
(*************************************************************************)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[      1064,         20]
NotebookDataLength[     20708,        505]
NotebookOptionsPosition[     21257,        500]
NotebookOutlinePosition[     21619,        516]
CellTagsIndexPosition[     21576,        513]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell[BoxData[{
 RowBox[{
  RowBox[{"Off", "[", 
   RowBox[{"InterpolatingFunction", "::", "dmval"}], "]"}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"CompDynStart", "[", 
    RowBox[{
    "rh_", ",", "rm_", ",", "ahm_", ",", "c_", ",", "k1_", ",", 
     "InitialHPower_", ",", "T_"}], "]"}], " ", ":=", "\[IndentingNewLine]", 
   "\[IndentingNewLine]", 
   RowBox[{"NDSolve", "[", 
    RowBox[{
     RowBox[{"{", "\[IndentingNewLine]", 
      RowBox[{
       RowBox[{
        RowBox[{
         RowBox[{"H", "'"}], "[", "t", "]"}], " ", "\[Equal]", 
        RowBox[{"rh", "*", 
         RowBox[{"H", "[", "t", "]"}], "*", 
         RowBox[{"(", 
          RowBox[{"1", "-", 
           FractionBox[
            RowBox[{
             RowBox[{"H", "[", "t", "]"}], "+", 
             RowBox[{"ahm", "*", 
              RowBox[{"M", "[", "t", "]"}]}]}], "k1"]}], ")"}]}]}], ",", 
       "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{
         RowBox[{"M", "'"}], "[", "t", "]"}], " ", "\[Equal]", 
        RowBox[{"rm", "*", 
         RowBox[{"M", "[", "t", "]"}], "*", 
         RowBox[{"(", 
          RowBox[{"1", "-", 
           FractionBox[
            RowBox[{
             RowBox[{"M", "[", "t", "]"}], "+", 
             RowBox[{
              RowBox[{"(", 
               FractionBox["ahm", "c"], ")"}], "*", 
              RowBox[{"H", "[", "t", "]"}]}]}], "k1"]}], ")"}]}]}], ",", 
       "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{"H", "[", "0", "]"}], "\[Equal]", 
        RowBox[{"10", "^", "InitialHPower"}]}], ",", 
       RowBox[{
        RowBox[{"M", "[", "0", "]"}], "\[Equal]", "2"}], ",", 
       "\[IndentingNewLine]", 
       RowBox[{"WhenEvent", "[", "\[IndentingNewLine]", 
        RowBox[{
         RowBox[{
          RowBox[{"M", "[", "t", "]"}], "\[Equal]", 
          RowBox[{"0.9", "*", 
           RowBox[{"(", "k1", ")"}]}]}], ",", "\[IndentingNewLine]", 
         RowBox[{
          RowBox[{"end", "=", "t"}], ";", "\[IndentingNewLine]", 
          RowBox[{"Hend", " ", "=", " ", 
           RowBox[{"H", "[", "t", "]"}]}], ";", "\[IndentingNewLine]", 
          RowBox[{"Mend", " ", "=", " ", 
           RowBox[{"M", "[", "t", "]"}]}], ";", "\[IndentingNewLine]", 
          "\"\<StopIntegration\>\""}]}], "\[IndentingNewLine]", "]"}]}], 
      "\[IndentingNewLine]", "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"H", ",", "M"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"t", ",", "0", ",", "T"}], "}"}], ",", 
     RowBox[{"MaxSteps", "\[Rule]", "Infinity"}], ",", 
     RowBox[{"InterpolationOrder", "\[Rule]", "All"}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"CompDynEnd", "[", 
    RowBox[{
    "rh_", ",", "rm_", ",", "ahm_", ",", "c_", ",", "k2_", ",", "T_", ",", 
     "Mstart_", ",", "Tstart_"}], "]"}], " ", ":=", "\[IndentingNewLine]", 
   "\[IndentingNewLine]", 
   RowBox[{"NDSolve", "[", 
    RowBox[{
     RowBox[{"{", "\[IndentingNewLine]", 
      RowBox[{
       RowBox[{
        RowBox[{
         RowBox[{"H", "'"}], "[", "t", "]"}], " ", "\[Equal]", 
        RowBox[{"rh", "*", 
         RowBox[{"H", "[", "t", "]"}], "*", 
         RowBox[{"(", 
          RowBox[{"1", "-", 
           FractionBox[
            RowBox[{
             RowBox[{"H", "[", "t", "]"}], "+", 
             RowBox[{"ahm", "*", 
              RowBox[{"M", "[", "t", "]"}]}]}], "k2"]}], ")"}]}]}], ",", 
       "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{
         RowBox[{"M", "'"}], "[", "t", "]"}], " ", "\[Equal]", 
        RowBox[{"rm", "*", 
         RowBox[{"M", "[", "t", "]"}], "*", 
         RowBox[{"(", 
          RowBox[{"1", "-", 
           FractionBox[
            RowBox[{
             RowBox[{"M", "[", "t", "]"}], "+", 
             RowBox[{
              RowBox[{"(", 
               FractionBox["ahm", "c"], ")"}], "*", 
              RowBox[{"H", "[", "t", "]"}]}]}], "k2"]}], ")"}]}]}], ",", 
       "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{"H", "[", "Tstart", "]"}], "\[Equal]", 
        RowBox[{"10", "*", 
         RowBox[{"10", "^", "9"}]}]}], ",", 
       RowBox[{
        RowBox[{"M", "[", "Tstart", "]"}], "\[Equal]", "Mstart"}]}], 
      "\[IndentingNewLine]", "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"H", ",", "M"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"t", ",", "Tstart", ",", "T"}], "}"}], ",", 
     RowBox[{"MaxSteps", "\[Rule]", "Infinity"}], ",", 
     RowBox[{"InterpolationOrder", "\[Rule]", "All"}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"Manipulate", "[", "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"With", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       RowBox[{"rh", "=", "0.0067"}], ",", 
       RowBox[{"rm", "=", 
        RowBox[{"0.0067", "*", "1.5"}]}], ",", 
       RowBox[{"ahm", "=", "8"}], ",", 
       RowBox[{"c", "=", "10"}], ",", 
       RowBox[{"InitialHPower", "=", "IH"}], ",", 
       RowBox[{"k1", "=", "K1"}], ",", 
       RowBox[{"k2", "=", 
        RowBox[{"10", "*", 
         RowBox[{"10", "^", "9"}]}]}], ",", 
       RowBox[{"T", "=", 
        RowBox[{"10", "^", "4"}]}], ",", 
       RowBox[{"Mstart", "=", "Mend"}], ",", 
       RowBox[{"Tstart", "=", "end"}]}], "}"}], ",", "\[IndentingNewLine]", 
     RowBox[{
      RowBox[{"Quiet", "[", 
       RowBox[{"Traj1", "=", 
        RowBox[{"Evaluate", "[", 
         RowBox[{
          RowBox[{"{", 
           RowBox[{
            RowBox[{"H", "[", "t", "]"}], ",", 
            RowBox[{"M", "[", "t", "]"}]}], "}"}], "/.", 
          RowBox[{"CompDynStart", "[", 
           RowBox[{
           "rh", ",", "rm", ",", "ahm", ",", "c", ",", "k1", ",", 
            "InitialHPower", ",", "T"}], "]"}]}], " ", "]"}]}], "]"}], ";", 
      "\[IndentingNewLine]", 
      RowBox[{"Quiet", "[", 
       RowBox[{"Traj2", "=", 
        RowBox[{"Evaluate", "[", 
         RowBox[{
          RowBox[{"{", 
           RowBox[{
            RowBox[{"H", "[", "t", "]"}], ",", 
            RowBox[{"M", "[", "t", "]"}]}], "}"}], "/.", 
          RowBox[{"CompDynEnd", "[", 
           RowBox[{
           "rh", ",", "rm", ",", "ahm", ",", "c", ",", "k2", ",", "T", ",", 
            "Mstart", ",", "Tstart"}], "]"}]}], " ", "]"}]}], "]"}], ";", 
      "\[IndentingNewLine]", 
      RowBox[{"Show", "[", 
       RowBox[{"{", "\[IndentingNewLine]", 
        RowBox[{
         RowBox[{"LogPlot", "[", 
          RowBox[{"Traj1", ",", 
           RowBox[{"{", 
            RowBox[{"t", ",", "0", ",", "end"}], "}"}], ",", 
           RowBox[{"PlotRange", "\[Rule]", 
            RowBox[{"{", 
             RowBox[{
              RowBox[{"{", 
               RowBox[{"0", ",", "T"}], "}"}], ",", 
              RowBox[{"{", 
               RowBox[{"1", ",", 
                RowBox[{"10", "*", 
                 RowBox[{"10", "^", "12"}]}]}], "}"}]}], "}"}]}], ",", 
           RowBox[{"AspectRatio", "\[Rule]", "1"}], ",", 
           RowBox[{"Frame", "\[Rule]", "True"}], ",", 
           RowBox[{"FrameLabel", "\[Rule]", 
            RowBox[{"{", 
             RowBox[{
             "\"\<Years after 1816\>\"", ",", "\"\<Population size\>\""}], 
             "}"}]}], ",", 
           RowBox[{"PlotStyle", "\[Rule]", 
            RowBox[{"{", 
             RowBox[{
              RowBox[{"ColorData", "[", 
               RowBox[{"97", ",", "1"}], "]"}], ",", 
              RowBox[{"ColorData", "[", 
               RowBox[{"97", ",", "3"}], "]"}]}], "}"}]}], ",", 
           RowBox[{"PlotLegends", "\[Rule]", 
            RowBox[{"LineLegend", "[", 
             RowBox[{"{", 
              RowBox[{"\"\<Humans\>\"", ",", "\"\<Creatures\>\""}], "}"}], 
             "]"}]}]}], "]"}], ",", "\[IndentingNewLine]", 
         RowBox[{"LogPlot", "[", 
          RowBox[{"Traj2", ",", 
           RowBox[{"{", 
            RowBox[{"t", ",", "end", ",", "T"}], "}"}], ",", 
           RowBox[{"PlotRange", "\[Rule]", 
            RowBox[{"{", 
             RowBox[{"1", ",", "All"}], "}"}]}], ",", 
           RowBox[{"PlotStyle", "\[Rule]", 
            RowBox[{"{", 
             RowBox[{
              RowBox[{"ColorData", "[", 
               RowBox[{"97", ",", "1"}], "]"}], ",", 
              RowBox[{"ColorData", "[", 
               RowBox[{"97", ",", "3"}], "]"}]}], "}"}]}]}], "]"}], ",", 
         "\[IndentingNewLine]", 
         RowBox[{"Graphics", "[", 
          RowBox[{"{", 
           RowBox[{
            RowBox[{"ColorData", "[", 
             RowBox[{"97", ",", "1"}], "]"}], ",", "Dashed", ",", 
            RowBox[{"Thickness", "[", "0.005", "]"}], ",", 
            RowBox[{"Line", "[", 
             RowBox[{"{", 
              RowBox[{
               RowBox[{"{", 
                RowBox[{"end", ",", 
                 RowBox[{"Log", "@", "Hend"}]}], "}"}], ",", 
               RowBox[{"{", 
                RowBox[{"end", ",", 
                 RowBox[{"Log", "@", 
                  RowBox[{"(", 
                   RowBox[{"10", "*", 
                    RowBox[{"10", "^", "9"}]}], ")"}]}]}], "}"}]}], "}"}], 
             "]"}]}], "}"}], "]"}], ",", "\[IndentingNewLine]", 
         RowBox[{"Graphics", "[", 
          RowBox[{"Text", "[", 
           RowBox[{"\"\<Global Comp.\>\"", ",", 
            RowBox[{"{", 
             RowBox[{"end", ",", 
              RowBox[{"Log", "@", 
               RowBox[{"(", 
                RowBox[{"10", "^", "11"}], ")"}]}]}], "}"}]}], "]"}], "]"}], 
         ",", "\[IndentingNewLine]", 
         RowBox[{"Graphics", "[", 
          RowBox[{"Arrow", "[", 
           RowBox[{"{", 
            RowBox[{
             RowBox[{"{", 
              RowBox[{"end", ",", 
               RowBox[{"Log", "@", 
                RowBox[{"(", 
                 RowBox[{"10", "^", "10.9"}], ")"}]}]}], "}"}], ",", 
             RowBox[{"{", 
              RowBox[{"end", ",", 
               RowBox[{"Log", "@", 
                RowBox[{"(", 
                 RowBox[{"10", "^", "10.1"}], ")"}]}]}], "}"}]}], "}"}], 
           "]"}], "]"}]}], "\[IndentingNewLine]", "}"}], "]"}]}]}], 
    "\[IndentingNewLine]", "]"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"K1", ",", 
       RowBox[{"600", "*", 
        RowBox[{"10", "^", "6"}]}], ",", 
       "\"\<Initial carrying capacity\>\""}], "}"}], ",", "10", ",", 
     RowBox[{"600", "*", 
      RowBox[{"10", "^", "6"}]}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"IH", ",", "8.8", ",", "\"\<Initial human population\>\""}], 
      "}"}], ",", "0", ",", "10"}], "}"}], ",", "\[IndentingNewLine]", 
   RowBox[{"(*", 
    RowBox[{"bookmark", " ", "control"}], "*)"}], "\[IndentingNewLine]", 
   RowBox[{"Grid", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"Style", "[", 
         RowBox[{"\"\<Scenario\>\"", ",", "Bold"}], "]"}], ",", 
        RowBox[{"RadioButtonBar", "[", 
         RowBox[{
          RowBox[{"Dynamic", "[", 
           RowBox[{"bookMark", ",", 
            RowBox[{
             RowBox[{"{", 
              RowBox[{
               RowBox[{"bookMark", "=", "#"}], ";", "\[IndentingNewLine]", 
               RowBox[{"Which", "[", "\[IndentingNewLine]", 
                RowBox[{
                 RowBox[{"bookMark", "\[Equal]", "\"\<Europe\>\""}], ",", 
                 RowBox[{
                  RowBox[{"K1", "=", 
                   RowBox[{"600", "*", 
                    RowBox[{"10", "^", "6"}]}]}], ";", 
                  RowBox[{"IH", "=", "8.25"}]}], ",", "\[IndentingNewLine]", 
                 RowBox[{"bookMark", "\[Equal]", "\"\<Amazon\>\""}], ",", 
                 RowBox[{
                  RowBox[{"K1", "=", 
                   RowBox[{"1.383", "*", 
                    RowBox[{"10", "^", "6"}]}]}], ";", 
                  RowBox[{"IH", "=", "4.64"}]}]}], "\[IndentingNewLine]", 
                "]"}]}], "}"}], "&"}]}], "]"}], ",", 
          RowBox[{"{", 
           RowBox[{"\"\<Europe\>\"", ",", "\"\<Amazon\>\""}], "}"}], ",", 
          RowBox[{"Appearance", "\[Rule]", "\"\<Horizontal\>\""}]}], "]"}]}], 
       "}"}], "}"}], ",", 
     RowBox[{"Spacings", "\[Rule]", 
      RowBox[{"{", 
       RowBox[{".5", ",", "1"}], "}"}]}], ",", 
     RowBox[{"Frame", "\[Rule]", "False"}]}], "]"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"bookMark", ",", "\"\<Europe\>\""}], "}"}], ",", "None"}], 
    "}"}]}], "\[IndentingNewLine]", "]"}]}], "Input",
 CellOpen->False,
 CellChangeTimes->{{3.684596333696332*^9, 3.684596345204221*^9}, {
   3.684596383117826*^9, 3.684596415427936*^9}, {3.6845965216867447`*^9, 
   3.684596625749658*^9}, {3.684596686739855*^9, 3.6845967308954277`*^9}, {
   3.6845967868873568`*^9, 3.684596825260557*^9}, {3.6845999694379673`*^9, 
   3.684599990393537*^9}, {3.684600051021625*^9, 3.684600051221957*^9}, {
   3.6846002076233997`*^9, 3.6846005384896584`*^9}, {3.684600591773375*^9, 
   3.684600593528693*^9}, {3.684600649462737*^9, 3.684600652037298*^9}, {
   3.684600870449545*^9, 3.684600872115103*^9}, {3.6846009492807503`*^9, 
   3.6846010217983713`*^9}, {3.6846013939097233`*^9, 3.684601489865144*^9}, {
   3.684601521548088*^9, 3.684601675859783*^9}, {3.6846017635671673`*^9, 
   3.684601903149778*^9}, 3.6846019737291517`*^9, {3.684602170475642*^9, 
   3.6846021727573357`*^9}, {3.6846047970944633`*^9, 3.684604814866477*^9}, {
   3.684604854614447*^9, 3.6846048849448547`*^9}, {3.684604916648595*^9, 
   3.684604920414679*^9}, {3.684604962009788*^9, 3.6846049626385803`*^9}, {
   3.6846052408471403`*^9, 3.684605242257503*^9}, {3.684608919268716*^9, 
   3.684608960673057*^9}, {3.6846090250817623`*^9, 3.684609025524975*^9}, {
   3.6846090846563177`*^9, 3.6846090905884933`*^9}, {3.6846122496059713`*^9, 
   3.684612270639015*^9}, {3.684672587133233*^9, 3.684672587227708*^9}, {
   3.684672636934908*^9, 3.684672784669301*^9}, {3.6846730760017767`*^9, 
   3.684673146493273*^9}, {3.684673355060133*^9, 3.68467341707698*^9}, {
   3.684673589375586*^9, 3.684673644763959*^9}, {3.684673717884748*^9, 
   3.6846737651060343`*^9}, {3.684673816221916*^9, 3.684673889977099*^9}, {
   3.684674362920289*^9, 3.684674371152234*^9}, {3.684674877909853*^9, 
   3.684674971217514*^9}, {3.6846754330704203`*^9, 3.68467544183438*^9}}],

Cell[BoxData[
 TagBox[
  StyleBox[
   DynamicModuleBox[{$CellContext`bookMark$$ = "Europe", $CellContext`IH$$ = 
    3.68, $CellContext`K1$$ = 4.41*^8, Typeset`show$$ = True, 
    Typeset`bookmarkList$$ = {}, Typeset`bookmarkMode$$ = "Menu", 
    Typeset`animator$$, Typeset`animvar$$ = 1, Typeset`name$$ = 
    "\"untitled\"", Typeset`specs$$ = {{{
       Hold[$CellContext`K1$$], 600000000, "Initial carrying capacity"}, 10, 
      600000000}, {{
       Hold[$CellContext`IH$$], 8.8, "Initial human population"}, 0, 10}, {
      Hold[
       Grid[{{
          Style["Scenario", Bold], 
          RadioButtonBar[
           
           Dynamic[$CellContext`bookMark$$, {$CellContext`bookMark$$ = #; 
             Which[$CellContext`bookMark$$ == 
               "Europe", $CellContext`K1$$ = 600 10^6; $CellContext`IH$$ = 
                8.25, $CellContext`bookMark$$ == 
               "Amazon", $CellContext`K1$$ = 1.383 10^6; $CellContext`IH$$ = 
                4.64]}& ], {"Europe", "Amazon"}, Appearance -> 
           "Horizontal"]}}, Spacings -> {0.5, 1}, Frame -> False]], 
      Manipulate`Dump`ThisIsNotAControl}, {{
       Hold[$CellContext`bookMark$$], "Europe"}, 0}}, Typeset`size$$ = {
    463., {170., 175.}}, Typeset`update$$ = 0, Typeset`initDone$$, 
    Typeset`skipInitDone$$ = True, $CellContext`K1$12705$$ = 
    0, $CellContext`IH$12706$$ = 0}, 
    DynamicBox[Manipulate`ManipulateBoxes[
     1, StandardForm, 
      "Variables" :> {$CellContext`bookMark$$ = "Europe", $CellContext`IH$$ = 
        8.8, $CellContext`K1$$ = 600000000}, "ControllerVariables" :> {
        Hold[$CellContext`K1$$, $CellContext`K1$12705$$, 0], 
        Hold[$CellContext`IH$$, $CellContext`IH$12706$$, 0]}, 
      "OtherVariables" :> {
       Typeset`show$$, Typeset`bookmarkList$$, Typeset`bookmarkMode$$, 
        Typeset`animator$$, Typeset`animvar$$, Typeset`name$$, 
        Typeset`specs$$, Typeset`size$$, Typeset`update$$, Typeset`initDone$$,
         Typeset`skipInitDone$$}, "Body" :> 
      With[{$CellContext`rh = 0.0067, $CellContext`rm = 
         0.0067 1.5, $CellContext`ahm = 8, $CellContext`c = 
         10, $CellContext`InitialHPower = $CellContext`IH$$, $CellContext`k1 = \
$CellContext`K1$$, $CellContext`k2 = 10 10^9, $CellContext`T = 
         10^4, $CellContext`Mstart = $CellContext`Mend, $CellContext`Tstart = \
$CellContext`end}, Quiet[$CellContext`Traj1 = Evaluate[
            ReplaceAll[{
              $CellContext`H[$CellContext`t], 
              $CellContext`M[$CellContext`t]}, 
             $CellContext`CompDynStart[$CellContext`rh, $CellContext`rm, \
$CellContext`ahm, $CellContext`c, $CellContext`k1, \
$CellContext`InitialHPower, $CellContext`T]]]]; 
        Quiet[$CellContext`Traj2 = Evaluate[
            ReplaceAll[{
              $CellContext`H[$CellContext`t], 
              $CellContext`M[$CellContext`t]}, 
             $CellContext`CompDynEnd[$CellContext`rh, $CellContext`rm, \
$CellContext`ahm, $CellContext`c, $CellContext`k2, $CellContext`T, \
$CellContext`Mstart, $CellContext`Tstart]]]]; Show[{
           
           LogPlot[$CellContext`Traj1, {$CellContext`t, 0, $CellContext`end}, 
            PlotRange -> {{0, $CellContext`T}, {1, 10 10^12}}, AspectRatio -> 
            1, Frame -> True, 
            FrameLabel -> {"Years after 1816", "Population size"}, PlotStyle -> {
              ColorData[97, 1], 
              ColorData[97, 3]}, PlotLegends -> 
            LineLegend[{"Humans", "Creatures"}]], 
           
           LogPlot[$CellContext`Traj2, {$CellContext`t, $CellContext`end, \
$CellContext`T}, PlotRange -> {1, All}, PlotStyle -> {
              ColorData[97, 1], 
              ColorData[97, 3]}], 
           Graphics[{
             ColorData[97, 1], Dashed, 
             Thickness[0.005], 
             Line[{{$CellContext`end, 
                Log[$CellContext`Hend]}, {$CellContext`end, 
                Log[10 10^9]}}]}], 
           Graphics[
            Text["Global Comp.", {$CellContext`end, 
              Log[10^11]}]], 
           Graphics[
            Arrow[{{$CellContext`end, 
               Log[10^10.9]}, {$CellContext`end, 
               Log[10^10.1]}}]]}]], 
      "Specifications" :> {{{$CellContext`K1$$, 600000000, 
          "Initial carrying capacity"}, 10, 
         600000000}, {{$CellContext`IH$$, 8.8, "Initial human population"}, 0,
          10}, 
        Grid[{{
           Style["Scenario", Bold], 
           RadioButtonBar[
            
            Dynamic[$CellContext`bookMark$$, {$CellContext`bookMark$$ = #; 
              Which[$CellContext`bookMark$$ == 
                "Europe", $CellContext`K1$$ = 600 10^6; $CellContext`IH$$ = 
                 8.25, $CellContext`bookMark$$ == 
                "Amazon", $CellContext`K1$$ = 1.383 10^6; $CellContext`IH$$ = 
                 4.64]}& ], {"Europe", "Amazon"}, Appearance -> 
            "Horizontal"]}}, Spacings -> {0.5, 1}, Frame -> 
         False], {{$CellContext`bookMark$$, "Europe"}, 0, ControlType -> 
         None}}, "Options" :> {}, "DefaultOptions" :> {}],
     ImageSizeCache->{508., {241., 246.}},
     SingleEvaluation->True],
    Deinitialization:>None,
    DynamicModuleValues:>{},
    SynchronousInitialization->True,
    UndoTrackedVariables:>{Typeset`show$$, Typeset`bookmarkMode$$},
    UnsavedVariables:>{Typeset`initDone$$},
    UntrackedVariables:>{Typeset`size$$}], "Manipulate",
   Deployed->True,
   StripOnInput->False],
  Manipulate`InterpretManipulate[1]]], "Output",
 CellChangeTimes->{3.684679489213717*^9}]
}, Open  ]]
},
WindowSize->{755, 821},
WindowMargins->{{-976, Automatic}, {Automatic, -387}},
FrontEndVersion->"10.0 for Mac OS X x86 (32-bit, 64-bit Kernel) (September 9, \
2014)",
StyleDefinitions->"Default.nb"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[CellGroupData[{
Cell[1486, 35, 14234, 348, 19, "Input",
 CellOpen->False],
Cell[15723, 385, 5518, 112, 504, "Output"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature vx02E3lZR@jvoAwaOYCBtvzp *)
