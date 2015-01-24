(* ::Package:: *)

(* Time-Stamp: <2015-01-24 13:32:06 misho> *)

(* :Context: SimpleMSSM` *)

(* :Author: Sho Iwamoto *)

(* :Summary:
    Subpackage of SimpleMSSM: Calculate MSSM mass spectrum at very tree-level.
*)

(* :Copyright: 2015 Sho Iwamoto *)

(* :Package Version: 1.1
   Version 1.0 [Jan. 2015, SI] Initial version.
   Version 1.1 [Jan. 2015, SI] Adopt to SLHA.m version 2.0.
*)

(* :Mathematica Version: 10.0 *)

(* :History: *)

(* :Keywords: *)

(* :Discussion: *)

(* :Warning: *)

BeginPackage["SimpleMSSM`"];

<<"MSSM.m";

(* Usage messages *)

CalculateSpectrum::usage = "\!\(\*RowBox[{\"CalculateSpectrum\",\"[\",StyleBox[\"slha\",\"TI\"],\",\",StyleBox[\"ofile\",\"TI\"],\"]\"}]\) receives a SLHA object \!\(\*RowBox[{StyleBox[\"slha\",\"TI\"]}]\) as input parameters and gives a SLHA object with full SLHA mass spectrum. When the output filename \!\(\*RowBox[{StyleBox[\"ofile\",\"TI\"]}]\) is specified, the result is also saved in the file.";

(* Error messages *)

MSSMparameter::tanb = "tanb(Minput) (EXTPAR 25) is used instead of MINPAR 3, but note that scale dependence is ignored.";
MSSMparameter::InvalidHiggsParameters = "Parameters for Higgs sector (EXTPAR 21-27) is invalid.";
MSSMparameter::InvalidMuSign = "Sign of mu is inconsistent. Note that Scale dependence is totally ignored.";
MSSMparameter::InvalidEWSB = "EWSB not occured.";
MSSMparameter::InvalidEWSBmasq = "EWSB not occured (Negative mA^2).";

CalculateSpectrum::SLHANotLoaded = "SLHA module not loaded.";
CalculateSpectrum::InputScaleIgnored = "Spectrum input scale (EXTPAR 0) is set but ignored.";
CalculateSpectrum::InvalidMODSEL = "MODSEL is assumed to be 1(SUGRA).";

(* Other messages *)


Begin["`Private`"];

If[Not[MemberQ[$Packages, "SLHA`"]], Message[CalculateSpectrum::SLHANotLoaded]; Abort[]];
IfMissing = SLHA`IfMissing;
IsMissing = SLHA`IsMissing;

(* SM parameters *)

SMparameterSetByHand[] :=
  Module[ {smbyhand},
  (* For simplicity, we use following values. *)
  smbyhand["v"]  = 246 // N;
  smbyhand["sw"] = Sqrt[0.23116];
  smbyhand["mm"] = 0.105658365;
  smbyhand];

SMparameter[slha_] :=
  Module[ {sm, smbyhand,
           ainv, gf, as, mz, mb, mt, mtau,
           sinw, w, v,
           ainvY0, ainvW0, ainvS0},
  {ainv, gf, as, mz, mb, mt, mtau} = slha["SMINPUTS"][#, IfMissing->"Abort"] & /@ {1,2,3,4,5,6,7};

  smbyhand = SMparameterSetByHand[];
  sm["v"]  = smbyhand["v"];
  sm["sw"] = smbyhand["sw"];
  sm["mm"] = smbyhand["mm"];

  sm["cw"] = Sqrt[1-sm["sw"]^2];
  {ainvY0, ainvW0, ainvS0} = {ainv * sm["cw"]^2, ainv * sm["sw"]^2, 1/as };
  sm["mZ"] = mz;
  sm["mW"] = mz * sm["cw"];
  sm["gY"] = Function[{Q}, Sqrt[ 4\[Pi] / (ainvY0 - (41/6) Log[Q/mz] / (2\[Pi]) )]];
  sm["g2"] = Function[{Q}, Sqrt[ 4\[Pi] / (ainvW0 + (19/6) Log[Q/mz] / (2\[Pi]) )]];
  sm["g3"] = Function[{Q}, Sqrt[ 4\[Pi] / (ainvS0 + (7)    Log[Q/mz] / (2\[Pi]) )]];
  sm["yt"] = Function[{Q}, Sqrt[2] * mt / sm["v"] * (1 - sm["g3"][Q]^2 / (3 \[Pi]^2))]; (* One-loop simple conversion *)
  sm["yb"] = Function[{Q}, Sqrt[2] * mb / sm["v"]];
  sm["ya"] = Function[{Q}, Sqrt[2] * mtau / sm["v"]]; (* calculated from pole mass; quantum correction ignored. *)
  sm["ym"] = Function[{Q}, Sqrt[2] * sm["mm"] / sm["v"]]; (* calculated from pole mass; quantum correction ignored. *)
  sm];

(* Masses *)
NeutralinoMass[sm_, mssm_] := NeutralinoMass[mssm["M1"], mssm["M2"], mssm["tb"], mssm["mu"], sm["mZ"], sm["sw"]];
CharginoMass  [sm_, mssm_] := CharginoMass  [mssm["M2"], mssm["tb"], mssm["mu"], sm["mW"]];

SUmassSQ[sm_, mssm_, flavor_] := With[
    { mF = If[flavor != 3, 0, sm["yt"][mssm["Q"]]*sm["v"]/Sqrt[2]],
      A  = If[flavor != 3, 0, mssm["At"]] },
    SUmassSQ[mF, mssm["QL", flavor], mssm["UR", flavor], mssm["mu"], A, ArcTan[mssm["tb"]], sm["mZ"], sm["sw"]]];

SDmassSQ[sm_, mssm_, flavor_] := With[
    { mF = If[flavor != 3, 0, sm["yb"][mssm["Q"]]*sm["v"]/Sqrt[2]],
      A  = If[flavor != 3, 0, mssm["Ab"]] },
    SDmassSQ[mF, mssm["QL", flavor], mssm["DR", flavor], mssm["mu"], A, ArcTan[mssm["tb"]], sm["mZ"], sm["sw"]]];

SLmassSQ[sm_, mssm_, flavor_] := With[
    { mF = If[flavor != 3, 0, sm["ya"][mssm["Q"]]*sm["v"]/Sqrt[2]],
      A  = If[flavor != 3, 0, mssm["Aa"]] },
    SLmassSQ[mF, mssm["LL", flavor], mssm["ER", flavor], mssm["mu"], A, ArcTan[mssm["tb"]], sm["mZ"], sm["sw"]]];

SNmassSQ[sm_, mssm_, flavor_] := SNmassSQ[mssm["LL", flavor], ArcTan[mssm["tb"]], sm["mZ"], sm["sw"]];


(* MSSM *)
MSSMparameter[slha_, sm_] := Module[
    {readminext, p,
     mu, m1sq, m2sq, masq, tb, mA0, mHP,
     c, mZsq, sgnmu, musq, beta, b},
    readminext[min_, ext_] := Module[{k}, k = slha["EXTPAR"][ext]; If[IsMissing[k], slha["MINPAR"][min, IfMissing->"abort"], k]] // N;

    p["tb"] = slha["MINPAR", IfMissing->"abort"][3, IfMissing->"abort"] // N;

    p["M1"] = readminext[2, 1];
    p["M2"] = readminext[2, 2];
    p["M3"] = readminext[2, 3];
    p["At"] = readminext[5, 11];
    p["Ab"] = readminext[5, 12];
    p["Aa"] = readminext[5, 13];
    p["Am"] = 0;
    Do[
        p["LL",i] = readminext[1, 30+i];
        p["ER",i] = readminext[1, 33+i];
        p["QL",i] = readminext[1, 40+i];
        p["UR",i] = readminext[1, 43+i];
        p["DR",i] = readminext[1, 46+i];
        ,  {i, 1, 3}];
    {mu, masq, tb, mA0, mHP} = slha["EXTPAR"][#] & /@ {23, 24, 25, 26, 27} // N;
    c = Total[If[IsMissing[#], 0, 1] & /@ {masq, mA0, mHP}];
    If[Not[IsMissing[tb]], Message[MSSMparameter::tanb]; p["tb"] = tb];
    beta = ArcTan[p["tb"]];
    mZsq = sm["mZ"]^2;

    If[IsMissing[mu],
       If[c != 0, Message[MSSMparameter::InvalidHiggsParameters]; Abort[]];
       m1sq = Module[{k}, k = slha["EXTPAR"][21]; If[IsMissing[k], slha["MINPAR"][1, IfMissing->"abort"]^2, k]];
       m2sq = Module[{k}, k = slha["EXTPAR"][22]; If[IsMissing[k], slha["MINPAR"][1, IfMissing->"abort"]^2, k]];
       (* Now m1sq and m2sq are set *)
       musq = (- m1sq - m2sq - mZsq + abs[(m1sq - m2sq)/Cos[2beta]]) / 2;
       If[musq < 0, Message[MSSMparameter::InvalidEWSB]; Abort[]];
       mu = Sign[slha["MINPAR"][4, IfMissing->"abort"]] * Sqrt[musq];
       masq = 2musq + m1sq + m2sq;
       b    = masq * Sin[2beta] / 2;
       If[masq < 0, Message[MSSMparameter::InvalidEWSBmasq]; Abort[]];
       ,
       If[c != 1, Message[MSSMparameter::InvalidHiggsParameters]; Abort[]];
       If[Not[IsMissing[mA0]], masq = mA0^2];
       If[Not[IsMissing[mHP]], masq = mHP^2 - sm["mW"]^2];
       sgnmu = slha["MINPAR"][4];
       If[And[Not[IsMissing[sgnmu]], sgnmu * mu < 0], Message[MSSMparameter::InvalidMuSign]; Abort[]];
       (* Now mu and masq are set *)
       b = masq * Sin[2beta] / 2;
       m2sq = - mu^2 + b Cot[beta] + mZsq Cos[2beta] / 2;
       m1sq = - mu^2 + b Tan[beta] - mZsq Cos[2beta] / 2;
      ];
    (* Now m1sq, m2sq, mu, b, and masq are prepared. *)

    (* ALPHA, HMIX and MSOFT; note that p["tb"] is already set, and p["v"] is not calculated but used sm["v"]. *)
    p["a"]    = ArcTan[Tan[2beta]*(masq + mZsq)/(masq-mZsq)]/2;
    p["mu"]   = mu;
    p["v"]    = sm["v"];
    p["mAsq"] = masq;
    p["m1sq"] = m1sq;
    p["m2sq"] = m2sq;

    (* MASS *)
    c = (masq - mZsq)^2 + 4*mZsq*masq*Sin[2beta]^2;
    If[c < 0, Message[MSSMparameter::InvalidEWSB]; Abort[]];
    p[24] = sm["mW"];
    p[25] = Sqrt[(masq + mZsq - Sqrt[c]) / 2];
    p[35] = Sqrt[(masq + mZsq + Sqrt[c]) / 2];
    p[36] = Sqrt[masq];
    p[37] = Sqrt[masq + sm["mW"]^2];
    p];

GenNewBlock[name_, headcomment_, list_, Q_:Null] := Module[
    {b = SLHA`Private`NewBlock[name, headcomment]},
    AddValue[b, #]& /@ list;
    If[NumberQ[Q], b["Q"] = Q];
    b];
AddValue[b_, {v_, c_}] := Module[
    {value = Last[v], key = Most[v]},
    b[Sequence@@key] = value;
    b["c", Sequence@@key] = c];

GenSPINFO[] := GenNewBlock["SPINFO", "# Program Information",
    {{{1, "TreeSpec"}, "# spectrum calculator"},
     {{2, "1.1"},      "# version number"}}];

GenMODSEL[] := GenNewBlock["MODSEL", "# model selection", {{{1, 1}, "# SUGRA-based general MSSM simulation"}}];

GenMass[mssm_] := GenNewBlock["MASS", "# mass spectrum",
    {{{     24, mssm[     24]}, "# W"},
     {{     25, mssm[     25]}, "# h (set by hand); tree level = " <> ToString[mssm[-25]]},
     {{     35, mssm[     35]}, "# H0"},
     {{     36, mssm[     36]}, "# A0"},
     {{     37, mssm[     37]}, "# H+"},
     {{1000022, mssm[1000022]}, "# N1"},
     {{1000023, mssm[1000023]}, "# N2"},
     {{1000025, mssm[1000025]}, "# N3"},
     {{1000035, mssm[1000035]}, "# N4"},
     {{1000024, mssm[1000024]}, "# C1"},
     {{1000037, mssm[1000037]}, "# C2"},
     {{1000011, mssm[1000011]}, "# ~e_L"},
     {{1000013, mssm[1000013]}, "# ~mu_L"},
     {{2000011, mssm[2000011]}, "# ~e_R"},
     {{2000013, mssm[2000013]}, "# ~mu_R"},
     {{1000015, mssm[1000015]}, "# ~tau_1"},
     {{2000015, mssm[2000015]}, "# ~tau_2"},
     {{1000012, mssm[1000012]}, "# ~nu_e"},
     {{1000014, mssm[1000014]}, "# ~nu_mu"},
     {{1000016, mssm[1000016]}, "# ~nu_tau"},
     {{1000021, mssm[1000021]}, "# gluino (tree!)"},
     {{1000001, mssm[1000001]}, "# ~d_L"},
     {{1000002, mssm[1000002]}, "# ~u_L"},
     {{1000003, mssm[1000003]}, "# ~s_L"},
     {{1000004, mssm[1000004]}, "# ~c_L"},
     {{2000001, mssm[2000001]}, "# ~d_R"},
     {{2000002, mssm[2000002]}, "# ~u_R"},
     {{2000003, mssm[2000003]}, "# ~s_R"},
     {{2000004, mssm[2000004]}, "# ~c_R"},
     {{1000005, mssm[1000005]}, "# ~b_1"},
     {{2000005, mssm[2000005]}, "# ~b_2"},
     {{1000006, mssm[1000006]}, "# ~t_1"},
     {{2000006, mssm[2000006]}, "# ~t_2"}}];

GenALPHA[mssm_] := GenNewBlock["ALPHA", "# Effective Higgs mixing parameter", {{{mssm["a"]}, "# alpha"}}];

GenHMIX[mssm_] := GenNewBlock["HMIX", "# Higgs mixing parameters",
    {{{1, mssm["mu"]},   "# mu (tree-level)"},
     {{2, mssm["tb"]},   "# tb (tree-level)"},
     {{3, mssm["v"]},    "# Higgs vev (input by hand)"},
     {{4, mssm["mAsq"]}, "# mA^2 (tree-level)"}}, mssm["Q"]];

GenYAsub[name_, Q_, value_, comment_] := GenNewBlock[name, "",
      {(*{{1, 1, 0.}, ""}, *)
       (*{{2, 2, 0.}, ""}, *)
         {{3, 3, value}, comment}},
      Q];

GenYtYbYa[sm_, mssm_] :=
    { GenYAsub["YU", mssm["Q"], sm["yt"][mssm["Q"]]/Sin[ArcTan[mssm["tb"]]], "# Yt (tree; 1loop for Pole->MSbar)"],
      GenYAsub["YD", mssm["Q"], sm["yb"][mssm["Q"]]/Cos[ArcTan[mssm["tb"]]], "# Yb (tree; from msbar mass)"],
      GenYAsub["YE", mssm["Q"], sm["ya"][mssm["Q"]]/Cos[ArcTan[mssm["tb"]]], "# Ytau (tree; from pole mass)"] };

GenAtAbAa[sm_, mssm_] :=
    { GenYAsub["AU", mssm["Q"], mssm["At"], "# At (tree; i.e., = input value)"],
      GenYAsub["AD", mssm["Q"], mssm["Ab"], "# Ab (tree; i.e., = input value)"],
      GenYAsub["AE", mssm["Q"], mssm["Aa"], "# Atau (tree; i.e., = input value)"] };

GenMix[name_, matrix_, comment_] := GenNewBlock[name, comment,
    MapIndexed[Function[{v, i}, MapIndexed[Function[{v2, j}, {Flatten[{i, j, v2}], ""}], v]], matrix]//Flatten[#,1] &];

GenGAUGE[sm_] := GenNewBlock["GAUGE", "# Gauge couplings (at mZ)",
    {{{1, sm["gY"][sm["mZ"]]}, "# gY(mZ) MSbar"},
     {{2, sm["g2"][sm["mZ"]]}, "# g2(mZ) MSbar"},
     {{3, sm["g3"][sm["mZ"]]}, "# g3(mZ) MSbar"}}, sm["mZ"]];

GenMSOFT[mssm_] := GenNewBlock["MSOFT", "# MSSM SUSY-breaking (tree-level, i.e., = input)", 
    Join[{ {{ 1, mssm["M1"  ]}, "# M_1"},
           {{ 2, mssm["M2"  ]}, "# M_2"},
           {{ 3, mssm["M3"  ]}, "# M_3"},
           {{21, mssm["m1sq"]}, "# mH1^2 = mHd^2"},
           {{22, mssm["m2sq"]}, "# mH2^2 = mHu^2"}
         },
         {{30+#, mssm["LL",#]}, "# m_LL["<>ToString[#]<>"]"} &/@{1,2,3},
         {{33+#, mssm["ER",#]}, "# m_eR["<>ToString[#]<>"]"} &/@{1,2,3},
         {{40+#, mssm["QL",#]}, "# m_QL["<>ToString[#]<>"]"} &/@{1,2,3},
         {{43+#, mssm["UR",#]}, "# m_uR["<>ToString[#]<>"]"} &/@{1,2,3},
         {{46+#, mssm["DR",#]}, "# m_dR["<>ToString[#]<>"]"} &/@{1,2,3}], mssm["Q"]];

Options[WriteBlock] = {Comment->False, Order->{}};
WriteBlock[ofs_, b_, OptionsPattern[]] := Module[
    {list = b["tostring", Order->OptionValue[Order]]},
    If[IsMissing[list], $Failed,
       If[OptionValue[Comment], list = StringReplace[list, RegularExpression["^ ?"]->"#"]];
       WriteString[ofs, #<>"\n"] &/@ Append[list, "#"]]];

CalculateSpectrum[slha_, OFile_] := Module[
    {sm, mssm,
     ofs,
     nmix, umix, vmix, stopmix, sbotmix, staumix, diag},
    sm = SMparameter[slha];

    (* MSSM parameters *)
    If[slha["MODSEL", 1] =!= 1, Message[CalculateSpectrum::InvalidMODSEL]];
    If[Not[SLHA`IsMissing[slha["EXTPAR", 0]]], Message[CalculateSpectrum::InputScaleIgnored]];
    mssm = MSSMparameter[slha, sm];
    mssm["Q"] = 200.; (* Simplicity *)

    (* masses *)
    {diag, nmix} = NeutralinoRealDiagonalize[NeutralinoMass[sm, mssm]];
    {mssm[1000022], mssm[1000023], mssm[1000025], mssm[1000035]} = Diagonal[diag];
    nmix = GenMix["NMIX", nmix, "# Neutralino mixing matrix"];

    {diag, umix, vmix} = CharginoDiagonalize[CharginoMass[sm, mssm]];
    {mssm[1000024], mssm[1000037]} = Diagonal[diag];
    umix = GenMix["UMIX", umix, "# Chargino U mixing matrix"];
    vmix = GenMix["VMIX", vmix, "# Chargino V mixing matrix"];

    {diag, stopmix} = SfermionDiagonalize[SUmassSQ[sm, mssm, 3]];
    {mssm[1000006], mssm[2000006]} = Diagonal[diag] // Sqrt;
    {diag, sbotmix} = SfermionDiagonalize[SDmassSQ[sm, mssm, 3]];
    {mssm[1000005], mssm[2000005]} = Diagonal[diag] // Sqrt;
    {diag, staumix} = SfermionDiagonalize[SLmassSQ[sm, mssm, 3]];
    {mssm[1000015], mssm[2000015]} = Diagonal[diag] // Sqrt;
    stopmix = GenMix["STOPMIX", stopmix, "# Stop mixing matrix"];
    sbotmix = GenMix["SBOTMIX", sbotmix, "# Sbottom mixing matrix"];
    staumix = GenMix["STAUMIX", staumix, "# Stau V mixing matrix"];

    {mssm[1000011], mssm[2000011]} = Diagonal[SLmassSQ[sm, mssm, 1]] // Sqrt; (* TODO: To meet MG5 requirement that selectron and smuon have the exactly same mass *)
    {mssm[1000013], mssm[2000013]} = Diagonal[SLmassSQ[sm, mssm, 2]] // Sqrt;
    {mssm[1000001], mssm[2000001]} = Diagonal[SDmassSQ[sm, mssm, 1]] // Sqrt;
    {mssm[1000003], mssm[2000003]} = Diagonal[SDmassSQ[sm, mssm, 2]] // Sqrt;
    {mssm[1000002], mssm[2000002]} = Diagonal[SUmassSQ[sm, mssm, 1]] // Sqrt;
    {mssm[1000004], mssm[2000004]} = Diagonal[SUmassSQ[sm, mssm, 2]] // Sqrt;
    {mssm[1000012], mssm[1000014], mssm[1000016]} = SNmassSQ[sm, mssm, #] &/@{1,2,3} // Sqrt;

    (* masses set by hand*)
    mssm[-25] = mssm[25];
    mssm[25] = 126.;
    mssm[1000021] = mssm["M3"];

    (* output *)
    ofs = OpenWrite[OFile, PageWidth->Infinity];
    WriteBlock[ofs, #] &/@ {GenSPINFO[], slha["SMINPUTS"]};
    WriteBlock[ofs, #, Comment->True] &/@ {GenMODSEL[], slha["MINPAR"], slha["EXTPAR"]};
    WriteBlock[ofs, GenMass[mssm],
               Order->{24, 25, 35, 36, 37, 1000021, 1000001, 1000002, 1000003, 1000004, 2000001, 2000002, 2000003, 2000004, 1000005, 2000005, 1000006, 2000006,
                       1000022, 1000023, 1000025, 1000035, 1000024, 1000037, 1000011, 1000013, 1000012, 1000014, 2000011, 2000013, 1000015, 2000015, 1000016}];
    WriteBlock[ofs, #] &/@ Flatten[{GenALPHA[mssm], GenHMIX[mssm],
                                    nmix, umix, vmix, stopmix, sbotmix, staumix,
                                    GenYtYbYa[sm, mssm], GenAtAbAa[sm, mssm], GenGAUGE[sm], GenMSOFT[mssm]}];
    Close[ofs]];

End[];
EndPackage[];

