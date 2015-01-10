(* ::Package:: *)

(* Time-Stamp: <2015-01-10 21:49:31 misho> *)

(* :Context: TreeSpec` *)

(* :Author: Sho Iwamoto *)

(* :Summary:
    Package to calculate MSSM mass spectrum at very treelevel.
*)

(* :Copyright: 2015 Sho Iwamoto *)

(* :Package Version: 1.0 $Revision: 0.1 $ *)

(* :Mathematica Version: 9.0 *)

(* :History: *)

(* :Keywords: *)

(* :Discussion: *)

(* :Warning: *)

BeginPackage["TreeSpec`"];

(* Usage messages *)

NeutralinoRealDiagonalize::usage = "\!\(\*RowBox[{\"NeutralinoRealDiagonalize\",\"[\",StyleBox[\"m\",\"TI\"],\"]\"}]\) diagonalizes a \*StyleBox[\"complex\",Bold] symmetric matrix \!\(\*RowBox[{StyleBox[\"m\",\"TI\"]}]\) by Automne-Takagi factorization, giving a list of matrices \!\(\*RowBox[{\"{\",RowBox[{StyleBox[\"d\",\"TI\"],\",\",StyleBox[\"N\",\"TI\"]}],\"}\"}]\), where \!\(\*RowBox[{StyleBox[\"d\",\"TI\"]}]\) is a \*StyleBox[\"real\",Bold] diagonal matrix with \!\(\*RowBox[{RowBox[{\"Abs\",\"[\",\(\*StyleBox[\"d\",\"TI\"]\_\*RowBox[{StyleBox[\"i\",\"TI\"],StyleBox[\"i\",\"TI\"]}]\),\"]\"}],\"<\",RowBox[{\"Abs\",\"[\",\(\*StyleBox[\"d\",\"TI\"]\_\*RowBox[{StyleBox[\"j\",\"TI\"],StyleBox[\"j\",\"TI\"]}]\),\"]\"}]}]\) for \!\(\*RowBox[{StyleBox[\"i\",\"TI\"],\"<\",StyleBox[\"j\",\"TI\"]}]\), and \!\(\*RowBox[{StyleBox[\"N\",\"TI\"]}]\) is the NMIX matrix in SLHA convention, i.e., \!\(\*RowBox[{StyleBox[\"d\",\"TI\"], \"\\[Equal]\", RowBox[{\(\*StyleBox[\"N\",\"TI\"]\^\"*\"\),StyleBox[\"m\",\"TI\"],\(\*StyleBox[\"N\",\"TI\"]\^\"\\[Dagger]\"\)}]}]\). The components of \!\(\*RowBox[{StyleBox[\"d\",\"TI\"]}]\) can be negative, while \!\(\*RowBox[{StyleBox[\"N\",\"TI\"]}]\) is real when \!\(\*RowBox[{StyleBox[\"m\",\"TI\"]}]\) is real.";

NeutralinoPositiveDiagonalize::usage = "\!\(\*RowBox[{\"NeutralinoPositiveDiagonalize\",\"[\",StyleBox[\"m\",\"TI\"],\"]\"}]\) diagonalizes a \*StyleBox[\"complex\",Bold] symmetric matrix \!\(\*RowBox[{StyleBox[\"m\",\"TI\"]}]\) by Automne-Takagi factorization, giving a list of matrices \!\(\*RowBox[{\"{\",RowBox[{StyleBox[\"d\",\"TI\"],\",\",StyleBox[\"N\",\"TI\"]}],\"}\"}]\), where \!\(\*RowBox[{StyleBox[\"d\",\"TI\"]}]\) is a \*StyleBox[\"positive-definite\",Bold] diagonal matrix with \!\(\*RowBox[{\(\*StyleBox[\"d\",\"TI\"]\_\*RowBox[{StyleBox[\"i\",\"TI\"],StyleBox[\"i\",\"TI\"]}]\),\"<\",\(\*StyleBox[\"d\",\"TI\"]\_\*RowBox[{StyleBox[\"j\",\"TI\"],StyleBox[\"j\",\"TI\"]}]\)}]\) for \!\(\*RowBox[{StyleBox[\"i\",\"TI\"],\"<\",StyleBox[\"j\",\"TI\"]}]\), and \!\(\*RowBox[{StyleBox[\"N\",\"TI\"]}]\) is the NMIX matrix in SLHA convention, i.e., \!\(\*RowBox[{StyleBox[\"d\",\"TI\"], \"\\[Equal]\", RowBox[{\(\*StyleBox[\"N\",\"TI\"]\^\"*\"\),StyleBox[\"m\",\"TI\"],\(\*StyleBox[\"N\",\"TI\"]\^\"\\[Dagger]\"\)}]}]\), but usually complex and cannot directly be written in SLHA file.";

CharginoDiagonalize::usage = "\!\(\*RowBox[{\"CharginoDiagonalize\",\"[\",StyleBox[\"m\",\"TI\"],\"]\"}]\) gives a list of matrices \!\(\*RowBox[{\"{\",RowBox[{StyleBox[\"d\",\"TI\"],\",\",StyleBox[\"U\",\"TI\"],\",\",StyleBox[\"V\",\"TI\"]}],\"}\"}]\) by singular value decomposition, where \!\(\*RowBox[{StyleBox[\"d\",\"TI\"]}]\) is a non-negative diagonal matrix with \!\(\*RowBox[{\(\*StyleBox[\"d\",\"TI\"]\_\*RowBox[{StyleBox[\"i\",\"TI\"],StyleBox[\"i\",\"TI\"]}]\),\"<\",\(\*StyleBox[\"d\",\"TI\"]\_\*RowBox[{StyleBox[\"j\",\"TI\"],StyleBox[\"j\",\"TI\"]}]\)}]\) for \!\(\*RowBox[{StyleBox[\"i\",\"TI\"],\"<\",StyleBox[\"j\",\"TI\"]}]\), and \!\(\*RowBox[{StyleBox[\"U\",\"TI\"]}]\) and \!\(\*RowBox[{StyleBox[\"V\",\"TI\"]}]\) are the UMIX and VMIX matrices in SLHA convention, i.e., \!\(\*RowBox[{StyleBox[\"d\",\"TI\"], \"\\[Equal]\", RowBox[{\(\*StyleBox[\"U\",\"TI\"]\^\"*\"\),StyleBox[\"m\",\"TI\"],\(\*StyleBox[\"V\",\"TI\"]\^\"\\[Dagger]\"\)}]}]\). Note that \!\(\*RowBox[{StyleBox[\"m\",\"TI\"]}]\) is assumed as \!\(\*RowBox[{\"\\[ScriptCapitalL]\",\"\\[ReverseElement]\",RowBox[{\((\(\\[Psi]\_\"-\"\))\^\"T\"\),StyleBox[\"m\", \"TI\"], \(\\[Psi]\_\"+\"\)}]}]\).";

SfermionDiagonalize::usage = "\!\(\*RowBox[{\"SfermionDiagonalize\",\"[\",StyleBox[\"m\",\"TI\"],\"]\"}]\) receives a two-by-two Hermitian matrix \!\(\*RowBox[{StyleBox[\"m\",\"TI\"]}]\) to yield a list of matrices \!\(\*RowBox[{\"{\",RowBox[{StyleBox[\"d\",\"TI\"],\",\",StyleBox[\"F\",\"TI\"]}],\"}\"}]\), where \!\(\*RowBox[{StyleBox[\"d\",\"TI\"]}]\) is a non-negative diagonal matrix with \!\(\*RowBox[{\(\*StyleBox[\"d\",\"TI\"]\_\"11\"\),\"<\",\(\*StyleBox[\"d\",\"TI\"]\_\"22\"\)}]\), and \!\(\*RowBox[{StyleBox[\"F\",\"TI\"]}]\) is the mixing matrix following SLHA convention, i.e., \!\(\*RowBox[{StyleBox[\"d\",\"TI\"], \"\\[Equal]\", RowBox[{StyleBox[\"F\",\"TI\"],StyleBox[\"m\",\"TI\"],\(\*StyleBox[\"F\",\"TI\"]\^\"\\[Dagger]\"\)}]}]\), or \!\(\*RowBox[{RowBox[{\"(\",GridBox[{{\(\"\\[Psi]\"\_1\)},{\(\"\\[Psi]\"\_2\)}}],\")\"}],\"\\[Equal]\",\"F\",RowBox[{\"(\",GridBox[{{\(\"\\[Psi]\"\_\"L\"\)},{\(\"\\[Psi]\"\_\"R\"\)}}],\")\"}]}]\).";

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
  {ainv, gf, as, mz, mb, mt, mtau} = SLHA`GetData[slha, "SMINPUTS", #] & /@ {1,2,3,4,5,6,7};

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

(* Neutralino mass diagonalization : N s.t. Diag = N^* M N^dagger *)
AutomneTakagi[m_] := Module[{v, p, x},
  v = Eigenvectors[Conjugate[Transpose[m]].m] // Orthogonalize;
  p = DiagonalMatrix[Exp[-I Arg[Diagonal[v.m.Transpose[v]]]/2]];
  x = Reverse[p.v];
  {x.m.Transpose[x], x}
]; (* Diag = x m x^T *)

NeutralinoPositiveDiagonalize[m_] := Module[{v, x},
  {v, x} = AutomneTakagi[m//N];
  {v, Conjugate[x]} // Chop];

NeutralinoRealDiagonalize[m_] := Module[{v, x},
  {v, x} = AutomneTakagi[m // N];
  (* Remove complex phase allowing negative mass *)
  Do[If[Re[x[[i,1]]]==Re[x[[i,2]]]==Re[x[[i,3]]]==Re[x[[i,4]]]==0, x[[i]] = x[[i]] * I], {i, 1, 4}];
  {x.m.Transpose[x], Conjugate[x] } // Chop];

(* Chargino mass diagonalization : (chi-)^T M (chi+),  (U, V) s.t. Diag = U^* M V^dagger *)
CharginoDiagonalize[m_] := Module[{u, d, v},
  {u, d, v} = SingularValueDecomposition[m // N]; (*Diag=u^dagger m v, but diag descending*)
  {u, v} = Transpose[Reverse[Transpose[#]]] & /@ {u, v}; (*Now Diag ascending*)
  u = Transpose[u];
  v = Transpose[Conjugate[v]];
  {Conjugate[u].m.Conjugate[Transpose[v]], u, v} // Chop];

(* Sfermion mass diagonalization : {{f1}, {f2}} = F_ij . {{fL}, {fR}}, or Diag = F M F^dagger *)
SfermionDiagonalize[m_] := Module[{a, b, f, diag},
  a = ArcTan[N[2 Abs[m[[1, 2]]]/(m[[1, 1]] - m[[2, 2]])]]/2;
  b = Arg[N[m[[1, 2]]]];
  f = { { Cos[a],           Sin[a] Exp[I b]},
        {-Sin[a] Exp[-I b], Cos[a]} };
  diag = f.m.Conjugate[Transpose[f]];
  If[Re[diag[[1, 1]] - diag[[2, 2]]] > 0, f = Reverse[f]];
  {f.m.Conjugate[Transpose[f]], f} // Chop]

(* Mass matrices *)
NeutralinoMass[sm_, mssm_] := Module[
    {sb = Sin[ArcTan[mssm["tb"]]], sw = sm["sw"],
     cb = Cos[ArcTan[mssm["tb"]]], cw = sm["cw"], mz = sm["mZ"]},
    { { mssm["M1"],          0,   -cb sw mz,    sb sw mz},
      {          0, mssm["M2"],    cb cw mz,   -sb cw mz},
      {  -cb sw mz,   cb cw mz,           0, -mssm["mu"]},
      {   sb sw mz,  -sb cw mz, -mssm["mu"],           0}}];

CharginoMass[sm_, mssm_] := Module[
    {sb = Sin[ArcTan[mssm["tb"]]], cb = Cos[ArcTan[mssm["tb"]]], mw = sm["mW"]},
    { {    mssm["M2"], Sqrt[2] sb mw },
      { Sqrt[2] cb mw,    mssm["mu"] }}];

SfermionMassSQ[mF_, mL_, mR_, DL_, DR_, m12_] := (* m12 = mF * Xterm, Xterm = Conjugate[A] - MuTan(Cot)Beta *)
    { {    mL^2 + mF^2 + DL, m12 },
      { Conjugate[m12], mR^2 + mF^2 + DR } };

Dterm[T3_, Q_, sm_, mssm_] := (T3 - Q sm["sw"]^2) Cos[2ArcTan[mssm["tb"]]] sm["mZ"]^2;

SUmassSQ[sm_, mssm_, flavor_, IgnoreDiagMF_:False] := Module[
    {DL = Dterm[1/2,  2/3, sm, mssm],
     DR = Dterm[ 0 , -2/3, sm, mssm],
     X  = If[flavor != 3, 0, Conjugate[mssm["At"]] - mssm["mu"]/mssm["tb"]],
     mF = If[flavor != 3, 0, sm["yt"][mssm["Q"]]*sm["v"]/Sqrt[2]] },
    SfermionMassSQ[If[IgnoreDiagMF,0,mF], mssm["QL", flavor], mssm["UR", flavor], DL, DR, X mF]];

SDmassSQ[sm_, mssm_, flavor_, IgnoreDiagMF_:False] := Module[
    {DL = Dterm[-1/2, -1/3, sm, mssm],
     DR = Dterm[  0 ,  1/3, sm, mssm],
     X  = If[flavor != 3, 0, Conjugate[mssm["Ab"]] - mssm["mu"]*mssm["tb"]],
     mF = If[flavor != 3, 0, sm["yb"][mssm["Q"]]*sm["v"]/Sqrt[2]] },
    SfermionMassSQ[If[IgnoreDiagMF,0,mF], mssm["QL", flavor], mssm["DR", flavor], DL, DR, X mF]];

SLmassSQ[sm_, mssm_, flavor_, IgnoreDiagMF_:False] := Module[
    {DL = Dterm[-1/2, -1, sm, mssm],
     DR = Dterm[  0 ,  1, sm, mssm],
     A  = Switch[flavor, 1, 0, 2, mssm["Am"], 3, mssm["Aa"], _, Abort[]],
     mF = Switch[flavor, 1, 0, 2, sm["ym"][mssm["Q"]], 3, sm["ya"][mssm["Q"]], _, Abort[]] * sm["v"]/Sqrt[2]
    },
    SfermionMassSQ[If[IgnoreDiagMF,0,mF], mssm["LL", flavor], mssm["ER", flavor], DL, DR, mF Conjugate[A - mssm["mu"] * mssm["tb"]]]];

SNmassSQ[sm_, mssm_, flavor_] := mssm["LL", flavor]^2 + Dterm[1/2, 0, sm, mssm];


(* MSSM *)
MSSMparameter[slha_, sm_] := Module[
    {p,
     mu, m1sq, m2sq, masq, tb, mA0, mHP,
     c, mZsq, sgnmu, musq, beta, b},
    p["M1"]   = Module[{k}, k = SLHA`Data[slha, "EXTPAR",  1]; If[k === Null, SLHA`GetData[slha, "MINPAR", 2], k]] // N;
    p["M2"]   = Module[{k}, k = SLHA`Data[slha, "EXTPAR",  2]; If[k === Null, SLHA`GetData[slha, "MINPAR", 2], k]] // N;
    p["M3"]   = Module[{k}, k = SLHA`Data[slha, "EXTPAR",  3]; If[k === Null, SLHA`GetData[slha, "MINPAR", 2], k]] // N;
    p["At"]   = Module[{k}, k = SLHA`Data[slha, "EXTPAR", 11]; If[k === Null, SLHA`GetData[slha, "MINPAR", 5], k]] // N;
    p["Ab"]   = Module[{k}, k = SLHA`Data[slha, "EXTPAR", 12]; If[k === Null, SLHA`GetData[slha, "MINPAR", 5], k]] // N;
    p["Aa"]   = Module[{k}, k = SLHA`Data[slha, "EXTPAR", 13]; If[k === Null, SLHA`GetData[slha, "MINPAR", 5], k]] // N;
    p["Am"]   = 0;
    Do[
        p["LL",i] = Module[{k}, k = SLHA`Data[slha, "EXTPAR", 30+i]; If[k === Null, SLHA`GetData[slha, "MINPAR", 1], k]] // N;
        p["ER",i] = Module[{k}, k = SLHA`Data[slha, "EXTPAR", 33+i]; If[k === Null, SLHA`GetData[slha, "MINPAR", 1], k]] // N;
        p["QL",i] = Module[{k}, k = SLHA`Data[slha, "EXTPAR", 40+i]; If[k === Null, SLHA`GetData[slha, "MINPAR", 1], k]] // N;
        p["UR",i] = Module[{k}, k = SLHA`Data[slha, "EXTPAR", 43+i]; If[k === Null, SLHA`GetData[slha, "MINPAR", 1], k]] // N;
        p["DR",i] = Module[{k}, k = SLHA`Data[slha, "EXTPAR", 46+i]; If[k === Null, SLHA`GetData[slha, "MINPAR", 1], k]] // N;
        ,  {i, 1, 3}];
    p["tb"]  = SLHA`GetData[slha, "MINPAR", 3] // N;
    {mu, masq, tb, mA0, mHP} = SLHA`Data[slha, "EXTPAR", #] & /@ {23, 24, 25, 26, 27} // N;
    c = Total[If[# === Null, 0, 1] & /@ {masq, mA0, mHP}];
    If[tb =!= Null, Message[MSSMparameter::tanb]; p["tb"] = tb];
    beta = ArcTan[p["tb"]];
    mZsq = sm["mZ"]^2;

    If[mu === Null,
       If[c != 0, Message[MSSMparameter::InvalidHiggsParameters]; Abort[]];
       m1sq = Module[{k}, k = SLHA`Data[slha, "EXTPAR", 21]; If[k === Null, SLHA`GetData[slha, "MINPAR", 1]^2, k]];
       m2sq = Module[{k}, k = SLHA`Data[slha, "EXTPAR", 22]; If[k === Null, SLHA`GetData[slha, "MINPAR", 1]^2, k]];
       (* Now m1sq and m2sq are set *)
       musq = (- m1sq - m2sq - mZsq + abs[(m1sq - m2sq)/Cos[2beta]]) / 2;
       If[musq < 0, Message[MSSMparameter::InvalidEWSB]; Abort[]];
       mu = Sign[GetSLHA`Data[slha, "MINPAR", 4]] * Sqrt[musq];
       masq = 2musq + m1sq + m2sq;
       b    = masq * Sin[2beta] / 2;
       If[masq < 0, Message[MSSMparameter::InvalidEWSBmasq]; Abort[]];
       ,
       If[c != 1, Message[MSSMparameter::InvalidHiggsParameters]; Abort[]];
       If[mA0 =!= Null, masq = mA0^2];
       If[mHP =!= Null, masq = mHP^2 - sm["mW"]^2];
       sgnmu = SLHA`Data[slha, "MINPAR", 4];
       If[And[sgnmu =!= Null, sgnmu * mu < 0], Message[MSSMparameter::InvalidMuSign]; Abort[]];
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

GenSPINFO[] :=
    {"SPINFO", {},
     {{"Block", "SPINFO"}, " Program Information"},
     {{1, "TreeSpec"},     " spectrum calculator"},
     {{2, "0.1"},          " version number"}};

GenMODSEL[] :=
    {"MODSEL", {},
     {{"Block", "MODSEL"}, " model selection"},
     {{1, 1}, " SUGRA-based general MSSM simulation"}};

GenMass[mssm_] :=
    {"MASS", {},
     {{"Block", "MASS"}, " mass spectrum"},
     {{     24, mssm[     24]}, " W"},
     {{     25, mssm[     25]}, " h (set by hand); tree level = " <> ToString[mssm[-25]]},
     {{     35, mssm[     35]}, " H0"},
     {{     36, mssm[     36]}, " A0"},
     {{     37, mssm[     37]}, " H+"},
     {{1000022, mssm[1000022]}, " N1"},
     {{1000023, mssm[1000023]}, " N2"},
     {{1000025, mssm[1000025]}, " N3"},
     {{1000035, mssm[1000035]}, " N4"},
     {{1000024, mssm[1000024]}, " C1"},
     {{1000037, mssm[1000037]}, " C2"},
     {{1000011, mssm[1000011]}, " ~e_L"},
     {{1000013, mssm[1000013]}, " ~mu_L"},
     {{2000011, mssm[2000011]}, " ~e_R"},
     {{2000013, mssm[2000013]}, " ~mu_R"},
     {{1000015, mssm[1000015]}, " ~tau_1"},
     {{2000015, mssm[2000015]}, " ~tau_2"},
     {{1000012, mssm[1000012]}, " ~nu_e"},
     {{1000014, mssm[1000014]}, " ~nu_mu"},
     {{1000016, mssm[1000016]}, " ~nu_tau"},
     {{1000021, mssm[1000021]}, " gluino (tree!)"},
     {{1000001, mssm[1000001]}, " ~d_L"},
     {{1000002, mssm[1000002]}, " ~u_L"},
     {{1000003, mssm[1000003]}, " ~s_L"},
     {{1000004, mssm[1000004]}, " ~c_L"},
     {{2000001, mssm[2000001]}, " ~d_R"},
     {{2000002, mssm[2000002]}, " ~u_R"},
     {{2000003, mssm[2000003]}, " ~s_R"},
     {{2000004, mssm[2000004]}, " ~c_R"},
     {{1000005, mssm[1000005]}, " ~b_1"},
     {{2000005, mssm[2000005]}, " ~b_2"},
     {{1000006, mssm[1000006]}, " ~t_1"},
     {{2000006, mssm[2000006]}, " ~t_2"}};

GenALPHA[mssm_] :=
    {"ALPHA", {},
     {{"Block", "ALPHA"}, " Effective Higgs mixing parameter"},
     {{mssm["a"]}, " alpha"}};

GenHMIX[mssm_] :=
    {"HMIX", {"Q" -> mssm["Q"]},
     {{"Block", "HMIX"}, " Higgs mixing parameters"},
     {{1, mssm["mu"]},   " mu (tree-level)"},
     {{2, mssm["tb"]},   " tb (tree-level)"},
     {{3, mssm["v"]},    " Higgs vev (input by hand)"},
     {{4, mssm["mAsq"]}, " mA^2 (tree-level)"}};

GenYAsub[name_, Q_, value_, comment_] :=
    {name, {"Q" -> Q},
     {{"Block", name}, ""},
(*   {{1, 1, 0.}, ""}, *)
(*   {{2, 2, 0.}, ""}, *)
     {{3, 3, value}, comment}};

GenYtYbYa[sm_, mssm_] :=
    { GenYAsub["YU", mssm["Q"], sm["yt"][mssm["Q"]]/Sin[ArcTan[mssm["tb"]]], " Yt (tree; 1loop for Pole->MSbar)"],
      GenYAsub["YD", mssm["Q"], sm["yb"][mssm["Q"]]/Cos[ArcTan[mssm["tb"]]], " Yb (tree; from msbar mass)"],
      GenYAsub["YE", mssm["Q"], sm["ya"][mssm["Q"]]/Cos[ArcTan[mssm["tb"]]], " Ytau (tree; from pole mass)"] };

GenAtAbAa[sm_, mssm_] :=
    { GenYAsub["AU", mssm["Q"], mssm["At"], " At (tree; i.e., = input value)"],
      GenYAsub["AD", mssm["Q"], mssm["Ab"], " Ab (tree; i.e., = input value)"],
      GenYAsub["AE", mssm["Q"], mssm["Aa"], " Atau (tree; i.e., = input value)"] };

GenMix[name_, matrix_, comment_] := Module[
    {k},
    k = MapIndexed[Function[{v, i}, MapIndexed[Function[{v2, j}, {Flatten[{i, j, v2}], ""}], v]], matrix]//Flatten[#,1] &;
    Join[{name, {}, {{"Block", name}, comment}}, k]];

GenGAUGE[sm_] :=
    {"GAUGE", {"Q"->sm["mZ"]},
     {{"Block", "GAUGE"}, " Gauge couplings (at mZ)"},
     {{1, sm["gY"][sm["mZ"]]}, " gY(mZ) MSbar"},
     {{2, sm["g2"][sm["mZ"]]}, " g2(mZ) MSbar"},
     {{3, sm["g3"][sm["mZ"]]}, " g3(mZ) MSbar"}};

GenMSOFT[mssm_] :=
    Join[{"MSOFT", {"Q" -> mssm["Q"]},
          {{"Block", "MSOFT"}, " MSSM SUSY-breaking (tree-level, i.e., = input)"},
          {{ 1, mssm["M1"  ]}, " M_1"},
          {{ 2, mssm["M2"  ]}, " M_2"},
          {{ 3, mssm["M3"  ]}, " M_3"},
          {{21, mssm["m1sq"]}, " mH1^2 = mHd^2"},
          {{22, mssm["m2sq"]}, " mH2^2 = mHu^2"}
         },
         {{30+#, mssm["LL",#]}, " m_LL["<>ToString[#]<>"]"} &/@{1,2,3},
         {{33+#, mssm["ER",#]}, " m_eR["<>ToString[#]<>"]"} &/@{1,2,3},
         {{40+#, mssm["QL",#]}, " m_QL["<>ToString[#]<>"]"} &/@{1,2,3},
         {{43+#, mssm["UR",#]}, " m_uR["<>ToString[#]<>"]"} &/@{1,2,3},
         {{46+#, mssm["DR",#]}, " m_dR["<>ToString[#]<>"]"} &/@{1,2,3}
        ];

CalculateSpectrum[slha_, OFile_] := Module[
    {sm, mssm,
     ofs,
     nmix, umix, vmix, stopmix, sbotmix, staumix, diag},
    sm = SMparameter[slha];

    (* MSSM parameters *)
    If[SLHA`Data[slha, "MODSEL", 1] =!= 1, Message[CalculateSpectrum::InvalidMODSEL]];
    If[SLHA`Data[slha, "EXTPAR", 0] =!= Null, Message[CalculateSpectrum::InputScaleIgnored]];
    mssm = MSSMparameter[slha, sm];
    mssm["Q"] = 200.; (* Simplicity *)

    (* masses *)
    {diag, nmix} = NeutralinoRealDiagonalize[NeutralinoMass[sm, mssm]];
    {mssm[1000022], mssm[1000023], mssm[1000025], mssm[1000035]} = Diagonal[diag];
    nmix = GenMix["NMIX", nmix, " Neutralino mixing matrix"];

    {diag, umix, vmix} = CharginoDiagonalize[CharginoMass[sm, mssm]];
    {mssm[1000024], mssm[1000037]} = Diagonal[diag];
    umix = GenMix["UMIX", umix, " Chargino U mixing matrix"];
    vmix = GenMix["VMIX", vmix, " Chargino V mixing matrix"];

    {diag, stopmix} = SfermionDiagonalize[SUmassSQ[sm, mssm, 3]];
    {mssm[1000006], mssm[2000006]} = Diagonal[diag] // Sqrt;
    {diag, sbotmix} = SfermionDiagonalize[SDmassSQ[sm, mssm, 3]];
    {mssm[1000005], mssm[2000005]} = Diagonal[diag] // Sqrt;
    {diag, staumix} = SfermionDiagonalize[SLmassSQ[sm, mssm, 3]];
    {mssm[1000015], mssm[2000015]} = Diagonal[diag] // Sqrt;
    stopmix = GenMix["STOPMIX", stopmix, " Stop mixing matrix"];
    sbotmix = GenMix["SBOTMIX", sbotmix, " Sbottom mixing matrix"];
    staumix = GenMix["STAUMIX", staumix, " Stau V mixing matrix"];

    {mssm[1000011], mssm[2000011]} = Diagonal[SLmassSQ[sm, mssm, 1]] // Sqrt; (* mixing is not considered (as required in MG5 mssm model) *)
    {mssm[1000013], mssm[2000013]} = Diagonal[SLmassSQ[sm, mssm, 2, True]] // Sqrt; (* Ignore diagonal fermion mass *)
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
    SLHA`WriteBlock         [ofs, #] & /@ {GenSPINFO[], SLHA`GetBlock[slha, "SMINPUTS"]};
    SLHA`WriteBlockAsComment[ofs, #] & /@ {GenMODSEL[], SLHA`GetBlock[slha, "MINPAR"], SLHA`GetBlock[slha, "EXTPAR"]};
    SLHA`WriteBlock         [ofs, #] & /@ {GenMass[mssm], GenALPHA[mssm], GenHMIX[mssm]};
    SLHA`WriteBlock         [ofs, #] & /@ {nmix, umix, vmix, stopmix, sbotmix, staumix};
    SLHA`WriteBlock         [ofs, #] & /@ GenYtYbYa[sm, mssm];
    SLHA`WriteBlock         [ofs, #] & /@ GenAtAbAa[sm, mssm];
    SLHA`WriteBlock         [ofs, #] & /@ {GenGAUGE[sm], GenMSOFT[mssm]};
    Close[ofs];
];

End[];
EndPackage[];

