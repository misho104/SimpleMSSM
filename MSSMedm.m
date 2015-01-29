(* ::Package:: *)

(* Time-Stamp: <2015-01-28 21:33:39 misho> *)

(* :Context: SimpleMSSM` *)

(* :Author: Sho Iwamoto *)

(* :Summary:
    Subpackage of SimpleMSSM: Return EDM values under specified mass spectrum.
    EDM formulae are obtained from 1406.0083 by T. Ibrahim, A. Itani, and P. Nath,
    and hep-ph/9906206 by S. Pokorski, J. Rosiek, and C. A. Savoy.
*)

(* :Copyright: 2015 Sho Iwamoto *)

(* :Package Version: 0.1
   Version 0.1 (beta version) hogehoge
*)

(* :Mathematica Version: 10.0 *)

(* :History: *)

(* :Keywords: *)

(* :Discussion: *)

(* :Warning: *)



SetDirectory[NotebookDirectory[]];
Needs["SLHA`", "vendor/slha-mathematica/SLHA.m"];
<<"MSSM.m";

BeginPackage["SimpleMSSM`"];

(* Usage messages *)
F := "\!\(\*StyleBox[\""<>#1<>"\", \""<>#2<>"\"]\)" &;
T := StringJoin[If[IntegerQ[#], F[ToString[#], "TR"], F[#, If[StringMatchQ[#,RegularExpression["[A-Za-z ]*"]], "TI", "TR"]]] & /@ List[##]]&;
B := "\!\(\*StyleBox[\""<>#<>"\", \"Bold\"]\)" &;
L := "\!\("<>T[#1]<>"\_"<>T[#2]<>"\)" &;
X := StringJoin[Riffle[List[##],", "]] &;

SUmassSQ::usage = "SUmassSQ["<>X[L["M","F"],L["M","L"],L["M","R"],T["\[Mu]"],T["A"],T["\[Beta]"],L["m","Z"],L["s","w"]]<>"] returns an up-type squark mass matrix with the given parameters, where "<>L["M","L"]<>" and "<>L["M","R"]<>" are squark soft masses, "<>L["M","F"]<>" is corresponding quark mass, and "<>T["A"]<>" is corresponding A-term.";
SDmassSQ::usage = "SDmassSQ["<>X[L["M","F"],L["M","L"],L["M","R"],T["\[Mu]"],T["A"],T["\[Beta]"],L["m","Z"],L["s","w"]]<>"] returns a down-type squark mass matrix with the given parameters, where "<>L["M","L"]<>" and "<>L["M","R"]<>" are squark soft masses, "<>L["M","F"]<>" is corresponding quark mass, and "<>T["A"]<>" is corresponding A-term.";

SLmassSQ::usage = "SLmassSQ["<>X[L["M","F"],L["M","L"],L["M","R"],T["\[Mu]"],T["A"],T["\[Beta]"],L["m","Z"],L["s","w"]]<>"] returns a charged slepton mass matrix with the given parameters, where "<>L["M","L"]<>" and "<>L["M","R"]<>" are slepton soft masses, "<>L["M","F"]<>" is corresponding lepton mass, and "<>T["A"]<>" is corresponding A-term.";

SNmassSQ::usage = "SNmassSQ["<>X[L["M","L"],T["\[Beta]"],L["m","Z"],L["s","w"]]<>"] returns sneutrino mass squared under the given parameters, where "<>L["M","L"]<>" is left-handed slepton soft masses. The Weinberg angle "<>L["s","w"]<>" is irrelevant to sneutrino mass but required for uniformity.";


Remove[SimpleMSSM`F, SimpleMSSM`T, SimpleMSSM`B, SimpleMSSM`L, SimpleMSSM`X];


(* Error messages *)

GetMatrices::NegativeSfermionMass = "Lighter sfermion mass becomes negative. The masses are: `1`, `2`.";

(* Other messages *)


Begin["`Private`"];


IfMissing=SLHA`IfMissing;


NewSM[slha_] := Module[{sm, ainv, gf, as, mz, mb, mt, mtau, ainvY0, ainvW0, ainvS0},
  {ainv, gf, as, mz, mb, mt, mtau} = slha["SMINPUTS"][#, IfMissing->"Abort"] & /@ {1,2,3,4,5,6,7};
  sm["v"]  = 246.;
  sm["sw"] = Sqrt[0.23116];
  sm["me"] = 0.510998928*^-3;
  sm["mm"] = 0.1056583715;
  sm["alpha0"] = 1/137.0356;

  {ainvY0, ainvW0, ainvS0} = {ainv * sm["cw"]^2, ainv * sm["sw"]^2, 1/as };
  sm["cw"] = Cos[ArcSin[sm["sw"]]];
  sm["tw"] = Tan[ArcSin[sm["sw"]]];
  sm["mZ"] = mz;
  sm["mW"] = mz * sm["cw"];
  sm["gY"] = Function[{Q}, Sqrt[ 4\[Pi] / (ainvY0 - (41/6) Log[Q/mz] / (2\[Pi]) )]];
  sm["g2"] = Function[{Q}, Sqrt[ 4\[Pi] / (ainvW0 + (19/6) Log[Q/mz] / (2\[Pi]) )]];
  sm["g3"] = Function[{Q}, Sqrt[ 4\[Pi] / (ainvS0 + (7)    Log[Q/mz] / (2\[Pi]) )]];
  sm["ye"] = Function[{Q}, Sqrt[2] * sm["me"] / sm["v"]]; (* calculated from pole mass; quantum correction ignored. *)
  sm];

NewMSSM[slha_] := Module[{mssm},
    mssm["M1"  ] = slha["MSOFT"][1, IfMissing->"Abort"];
    mssm["M2"]   = slha["MSOFT"][2, IfMissing->"Abort"];
    mssm["Ae",1] = slha["AE"][1, 1, IfMissing->0];
    mssm["mu"]   = slha["HMIX"][1, IfMissing->"Abort"];
    mssm["tb"]   = slha["HMIX"][2, IfMissing->"Abort"];
    mssm["beta"] = ArcTan[mssm["tb"]];
    mssm["LL",1] = slha["MSOFT"][31, IfMissing->"Abort"];
    mssm["ER",1] = slha["MSOFT"][34, IfMissing->"Abort"];
    Do[mssm[pid] = slha["MASS"][pid], {pid, Flatten[slha["MASS"]["keys"]]}];
    mssm];

NewPhase[] := Module[{ph},
    ph["B"]    = 0; (* theta_b *)
    ph["mu"]   = 0; (* theta_mu *)
    ph["vd"]   = 0; (* chi_1 *)
    ph["vu"]   = 0; (* chi_2 *)
    ph["Ae",1] = 0; (* xi_e *)
    ph["M1"]   = 0; (* xi_1 *)
    ph["M2"]   = 0; (* xi_2 *)
    ph["LL",1] = 0;
    ph["ER",1] = 0;
    ph[___]    = Missing[];
    ph]
(* Phases are (cf. hep-ph/9807501):
     B,   mu,   vd=v1, vu=v2, Ae, Ad, Au, Ma, sfermion-softmass;
     thB, thmu, chi1,  chi2,  xi_eud, xi_123,  IGNORED *)

NewParams[sm_, mssm_, phase_] := Module[{p},
    p[SM]    = sm;
    p[MSSM]  = mssm;
    p[PHASE] = phase;
    p[k__]  := With[
        {a = p[MSSM][k], b = p[PHASE][k]},
        If[AnyTrue[{a,b}, MatchQ[#, _Missing]&]===True, Missing[], a Exp[I b]]];
    p]


NeutralinoMass[p_]:=With[
    {m = NeutralinoMass[p["M1"], p["M2"], p[MSSM]["tb"], p["mu"], p[SM]["mZ"], p[SM]["sw"]],
     pd = Exp[-I p[PHASE]["vd"]], pu = Exp[-I p[PHASE]["vu"]]},
    {{m[[1,1]],   m[[1,2]],   m[[1,3]]pd, m[[1,4]]pu},
     {m[[2,1]],   m[[2,2]],   m[[2,3]]pd, m[[2,4]]pu},
     {m[[3,1]]pd, m[[3,2]]pd, m[[3,3]],   m[[3,4]]},
     {m[[4,1]]pu, m[[4,2]]pu, m[[4,3]],   m[[4,4]]}}]

CharginoMass[p_]:=With[
    {m = CharginoMass[p["M2"], p[MSSM]["tb"], p["mu"], p[SM]["mW"]]},
    {{m[[1,1]], m[[1,2]]Exp[-I p[PHASE]["vu"]]},
     {m[[2,1]]Exp[-I p[PHASE]["vd"]], m[[2,2]]}}];

SEmassSQ[p_]:=With[
    {m = SLmassSQ[p[SM]["me"], p["LL",1], p["ER",1], p["mu"], p["Ae",1], p[MSSM]["beta"], p[SM]["mZ"], p[SM]["sw"]],
     m12 = p[SM]["me"](Conjugate[p["Ae",1]Exp[I p[PHASE]["vd"]]] - p["mu"]Exp[I p[PHASE]["vu"]+I p[PHASE]["vd"]]p[MSSM]["tb"])},
    {{m[[1,1]],m12},{Conjugate[m12],m[[2,2]]}}];

SNmassSQ[p_]:=SNmassSQ[p["LL",1],p[MSSM]["beta"], p[SM]["mZ"], p[SM]["sw"]];


(* Notation in hep-ph/9807501 is: U = UMIX, V = VMIX, X = Dagger[NMIX], D = Dagger[SFERMIX].
   Note that MSSM.m always returns matrices in SLHA-convention. *)
GetMatrices[p_]:=Module[{m, tmp},
    m["N"]    = NeutralinoMass[p];
    m["C"]    = CharginoMass[p];
    m["L", 1] = SEmassSQ[p];
    m["V"]    = SNmassSQ[p];
    {m["n"], m["NMIX"]} = NeutralinoPositiveDiagonalize[m["N"]];
    {m["c"], m["UMIX"], m["VMIX"]} = CharginoDiagonalize[m["C"]];
    {m["l", 1], m["STAUMIX", 1]} = SfermionDiagonalize[m["L", 1]];
    m["n"] = Diagonal[Re[m["n"]]];
    m["c"] = Diagonal[Re[m["c"]]];
    m["l", 1] = Diagonal[Re[m["l", 1]]];
    m["l", 1] = If[AnyTrue[m["l", 1], #<0&], Message[GetMatrices::NegativeSfermionMass, Sequence@@(m["l",1])]; Abort[], Sqrt[m["l",1]]];
    m["v"] = Sqrt[m["V"]];
    m];

N2[x_]:=(1-x^2+2x Log[x])/(1-x)^3;
C2[x_]:=(3-4x+x^2+2Log[x])/(1-x)^3;
Fa[x_,y_]:= (C2[x]-C2[y])/(2(x-y));
Fb[x_,y_]:=-(N2[x]-N2[y])/(2(x-y));

(* EDM/e in GeV^{-1}. 1406.0083 by T. Ibrahim, A. Itani, and P. Nath.
   Note that the original reference, hep-ph/9807501, contains an error in the formulae
   even after the erratum in PRDx60. *)
ElectronEDM[p_]:=Module[
    {sm = p[SM], mssm = p[MSSM], ph = p[PHASE], m = GetMatrices[p],
     b = 3, Q = -1, T3 = -1/2, a0, b0, c0,
     diag, u, v, x, d, ke, Gammae, Etae, chargino, neutralino, prefactor,
     A = (C2[#]/2)&, B = (N2[#]/2)&},
    {a0, b0, c0} = {-Sqrt[2]sm["tw"](Q-T3), -Sqrt[2]T3, Sqrt[2]sm["tw"]Q};
    {u, v, x, d} = {m["UMIX"], m["VMIX"], Conjugate[Transpose[m["NMIX"]]], Conjugate[Transpose[m["STAUMIX",1]]]};
    ke=(sm["me"]Exp[-I ph["vd"]])/(Sqrt[2]sm["mW"]Cos[mssm["beta"]]);
    Gammae[i_] := ke Conjugate[u[[i,2]]v[[i,1]]];
    Etae[i_,k_] := (a0 x[[1,i]]Conjugate[d[[1,k]]]+b0 x[[2,i]]Conjugate[d[[1,k]]] - ke x[[b,i]]Conjugate[d[[2,k]]])(c0 x[[1,i]]d[[2,k]]-ke x[[b,i]]d[[1,k]]);
    prefactor = sm["alpha0"]/(4\[Pi] sm["sw"]^2);
    chargino[i_]      := prefactor m["c"][[i]]/(m["v"]^2) A[(m["c"][[i]]/m["v"])^2] Im[Gammae[i]];
    neutralino[i_,k_] := prefactor m["n"][[i]]/(m["l",1][[k]]^2) B[(m["n"][[i]]/m["l",1][[k]])^2] Im[Etae[i,k]] Q;
    {neutralino[#, 1]&/@{1,2,3,4}, neutralino[#, 2]&/@{1,2,3,4}, chargino[#]&/@{1,2}}]


(* EDM/e in GeV^{-1}. hep-ph/9906206 by S. Pokorski, J. Rosiek, and C. A. Savoy. *)
ElectronEDMins[p_]:=Module[
    {sm = p[SM], mssm = p[MSSM], ph = p[PHASE],
     M1, M2, mu, tb, mL, mE, mN,
     prefactor},
    (* Selectron masses are read from MSOFT because MASS can be in mass-eigenstates in SLHA2 convention. *)
    {M1, M2, mu, tb, mL, mE, mN} = {mssm["M1"], mssm["M2"], mssm["mu"], mssm["tb"], mssm["LL",1], mssm["ER",1], mssm[1000012]};
    prefactor = sm["me"] tb sm["alpha0"] / (4\[Pi] sm["sw"]^2);
    { prefactor            Fa[M2^2/mN^2,mu^2/mN^2] Im[p["M2"]p["mu"]] / mN^4,
      prefactor sm["tw"]^2 Fb[mL^2/M1^2,mE^2/M1^2] Im[p["M1"](p["mu"] - Conjugate[p["Ae",1]/tb])] / M1^4,
      prefactor sm["tw"]^2 Fb[M1^2/mL^2,mu^2/mL^2] Im[p["M1"]p["mu"]] / (2 mL^4),
     -prefactor            Fb[M2^2/mL^2,mu^2/mL^2] Im[p["M2"]p["mu"]] / (2 mL^4),
     -prefactor sm["tw"]^2 Fb[M1^2/mE^2,mu^2/mE^2] Im[p["M1"]p["mu"]] / mE^4}]


End[];
EndPackage[];
