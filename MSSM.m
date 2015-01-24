(* ::Package:: *)

(* Time-Stamp: <2015-01-24 13:28:04 misho> *)

(* :Context: SimpleMSSM` *)

(* :Author: Sho Iwamoto *)

(* :Summary:
    Tools related to the MSSM.
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



BeginPackage["SimpleMSSM`"];


(* Usage messages *)
F := "\!\(\*StyleBox[\""<>#1<>"\", \""<>#2<>"\"]\)" &;
T := StringJoin[If[IntegerQ[#], F[ToString[#], "TR"], F[#, If[StringMatchQ[#,RegularExpression["[A-Za-z ]*"]], "TI", "TR"]]] & /@ List[##]]&;
B := "\!\(\*StyleBox[\""<>#<>"\", \"Bold\"]\)" &;
L := "\!\("<>T[#1]<>"\_"<>T[#2]<>"\)" &;
X := StringJoin[Riffle[List[##],", "]] &;

NeutralinoRealDiagonalize::usage = "\!\(\*RowBox[{\"NeutralinoRealDiagonalize\",\"[\",StyleBox[\"m\",\"TI\"],\"]\"}]\) diagonalizes a \*StyleBox[\"complex\",Bold] symmetric matrix \!\(\*RowBox[{StyleBox[\"m\",\"TI\"]}]\) by Autonne-Takagi factorization, giving a list of matrices \!\(\*RowBox[{\"{\",RowBox[{StyleBox[\"d\",\"TI\"],\",\",StyleBox[\"N\",\"TI\"]}],\"}\"}]\), where \!\(\*RowBox[{StyleBox[\"d\",\"TI\"]}]\) is a \*StyleBox[\"real\",Bold] diagonal matrix with \!\(\*RowBox[{RowBox[{\"Abs\",\"[\",\(\*StyleBox[\"d\",\"TI\"]\_\*RowBox[{StyleBox[\"i\",\"TI\"],StyleBox[\"i\",\"TI\"]}]\),\"]\"}],\"<\",RowBox[{\"Abs\",\"[\",\(\*StyleBox[\"d\",\"TI\"]\_\*RowBox[{StyleBox[\"j\",\"TI\"],StyleBox[\"j\",\"TI\"]}]\),\"]\"}]}]\) for \!\(\*RowBox[{StyleBox[\"i\",\"TI\"],\"<\",StyleBox[\"j\",\"TI\"]}]\), and \!\(\*RowBox[{StyleBox[\"N\",\"TI\"]}]\) is the NMIX matrix in SLHA convention, i.e., \!\(\*RowBox[{StyleBox[\"d\",\"TI\"], \"\\[Equal]\", RowBox[{\(\*StyleBox[\"N\",\"TI\"]\^\"*\"\),StyleBox[\"m\",\"TI\"],\(\*StyleBox[\"N\",\"TI\"]\^\"\\[Dagger]\"\)}]}]\). The components of \!\(\*RowBox[{StyleBox[\"d\",\"TI\"]}]\) can be negative, while \!\(\*RowBox[{StyleBox[\"N\",\"TI\"]}]\) is real when \!\(\*RowBox[{StyleBox[\"m\",\"TI\"]}]\) is real.";

NeutralinoPositiveDiagonalize::usage = "\!\(\*RowBox[{\"NeutralinoPositiveDiagonalize\",\"[\",StyleBox[\"m\",\"TI\"],\"]\"}]\) diagonalizes a \*StyleBox[\"complex\",Bold] symmetric matrix \!\(\*RowBox[{StyleBox[\"m\",\"TI\"]}]\) by Autonne-Takagi factorization, giving a list of matrices \!\(\*RowBox[{\"{\",RowBox[{StyleBox[\"d\",\"TI\"],\",\",StyleBox[\"N\",\"TI\"]}],\"}\"}]\), where \!\(\*RowBox[{StyleBox[\"d\",\"TI\"]}]\) is a \*StyleBox[\"positive-definite\",Bold] diagonal matrix with \!\(\*RowBox[{\(\*StyleBox[\"d\",\"TI\"]\_\*RowBox[{StyleBox[\"i\",\"TI\"],StyleBox[\"i\",\"TI\"]}]\),\"<\",\(\*StyleBox[\"d\",\"TI\"]\_\*RowBox[{StyleBox[\"j\",\"TI\"],StyleBox[\"j\",\"TI\"]}]\)}]\) for \!\(\*RowBox[{StyleBox[\"i\",\"TI\"],\"<\",StyleBox[\"j\",\"TI\"]}]\), and \!\(\*RowBox[{StyleBox[\"N\",\"TI\"]}]\) is the NMIX matrix in SLHA convention, i.e., \!\(\*RowBox[{StyleBox[\"d\",\"TI\"], \"\\[Equal]\", RowBox[{\(\*StyleBox[\"N\",\"TI\"]\^\"*\"\),StyleBox[\"m\",\"TI\"],\(\*StyleBox[\"N\",\"TI\"]\^\"\\[Dagger]\"\)}]}]\), but usually complex and cannot directly be written in SLHA file.";

CharginoDiagonalize::usage = "\!\(\*RowBox[{\"CharginoDiagonalize\",\"[\",StyleBox[\"m\",\"TI\"],\"]\"}]\) gives a list of matrices \!\(\*RowBox[{\"{\",RowBox[{StyleBox[\"d\",\"TI\"],\",\",StyleBox[\"U\",\"TI\"],\",\",StyleBox[\"V\",\"TI\"]}],\"}\"}]\) by singular value decomposition, where \!\(\*RowBox[{StyleBox[\"d\",\"TI\"]}]\) is a non-negative diagonal matrix with \!\(\*RowBox[{\(\*StyleBox[\"d\",\"TI\"]\_\*RowBox[{StyleBox[\"i\",\"TI\"],StyleBox[\"i\",\"TI\"]}]\),\"<\",\(\*StyleBox[\"d\",\"TI\"]\_\*RowBox[{StyleBox[\"j\",\"TI\"],StyleBox[\"j\",\"TI\"]}]\)}]\) for \!\(\*RowBox[{StyleBox[\"i\",\"TI\"],\"<\",StyleBox[\"j\",\"TI\"]}]\), and \!\(\*RowBox[{StyleBox[\"U\",\"TI\"]}]\) and \!\(\*RowBox[{StyleBox[\"V\",\"TI\"]}]\) are the UMIX and VMIX matrices in SLHA convention, i.e., \!\(\*RowBox[{StyleBox[\"d\",\"TI\"], \"\\[Equal]\", RowBox[{\(\*StyleBox[\"U\",\"TI\"]\^\"*\"\),StyleBox[\"m\",\"TI\"],\(\*StyleBox[\"V\",\"TI\"]\^\"\\[Dagger]\"\)}]}]\). Note that \!\(\*RowBox[{StyleBox[\"m\",\"TI\"]}]\) is assumed as \!\(\*RowBox[{\"\\[ScriptCapitalL]\",\"\\[ReverseElement]\",RowBox[{\((\(\\[Psi]\_\"-\"\))\^\"T\"\),StyleBox[\"m\", \"TI\"], \(\\[Psi]\_\"+\"\)}]}]\).";

SfermionDiagonalize::usage = "\!\(\*RowBox[{\"SfermionDiagonalize\",\"[\",StyleBox[\"m\",\"TI\"],\"]\"}]\) receives a two-by-two Hermitian matrix \!\(\*RowBox[{StyleBox[\"m\",\"TI\"]}]\) to yield a list of matrices \!\(\*RowBox[{\"{\",RowBox[{StyleBox[\"d\",\"TI\"],\",\",StyleBox[\"F\",\"TI\"]}],\"}\"}]\), where \!\(\*RowBox[{StyleBox[\"d\",\"TI\"]}]\) is a non-negative diagonal matrix with \!\(\*RowBox[{\(\*StyleBox[\"d\",\"TI\"]\_\"11\"\),\"<\",\(\*StyleBox[\"d\",\"TI\"]\_\"22\"\)}]\), and \!\(\*RowBox[{StyleBox[\"F\",\"TI\"]}]\) is the mixing matrix following SLHA convention, i.e., \!\(\*RowBox[{StyleBox[\"d\",\"TI\"], \"\\[Equal]\", RowBox[{StyleBox[\"F\",\"TI\"],StyleBox[\"m\",\"TI\"],\(\*StyleBox[\"F\",\"TI\"]\^\"\\[Dagger]\"\)}]}]\), or \!\(\*RowBox[{RowBox[{\"(\",GridBox[{{\(\"\\[Psi]\"\_1\)},{\(\"\\[Psi]\"\_2\)}}],\")\"}],\"\\[Equal]\",\"F\",RowBox[{\"(\",GridBox[{{\(\"\\[Psi]\"\_\"L\"\)},{\(\"\\[Psi]\"\_\"R\"\)}}],\")\"}]}]\).";


NeutralinoMass::usage = "NeutralinoMass["<>L["M",1]<>", "<>L["M",2]<>", "<>T["tan","\[Beta]"]<>", "<>T["\[Mu]"]<>", "<>L["m","Z"]<>", "<>L["s","w"]<>"] returns a neutralino mass matrix with the given parameters.";

CharginoMass::usage = "CharginoMass["<>L["M",2]<>", "<>T["tan","\[Beta]"]<>", "<>T["\[Mu]"]<>", "<>L["m","W"]<>"] returns a chargino mass matrix with the given parameters.";


SUmassSQ::usage = "SUmassSQ["<>X[L["M","F"],L["M","L"],L["M","R"],T["\[Mu]"],T["A"],T["\[Beta]"],L["m","Z"],L["s","w"]]<>"] returns an up-type squark mass matrix with the given parameters, where "<>L["M","L"]<>" and "<>L["M","R"]<>" are squark soft masses, "<>L["M","F"]<>" is corresponding quark mass, and "<>T["A"]<>" is corresponding A-term.";
SDmassSQ::usage = "SDmassSQ["<>X[L["M","F"],L["M","L"],L["M","R"],T["\[Mu]"],T["A"],T["\[Beta]"],L["m","Z"],L["s","w"]]<>"] returns a down-type squark mass matrix with the given parameters, where "<>L["M","L"]<>" and "<>L["M","R"]<>" are squark soft masses, "<>L["M","F"]<>" is corresponding quark mass, and "<>T["A"]<>" is corresponding A-term.";

SLmassSQ::usage = "SLmassSQ["<>X[L["M","F"],L["M","L"],L["M","R"],T["\[Mu]"],T["A"],T["\[Beta]"],L["m","Z"],L["s","w"]]<>"] returns a charged slepton mass matrix with the given parameters, where "<>L["M","L"]<>" and "<>L["M","R"]<>" are slepton soft masses, "<>L["M","F"]<>" is corresponding lepton mass, and "<>T["A"]<>" is corresponding A-term.";

SNmassSQ::usage = "SNmassSQ["<>X[L["M","L"],T["\[Beta]"],L["m","Z"],L["s","w"]]<>"] returns sneutrino mass squared under the given parameters, where "<>L["M","L"]<>" is left-handed slepton soft masses. The Weinberg angle "<>L["s","w"]<>" is irrelevant to sneutrino mass but required for uniformity.";


Remove[SimpleMSSM`F, SimpleMSSM`T, SimpleMSSM`B, SimpleMSSM`L, SimpleMSSM`X];


(* Error messages *)

(* Other messages *)


Begin["`Private`"];


(* Neutralino mass diagonalization : N s.t. Diag = N^* M N^dagger *)
AutonneTakagi[m_] := Module[{v, p, x},
  v = Eigenvectors[Conjugate[Transpose[m]].m] // Orthogonalize;
  p = DiagonalMatrix[Exp[-I Arg[Diagonal[v.m.Transpose[v]]]/2]];
  x = Reverse[p.v];
  {x.m.Transpose[x], x}
]; (* Diag = x m x^T *)

NeutralinoPositiveDiagonalize[m_] := Module[{v, x},
  {v, x} = AutonneTakagi[m//N];
  {v, Conjugate[x]} // Chop];

NeutralinoRealDiagonalize[m_] := Module[{v, x},
  {v, x} = AutonneTakagi[m // N];
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

NeutralinoMass[M1_, M2_, tb_, mu_, mz_, sw_] := With[
    {sb = Sin[ArcTan[tb]], cb = Cos[ArcTan[tb]], cw = Sqrt[1-sw^2]},
    {{        M1,         0, -cb sw mz,  sb sw mz },
     {         0,        M2,  cb cw mz, -sb cw mz },
     { -cb sw mz,  cb cw mz,         0,       -mu },
     {  sb sw mz, -sb cw mz,       -mu,         0 }}];

CharginoMass[M2_, tb_, mu_, mw_] := With[
    {sb = Sin[ArcTan[tb]], cb = Cos[ArcTan[tb]]},
    {{            M2, Sqrt[2] sb mw },
     { Sqrt[2] cb mw,            mu }}];

Dterm[T3_, Q_, beta_, mZ_, sw_] := (T3 - Q sw^2) Cos[2beta] mZ^2;

SfermionMassSQ[mF_, mL_, mR_, mu_, A_, DL_, DR_, tbcb_] := With[
    {m12 = mF (Conjugate[A] - mu tbcb) },
    {{ mL^2 + mF^2 + DL,              m12 },
     {   Conjugate[m12], mR^2 + mF^2 + DR }}];

SUmassSQ[mF_, mL_, mR_, mu_, A_, beta_, mZ_, sw_] :=
    SfermionMassSQ[mF, mL, mR, mu, A,
                   Dterm[1/2,  2/3, beta, mZ, sw],
                   Dterm[  0, -2/3, beta, mZ, sw],
                   Cot[beta]];

SDmassSQ[mF_, mL_, mR_, mu_, A_, beta_, mZ_, sw_] :=
    SfermionMassSQ[mF, mL, mR, mu, A,
                   Dterm[-1/2, -1/3, beta, mZ, sw],
                   Dterm[   0,  1/3, beta, mZ, sw],
                   Tan[beta]];

SLmassSQ[mF_, mL_, mR_, mu_, A_, beta_, mZ_, sw_] :=
    SfermionMassSQ[mF, mL, mR, mu, A,
                   Dterm[-1/2, -1, beta, mZ, sw],
                   Dterm[   0,  1, beta, mZ, sw],
                   Tan[beta]];

SNmassSQ[mL_, beta_, mZ_, sw_] := mL^2 + Dterm[1/2, 0, beta, mZ, sw];


End[];
EndPackage[];

