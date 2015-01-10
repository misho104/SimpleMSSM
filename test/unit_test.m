#!/bin/bash
# -*- Mathematica -*-
# MathematicaScript and MathKernel 9 have a bug that prevents output sent
# from a Mathematica script to stdout from being catched on a pipe or
# redirected to a file. This workaround runs the script inline into a
# MathKernel session.
# See http://mathematica.stackexchange.com/questions/20954/
ARG=`perl -e '$,=",";print map{qq/"$_"/}@ARGV' -- "$0" "$@"`
MathKernel -noinit -noprompt -run "\$MyCommandLine={$ARG}; $(sed '1,/^exit/d' $0) ; Exit[]"
exit $?

Needs["SLHA`", "../../SLHA.m"];
Needs["TreeSpec`", "../TreeSpec.m"];

UnitTest::Pass = "";
UnitTest::Fail = "";
UnitTest::Undetermined = "Assertion result not determined.";

Off[General::stop];

BeginUnitTest[module_] := ($AssertErrors = {module, {}});
EndUnitTest  [obj___]  := If[Length[$AssertErrors[[2]]] > 0,
                             (UnitTest::Fail = "\""<>$AssertErrors[[1]] <> "\" fails \"" <> # <> "\"."; Message[UnitTest::Fail]) &/@ $AssertErrors[[2]]; If[Length[List[obj]] > 0, Print[obj]]; Abort[],
                             (UnitTest::Pass = "\""<>$AssertErrors[[1]] <> "\" passed."); Message[UnitTest::Pass]
                            ];
Assert        [a_,     message_] := If[Not[a], AppendTo[$AssertErrors[[2]], message], True, Message[UnitTest::Undetermined]; AppendTo[$AssertErrors[[2]], "internal error"]];
AssertReal    [a_,     message_] := Assert[And @@ (Element[#, Reals] & /@ Flatten[{a}]), message];
AssertZero    [a_,     message_] := Module[{zero}, zero = {a}*0; Assert[Chop[{a}] == zero, message]];
AssertEqual   [a_, b_, message_] := Module[{zero}, zero = {a}*0; Assert[And[Dimensions[a] === Dimensions[b], Chop[{a}-{b}] == zero], message]];
AssertDiagonal[a_,     message_] := Module[{zero}, zero = {a}*0; Assert[{a - DiagonalMatrix[Diagonal[a]]} == zero, message]];
AssertUnitary [a_,     message_] := Module[{one},  one = IdentityMatrix[Length[a]]; Assert[And[Dimensions[one] === Dimensions[a], Chop[a.Conjugate[Transpose[a]] - one] == 0*one], message]];

RandomMatrix          [n_]:=RandomReal[{-1000,1000}, {n,n}];
RandomSymmetric       [n_]:=Module[{m = RandomMatrix[n]}, m+Transpose[m]];
RandomComplexMatrix   [n_]:=RandomComplex[{-1000-1000I, 1000+1000I}, {n,n}];
RandomComplexSymmetric[n_]:=RandomSymmetric[n] + I * RandomSymmetric[n];
RandomHermitian       [n_]:=Module[{m = RandomComplexMatrix[n]}, m+Transpose[Conjugate[m]]];

UNeutralinoRealDiagonalize[] := Module[
    { message = {},
      d, n,
      m = RandomSymmetric[4]
    },
    BeginUnitTest["NeutralinoRealDiagonalize"];

    {d, n} = NeutralinoRealDiagonalize[m];
    Assert[Dimensions[{d,n}] === {2,4,4}, "returns two 4x4 matrices"];
    (* assertion for masses *)
    AssertDiagonal[d, "returns a diagonal mass matrix."];
    AssertReal[d, "returns real masses."];
    Do[Assert[Abs[d[[i,i]]] <= Abs[d[[i+1,i+1]]], "returns ascending masses"], {i, 1, 3}];
    (* assertion for mixings *)
    AssertUnitary[n, "returns unitary NMIX"];
    AssertEqual[d, Conjugate[n].m.Conjugate[Transpose[n]], "returns NMIX in SLHA convention"];
    AssertReal[n, "returns real NMIX for real input"];
    EndUnitTest[m, d, n]];

UNeutralinoRealDiagonalizeComplex[] := Module[
    { message = {},
      d, n,
      m = RandomComplexSymmetric[4]
    },
    BeginUnitTest["NeutralinoRealDiagonalize(Complex)"];

    {d, n} = NeutralinoRealDiagonalize[m];
    Assert[Dimensions[{d,n}] === {2,4,4}, "returns two 4x4 matrices"];
    (* assertion for masses *)
    AssertDiagonal[d, "returns a diagonal mass matrix."];
    AssertReal[d, "returns real masses."];
    Do[Assert[Abs[d[[i,i]]] <= Abs[d[[i+1,i+1]]], "returns ascending masses"], {i, 1, 3}];
    (* assertion for mixings *)
    AssertUnitary[n, "returns unitary NMIX"];
    AssertEqual[d, Conjugate[n].m.Conjugate[Transpose[n]], "returns NMIX in SLHA convention"];
    EndUnitTest[m, d, n]];

UNeutralinoPositiveDiagonalize[] := Module[
    { message = {},
      d, n,
      m = RandomComplexSymmetric[4]
    },
    BeginUnitTest["NeutralinoPositiveDiagonalize"];

    {d, n} = NeutralinoPositiveDiagonalize[m];
    Assert[Dimensions[{d,n}] === {2,4,4}, "returns two 4x4 matrices"];
    (* assertion for masses *)
    AssertDiagonal[d, "returns a diagonal mass matrix."];
    Do[Assert[d[[i,i]] >= 0, "returns a positive-definite diagonal mass matrix."], {i, 1, 4}];
    Do[Assert[d[[i,i]] <= d[[i+1,i+1]], "returns ascending masses"], {i, 1, 3}];
    (* assertion for mixings *)
    AssertUnitary[n, "returns unitary NMIX"];
    AssertEqual[d, Conjugate[n].m.Conjugate[Transpose[n]], "returns NMIX in SLHA convention"];
    EndUnitTest[m, d, n]];

UCharginoDiagonalize[] := Module[
    { message = {},
      d, u, v,
      m = RandomMatrix[2]
    },
    BeginUnitTest["CharginoDiagonalize"];

    {d, u, v} = CharginoDiagonalize[m];
    Assert[Dimensions[{d,u,v}] === {3,2,2}, "returns three 2x2 matrices"];
    (* assertion for masses *)
    AssertDiagonal[d, "returns a diagonal mass matrix."];
    AssertReal[d, "returns real masses."];
    Assert[d[[1,1]] <= d[[2,2]], "returns ascending masses"];
    (* assertion for mixings *)
    AssertUnitary[u, "returns unitary UMIX"];
    AssertUnitary[v, "returns unitary VMIX"];
    AssertEqual[d, Conjugate[u].m.Conjugate[Transpose[v]], "returns UMIX and VMIX in SLHA convention"];
    AssertReal[u, "returns real UMIX for real input"];
    AssertReal[v, "returns real VMIX for real input"];
    EndUnitTest[m, d, u, v]];

USfermionDiagonalize[] := Module[
    { message = {},
      d, f,
      m = RandomSymmetric[2]
    },
    BeginUnitTest["SfermionDiagonalize"];

    {d, f} = SfermionDiagonalize[m];
    Assert[Dimensions[{d,f}] === {2,2,2}, "returns two 2x2 matrices"];
    (* assertion for masses *)
    AssertDiagonal[d, "returns a diagonal mass matrix."];
    AssertReal[d, "returns real masses."];
    Assert[d[[1,1]] <= d[[2,2]], "returns ascending masses"];
    (* assertion for mixings *)
    AssertUnitary[f, "returns unitary MIXING"];
    AssertEqual[d, f.m.Conjugate[Transpose[f]], "returns MIXING in SLHA convention"];
    AssertReal[f, "returns real MIXING for real input"];
    EndUnitTest[m, d, f]];

UNeutralinoRealDiagonalize[];
UNeutralinoRealDiagonalizeComplex[];
UNeutralinoPositiveDiagonalize[];
UCharginoDiagonalize[];
USfermionDiagonalize[];