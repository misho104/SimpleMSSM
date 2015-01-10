(* ::Package:: *)

(* Time-Stamp: <2015-01-10 19:34:16 misho> *)

(* :Context: MSSMgm2` *)

(* :Author: Sho Iwamoto *)

(* :Summary:
    Package to calculate MSSM muon g-2 at one-loop level.
    Notation follows my dissertation.
*)

(* :Copyright: 2012-2015 Sho Iwamoto *)

(* :Package Version: 2.0 $Revision: 0.1 $ *)

(* :Mathematica Version: 9.0 *)

(* :History:
   Version 1.0 [Oct 2013, SI] Initial version.
   Version 2.0 [Jan 2015, SI] Connected to SLHA module.
*)

(* :Keywords: *)

(* :Discussion: *)

(* :Warning: *)

BeginPackage["MSSMgm2`"];

(* Usage messages *)

a\[Mu]MSSM::usage = "a\[Mu]MSSM[params] or a\[Mu]MSSM[M1, M2, tb, \[Mu], MsmuL, MsmuR, Msnumu] returns the list of the SUSY contribution to the muon g-2 calculated with mass eigenstates.";

a\[Mu]MSSMt::usage = "a\[Mu]MSSMt[params] returns total value of the SUSY contribution to the muon g-2 calculated with mass eigenstates.";

a\[Mu]MSSMi::usage = "a\[Mu]MSSMi[params] or a\[Mu]MSSMi[M1, M2, tb, \[Mu], MsmuL, MsmuR, Msnumu] returns the list of the SUSY contribution to the muon g-2 calculated with mass insertion method.";

a\[Mu]MSSMit::usage = "a\[Mu]MSSMit[params] returns total value of the SUSY contribution to the muon g-2 calculated with mass insertion method.";

values::usage = "values have Standard Model parameters.";

PrintMassMatrices::usage = "PrintMassMatrices[params] shows mass matrices for given params.";

a\[Mu]MSSMiSLHA::usage  = "SLHA version of a\[Mu]MSSMi";

a\[Mu]MSSMitSLHA::usage = "SLHA version of a\[Mu]MSSMit";

a\[Mu]MSSMSLHA::usage   = "SLHA version of a\[Mu]MSSM";

(* Error messages *)

GetSoftParams::SLHAnotloaded = "SLHA module not loaded.";

(* Other messages *)


Protect[m\[Mu], MZ, MW, \[Theta]W, gY, gW, y\[Mu]];
Protect[M1, M2, tb, \[Mu], MsmuL, MsmuR, A\[Mu], Msnumu];

Begin["`Private`"];

values={
  m\[Mu] ->  0.105658367,
  MZ -> 91.1876,
  MW -> 80.399,
  \[Theta]W -> ArcSin[Sqrt[0.23116]],
  gY -> 0.358,
  gW -> 0.651,
  y\[Mu] -> 0.105658367/(246.67/Sqrt[2])
};



(* ::Text:: *)
(*Mass Matrices*)


Module[{s\[Beta], c\[Beta]},
  {s\[Beta], c\[Beta]} = {tb/Sqrt[1+tb^2],1/Sqrt[1+tb^2]};
  matrix$smu = {
    { MsmuL^2, m\[Mu] ( A\[Mu] / ( y\[Mu] / c\[Beta] ) - \[Mu] tb ) },
    {          m\[Mu] ( A\[Mu] / ( y\[Mu] / c\[Beta] ) - \[Mu] tb ), MsmuR^2 }
  };
  matrix$neut = {
    {            M1,              0, -MZ c\[Beta] Sin[\[Theta]W],  MZ s\[Beta] Sin[\[Theta]W]},
    {             0,             M2,  MZ c\[Beta] Cos[\[Theta]W], -MZ s\[Beta] Cos[\[Theta]W]},
    {-MZ c\[Beta] Sin[\[Theta]W],  MZ c\[Beta] Cos[\[Theta]W],              0,             -\[Mu]},
    { MZ s\[Beta] Sin[\[Theta]W], -MZ s\[Beta] Cos[\[Theta]W],             -\[Mu],              0}
  };
  matrix$char = {
    {           M2, Sqrt[2] MW s\[Beta]},
    {Sqrt[2] MW c\[Beta],             \[Mu]}
  };
];


(* ::Text:: *)
(*Mass Diagonalization*)


mass$SMU[params_]:=Module[{mass,R},
  {mass,R}=Eigensystem[matrix$smu//.params//.values];
  R=Orthogonalize[R];
  If[mass[[#]]<0,Abort[]]&/@{1,2};
  {Sqrt[mass[[#]]],R[[#,1]],R[[#,2]]}&/@{1,2}
]
mass$BWH[params_]:=Module[{mass,R,ss},
  {mass,R}=Eigensystem[matrix$neut//.params//.values];
  R=Orthogonalize[R];
  (R[[#,All]]=R[[#,All]]*If[mass[[#]]<0,I,1])&/@{1,2,3,4};
  ss={Abs[mass[[#]]],R[[#,1]],R[[#,2]],R[[#,3]]}&/@{1,2,3,4};
  Conjugate[ss]
]
mass$WH[params_]:=Module[{U,W,V},
  {U,W,V}=SingularValueDecomposition[matrix$char//.params//.values];
  {U,V}=Orthogonalize/@{U,V};
  {W[[#,#]],V[[1,#]],Conjugate[U[[#,2]]]}&/@{1,2}
]

PrintMassMatrices[params_] := Module[{},
  Print[matrix$smu//.params//.values];
  Print[matrix$neut//.params//.values];
  Print[matrix$char//.params//.values];
]

(*Full Method*)

N1[x_]:=(1-6x+3x^2+2x^3-6x^2 Log[x])/(1-x)^4
N2[x_]:=(1-x^2+2x Log[x])/(1-x)^3
C1[x_]:=(2+3x-6x^2+x^3+6x Log[x])/(1-x)^4
C2[x_]:=(3-4x+x^2+2Log[x])/(1-x)^3
a\[Mu]NEUTsub[MMU_, MCHI_, GL_, GR_]:=1/(16\[Pi]^2) (-Re[GL Conjugate[GR]] (m\[Mu] MCHI)/MMU^2 N2[MCHI^2/MMU^2]-(Abs[GL]^2+Abs[GR]^2)/6 m\[Mu]^2/MMU^2 N1[MCHI^2/MMU^2])/.values
a\[Mu]CHARsub[MMU_, MCHI_, GL_, GR_]:=1/(16\[Pi]^2) (-Re[GL Conjugate[GR]] (m\[Mu] MCHI)/MMU^2 C2[MCHI^2/MMU^2]+(Abs[GL]^2+Abs[GR]^2)/6 m\[Mu]^2/MMU^2 C1[MCHI^2/MMU^2])/.values
a\[Mu]NN[y\[Mu]MSSM_, Msmu_, Mneut_, bino_, wino_, higgsino_, \[Mu]L_, \[Mu]R_]:=a\[Mu]NEUTsub[Msmu,Mneut,{Conjugate[(gY bino +gW wino)/Sqrt[2]], -Conjugate[higgsino] y\[Mu]MSSM}.{\[Mu]L,\[Mu]R},{-higgsino y\[Mu]MSSM,-Sqrt[2]gY bino}.{\[Mu]L,\[Mu]R}]
a\[Mu]CN[y\[Mu]MSSM_, Msnu_, Mchar_, wino_, higgsino_]:=a\[Mu]CHARsub[Msnu,Mchar,-gW wino,y\[Mu]MSSM Conjugate[higgsino]]

a\[Mu]MSSM[$M1_, $M2_, $tb_, $\[Mu]_, $smuL_, $smuR_, $snumu_, $A\[Mu]_] := a\[Mu]MSSM[{M1->$M1, M2->$M2, tb->$tb, \[Mu]->$\[Mu], MsmuL->$smuL, MsmuR->$smuR, Msnumu->$snumu, A\[Mu]->$A\[Mu]}]
a\[Mu]MSSM[$M1_, $M2_, $tb_, $\[Mu]_, $smuL_, $smuR_, $snumu_]       := a\[Mu]MSSM[{M1->$M1, M2->$M2, tb->$tb, \[Mu]->$\[Mu], MsmuL->$smuL, MsmuR->$smuR, Msnumu->$snumu, A\[Mu]->0}]
a\[Mu]MSSM[params_]:= Module[{y, m, n, c},
  y = y\[Mu] / (1/Sqrt[1+tb^2])//.params/.values;
  m = mass$SMU[params];
  n = mass$BWH[params];
  c = mass$WH[params];
  {a\[Mu]NN[y, m[[1,1]], #[[1]], #[[2]], #[[3]] ,#[[4]], m[[1,2]],m [[1,3]]]&/@n,
   a\[Mu]NN[y, m[[2,1]], #[[1]], #[[2]], #[[3]] ,#[[4]], m[[2,2]],m [[2,3]]]&/@n,
   a\[Mu]CN[y, Msnumu//.params, #[[1]], #[[2]], #[[3]]]&/@c
  }
]
a\[Mu]MSSMt[params__] := Total/@a\[Mu]MSSM[params]//Total;


(* ::Text:: *)
(*Mass Insertion Method*)


Fa[x_,y_]:= (C2[x]-C2[y])/(2(x-y))
Fb[x_,y_]:=-(N2[x]-N2[y])/(2(x-y))
F1[M1_,M2_,tb_,\[Mu]_,smuL_,smuR_,snumu_]:=  gW^2/( 8\[Pi]^2) (m\[Mu]^2 M2 \[Mu] tb)/snumu^4 Fa[  M2^2/snumu^2,    \[Mu]^2/snumu^2]//.values
F2[M1_,M2_,tb_,\[Mu]_,smuL_,smuR_,snumu_]:=  gY^2/( 8\[Pi]^2) (m\[Mu]^2    \[Mu] tb)/M1^3    Fb[smuL^2/   M1^2, smuR^2/   M1^2]//.values
F3[M1_,M2_,tb_,\[Mu]_,smuL_,smuR_,snumu_]:=  gY^2/(16\[Pi]^2) (m\[Mu]^2 M1 \[Mu] tb)/smuL^4  Fb[  M1^2/ smuL^2,    \[Mu]^2/ smuL^2]//.values
F4[M1_,M2_,tb_,\[Mu]_,smuL_,smuR_,snumu_]:= -gW^2/(16\[Pi]^2) (m\[Mu]^2 M2 \[Mu] tb)/smuL^4  Fb[  M2^2/ smuL^2,    \[Mu]^2/ smuL^2]//.values
F5[M1_,M2_,tb_,\[Mu]_,smuL_,smuR_,snumu_]:= -gY^2/( 8\[Pi]^2) (m\[Mu]^2 M1 \[Mu] tb)/smuR^4  Fb[  M1^2/ smuR^2,    \[Mu]^2/ smuR^2]//.values
a\[Mu]MSSMi[params_]:=a\[Mu]MSSMi[Sequence@@({M1,M2,tb,\[Mu],MsmuL,MsmuR,Msnumu}//.params)]
a\[Mu]MSSMit[params_]:=Total[a\[Mu]MSSMi[params]]
a\[Mu]MSSMi[params__]:=#[params]&/@{F1,F2,F3,F4,F5}
a\[Mu]MSSMit[params__]:=Total[a\[Mu]MSSMi[params]]


(* ::Text:: *)
(*With SLHA package*)


GetSoftParams[slha_]:=If[MemberQ[$ContextPath,"SLHA`"],
  { M1->GetData[slha,"MSOFT",1],M2->GetData[slha,"MSOFT",2],tb->GetData[slha,"HMIX",2],\[Mu]->GetData[slha,"HMIX",1],
    MsmuL->GetData[slha,"MASS",1000013],MsmuR->GetData[slha,"MASS",2000013],Msnumu->GetData[slha,"MASS",1000014],A\[Mu]->Data[slha,"AE",2,2]},
  Message[GetSoftParams::SLHAnotloaded];Abort[]];
a\[Mu]MSSMiSLHA[slha_]:=a\[Mu]MSSMi[GetSoftParams[slha]];
a\[Mu]MSSMitSLHA[slha_]:=a\[Mu]MSSMit[GetSoftParams[slha]];
a\[Mu]MSSMSLHA[slha_]:=a\[Mu]MSSM[Replace[GetSoftParams[slha],(A\[Mu]->Null)->(A\[Mu]->0),{1}]];


(* ::Text:: *)
(*Cho-Hagiwara-Matsumoto-Nomura Benchmark*)


(* ::Input:: *)
(*{t,u,m,e}={10,200,150,300};params=Sequence[0.5m,m,t,u,e+0.1,e,e-0.1];{a\[Mu]MSSMt[params]*10^10,a\[Mu]MSSMit[params]*10^10,25.0}*)
(*{t,u,m,e}={10,200,450,120};params=Sequence[0.5m,m,t,u,e+0.1,e,e-0.1];{a\[Mu]MSSMt[params]*10^10,a\[Mu]MSSMit[params]*10^10,25.9}*)
(*{t,u,m,e}={10,800,150,200};params=Sequence[0.5m,m,t,u,e+0.1,e,e-0.1];{a\[Mu]MSSMt[params]*10^10,a\[Mu]MSSMit[params]*10^10,27.1}*)
(*{t,u,m,e}={10,800,500,150};params=Sequence[0.5m,m,t,u,e+0.1,e,e-0.1];{a\[Mu]MSSMt[params]*10^10,a\[Mu]MSSMit[params]*10^10,24.7}*)
(*{t,u,m,e}={50,800,150,550};params=Sequence[0.5m,m,t,u,e+0.1,e,e-0.1];{a\[Mu]MSSMt[params]*10^10,a\[Mu]MSSMit[params]*10^10,26.3}*)
(*{t,u,m,e}={50,800,900,280};params=Sequence[0.5m,m,t,u,e+0.1,e,e-0.1];{a\[Mu]MSSMt[params]*10^10,a\[Mu]MSSMit[params]*10^10,27.7}*)
(**)


End[];
EndPackage[];
