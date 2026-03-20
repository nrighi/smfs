(* ::Package:: *)

(* ::Input:: *)
(**)


BeginPackage["SMFs`"];


(* ::Input::Initialization:: *)
PrecomputedRules::usage="PrecomputedRules[n] returns replacement rules that map the symbolic modular forms E4, E6, E10, E12, \[Chi]10, \[Chi]12 and their derivatives in variables T,Z,U to truncated Fourier expansions of order n. Useful for numerical and plotting purposes.";

Ek::usage="Ek[k, n, q1, r, q2] gives the truncated genus-2 Eisenstein series of weight k up to order n in q_1=e^{2\[Pi]iT} and q_2=e^{2\[Pi]iU}.";
Chi10::usage="Chi10[n, q1, r, q2] gives the truncated Fourier expansion of the Igusa cusp form \[Chi]_10 up to order n in the nome variables q_1=e^{2\[Pi]iT} and q_2=e^{2\[Pi]iU}.";
Chi12::usage="Chi12[n, q1, r, q2] gives the truncated Fourier expansion of the Igusa cusp form \[Chi]_12 up to order n in the nome variables q_1=e^{2\[Pi]iT} and q_2=e^{2\[Pi]iU}.";

E4::usage="E4[T,Z,U] denotes the genus-2 Siegel Eisenstein series of weight 4.";
E6::usage="E6[T,Z,U] denotes the genus-2 Siegel Eisenstein series of weight 6.";
E10::usage="E10[T,Z,U] denotes the genus-2 Siegel Eisenstein series of weight 10.";
E12::usage="E12[T,Z,U] denotes the genus-2 Siegel Eisenstein series of weight 12.";
\[Chi]10::usage="\[Chi]10[T,Z,U] denotes the symbolic Igusa cusp form of weight 10 with its truncated Fourier expansion after substituting the nomes q_1=e^{2\[Pi]iT}, r=e^{2\[Pi]iZ}, q_2=e^{2\[Pi]iU}. Apply PrecomputedRules[n] to replace it by its truncated Fourier expansion.";
\[Chi]12::usage="\[Chi]12[T,Z,U] denotes the symbolic Igusa cusp form of weight 12 with its truncated Fourier expansion after substituting the nomes q_1=e^{2\[Pi]iT}, r=e^{2\[Pi]iZ}, q_2=e^{2\[Pi]iU}. Apply PrecomputedRules[n] to replace it by its truncated Fourier expansion.";


Begin["`Private`"];


ClearAll[fundDiscDecompose,bernoulliDirichletB,cohenC,alphaTerm,fourierCoeff,SiegelEisensteinGenus2Expansion];

(*Decompose D=D0 f^2 with D0 a fundamental discriminant (D\[LessEqual]0)*)
fundDiscDecompose[D_Integer]:=Module[{abs=Abs[D],mAbs,m,D0,f},If[D==0,Return[{0,1}]];
mAbs=Times@@(First/@Select[FactorInteger[abs],OddQ[Last[#]]&]);
m=Sign[D]*mAbs;
(*squarefree up to the 4-factor rule*)D0=If[Mod[m,4]==1,m,4 m];
(*fundamental discriminant*)f=IntegerPart[Sqrt[abs/Abs[D0]]];
(*since D=D0 f^2*){D0,f}];
(*Generalized Bernoulli B_{k,\[Chi]_D0};\[Chi]_D0(n)=KroneckerSymbol[D0,n]*)
bernoulliDirichletB[D0_Integer,k_Integer]:=Module[{f=Abs[D0]},f^(k-1) Sum[KroneckerSymbol[D0,a]*BernoulliB[k,a/f],{a,1,f}]];

(*Cohen function C(s-1,D) with s=w (weights 4,6,10,12)*)
cohenC[w_Integer,D_Integer]:=Module[{D0,f,Lval},If[D==0,Return[1]];
(*conventionally not used; \[Alpha](0) handles this*){D0,f}=fundDiscDecompose[D];
(*L_{D0}(2-w) at negative integer via generalized Bernoulli numbers:L(1-k,\[Chi])=-B_{k,\[Chi]}/k with k=w-1*)Lval=-bernoulliDirichletB[D0,w-1]/(w-1);
Lval*Sum[MoebiusMu[d]*KroneckerSymbol[D0,d]*d^(w-2)*DivisorSigma[2 w-3,Quotient[f,d]],{d,Divisors[f]}]];

(*\[Alpha](D)=(1/\[Zeta](3-2w))*C(w-1,D)*)
alphaTerm[D_Integer,w_Integer]:=alphaTerm[D,w]=Module[{},If[D==0,1,(1/Zeta[3-2 w])*cohenC[w,D]]];

(*Fourier coefficient a_w(T) for T=[[a,b/2],[b/2,c]]*)
fourierCoeff[a_Integer,b_Integer,c_Integer,w_Integer]:=Module[{g,D,s},If[a==0&&b==0&&c==0,Return[1]];
If[a<0||c<0||b^2>4 a c,Return[0]];
g=GCD[a,b,c];
D=b^2-4 a c;
s=Sum[d^(w-1)*alphaTerm[Quotient[D,d^2],w],{d,Divisors[g]}];
-(2 w/BernoulliB[w])*s];

(* \[CapitalSigma] a_w(a,b,c) q1^a q2^c r^b with 0\[LessEqual]a,c\[LessEqual]nMax and|b|bounded by b^2\[LessEqual]4ac *)
SiegelEisensteinGenus2Expansion[weight_Integer,nMax_Integer,q1_,r_,q2_]:=Module[{w=weight,terms},If[!MemberQ[{4,6,10,12},w],Return[$Failed]];
terms=Flatten@Table[If[b^2<=4 a c,With[{coef=fourierCoeff[a,b,c,w]},If[coef===0,Nothing,coef*q1^a*q2^c*r^b]],Nothing],{a,0,nMax},{c,0,nMax},{b,-2 nMax,2 nMax}];
Total[terms]/. {q1^0->1,q2^0->1,r^0->1}];

(*truncation: keep only q1^a q2^c with 0<=a,c<=nMax (leave r free)*)
Clear[TruncateQ1Q2With];
TruncateQ1Q2With[expr_,nMax_Integer?NonNegative,q1sym_Symbol,q2sym_Symbol]:=Sum[Sum[With[{coeff=Coefficient[Coefficient[expr,q1sym,a],q2sym,c]},coeff*q1sym^a*q2sym^c],{c,0,nMax}],{a,0,nMax}]/. {q1sym^0->1,q2sym^0->1};

(*Siegel-Eisenstein expansions of weight w*)
Clear[Ek];
Ek[w_,n_,q1_,r_,q2_]:=SiegelEisensteinGenus2Expansion[w,n,q1,r,q2];

(*\[Chi]_10*)
igusaK10=-43867/(2^12*3^5*5^2*7*53);
Clear[Chi10];
Chi10[nMax_Integer?NonNegative,q1_,r_,q2_]:=Module[{raw},raw=igusaK10*(Ek[4,nMax,q1,r,q2]*Ek[6,nMax,q1,r,q2]-Ek[10,nMax,q1,r,q2]);
Expand@TruncateQ1Q2With[raw,nMax,q1,q2]];

(*\[Chi]_12*)
igusaK12=(131*593)/(2^13*3^7*5^3*7^2*337);
Clear[Chi12];
Chi12[nMax_Integer?NonNegative,q1_,r_,q2_]:=Module[{raw},raw=igusaK12*((3^2*7^2)*Ek[4,nMax,q1,r,q2]^3+(2*5^3)*Ek[6,nMax,q1,r,q2]^2-691*Ek[12,nMax,q1,r,q2]);
Expand@TruncateQ1Q2With[raw,nMax,q1,q2]];


Clear[SeriesWithExp];
SeriesWithExp[seriesFun_,nMax_:order][q1expr_,rexpr_,q2expr_]:=Module[{Q1=Unique["q1$"],R=Unique["r$"],Q2=Unique["q2$"]},seriesFun[nMax,Q1,R,Q2]/. {Q1->q1expr,R->rexpr,Q2->q2expr}];

(* \[Chi]12 *)
Clear[chi12s,chi12sb];
chi12s[q1_,r_,q2_,nMax_:order]:=SeriesWithExp[Chi12,nMax][q1,r,q2];
chi12sb[q1b_,rb_,q2b_,nMax_:order]:=SeriesWithExp[Chi12,nMax][q1b,rb,q2b];

(* \[Chi]10 *)
Clear[chi10s,chi10sb];
chi10s[q1_,r_,q2_,nMax_:order]:=SeriesWithExp[Chi10,nMax][q1,r,q2];
chi10sb[q1b_,rb_,q2b_,nMax_:order]:=SeriesWithExp[Chi10,nMax][q1b,rb,q2b];

(* E4 and E6 via Ek *)
Clear[E4n,E6n, E10n,E12n,e4s,e4sb,e6s,e6sb,e10s,e10sb,e12s,e12sb];

E4n[n_,q1_,r_,q2_]:=Ek[4,n,q1,r,q2];
E6n[n_,q1_,r_,q2_]:=Ek[6,n,q1,r,q2];
E10n[n_,q1_,r_,q2_]:=Ek[10,n,q1,r,q2];
E12n[n_,q1_,r_,q2_]:=Ek[12,n,q1,r,q2];

e4s[q1_,r_,q2_,nMax_:order]:=SeriesWithExp[E4n,nMax][q1,r,q2];
e4sb[q1b_,rb_,q2b_,nMax_:order]:=SeriesWithExp[E4n,nMax][q1b,rb,q2b];

e6s[q1_,r_,q2_,nMax_:order]:=SeriesWithExp[E6n,nMax][q1,r,q2];
e6sb[q1b_,rb_,q2b_,nMax_:order]:=SeriesWithExp[E6n,nMax][q1b,rb,q2b];

e10s[q1_,r_,q2_,nMax_:order]:=SeriesWithExp[E10n,nMax][q1,r,q2];
e10sb[q1b_,rb_,q2b_,nMax_:order]:=SeriesWithExp[E10n,nMax][q1b,rb,q2b];

e12s[q1_,r_,q2_,nMax_:order]:=SeriesWithExp[E12n,nMax][q1,r,q2];
e12sb[q1b_,rb_,q2b_,nMax_:order]:=SeriesWithExp[E12n,nMax][q1b,rb,q2b];

(* formal q-variables used for precomputation *)
ClearAll[$q1,$r,$q2];
SetAttributes[{chi12Base,chi10Base,e4Base,e6Base,chi12Deriv,chi10Deriv,e4Deriv,e6Deriv},HoldAll];
ClearAll[op1,opAll,diffOp];
op1[expr_,sym_,k_Integer?NonNegative]:=Nest[sym*D[#,sym]&,expr,k];
opAll[expr_,a_,b_,c_]:=op1[op1[op1[expr,$q1,a],$r,b],$q2,c];
diffOp[expr_,a_,b_,c_]:=(2 Pi I)^(a+b+c)*opAll[expr,a,b,c];

(* bases *)
chi12Base[ord_Integer]:=chi12Base[ord]=chi12s[$q1,$r,$q2,ord];
chi10Base[ord_Integer]:=chi10Base[ord]=chi10s[$q1,$r,$q2,ord];
e4Base[ord_Integer]:=e4Base[ord]=e4s[$q1,$r,$q2,ord];
e6Base[ord_Integer]:=e6Base[ord]=e6s[$q1,$r,$q2,ord];
e12Base[ord_Integer]:=e12Base[ord]=e12s[$q1,$r,$q2,ord];
e10Base[ord_Integer]:=e10Base[ord]=e10s[$q1,$r,$q2,ord];

(* derivatives *)
chi12Deriv[ord_,a_,b_,c_]:=chi12Deriv[ord,a,b,c]=diffOp[chi12Base[ord],a,b,c];
chi10Deriv[ord_,a_,b_,c_]:=chi10Deriv[ord,a,b,c]=diffOp[chi10Base[ord],a,b,c];
e4Deriv[ord_,a_,b_,c_]:=e4Deriv[ord,a,b,c]=diffOp[e4Base[ord],a,b,c];
e6Deriv[ord_,a_,b_,c_]:=e6Deriv[ord,a,b,c]=diffOp[e6Base[ord],a,b,c];
e10Deriv[ord_,a_,b_,c_]:=e10Deriv[ord,a,b,c]=diffOp[e10Base[ord],a,b,c];
e12Deriv[ord_,a_,b_,c_]:=e12Deriv[ord,a,b,c]=diffOp[e12Base[ord],a,b,c];

(* substitution for \[PlusMinus]sign in the exponent when taking derivatives *)
subExp[T_,Z_,U_,sgn_:1]:={$q1->Exp[sgn*2 Pi I T],$r->Exp[sgn*2 Pi I Z],$q2->Exp[sgn*2 Pi I U]};

(*build rule set at a fixed order*)
ClearAll[PrecomputedRules];
PrecomputedRules[ord_Integer:16] := Module[{phase},
  phase[a_, b_, c_, sgn_] := If[sgn === -1, (-1)^(a + b + c), 1];
  {HoldPattern[\[Chi]12[T_, Z_, U_]] :>(chi12Base[ord] /. subExp[T, Z, U, 1]),
   HoldPattern[\[Chi]10[T_, Z_, U_]] :>(chi10Base[ord] /. subExp[T, Z, U, 1]),
   HoldPattern[E4[T_, Z_, U_]] :>(e4Base[ord] /. subExp[T, Z, U, 1]),
   HoldPattern[E6[T_, Z_, U_]] :>(e6Base[ord] /. subExp[T, Z, U, 1]),
   HoldPattern[E10[T_, Z_, U_]] :>(e10Base[ord] /. subExp[T, Z, U, 1]),
   HoldPattern[E12[T_, Z_, U_]] :>(e12Base[ord] /. subExp[T, Z, U, 1]),
   HoldPattern[Derivative[a_, b_, c_][\[Chi]12][T_, Z_, U_]] :>(chi12Deriv[ord, a, b, c] /. subExp[T, Z, U, 1]),
   HoldPattern[Derivative[a_, b_, c_][\[Chi]10][T_, Z_, U_]] :>(chi10Deriv[ord, a, b, c] /. subExp[T, Z, U, 1]),
   HoldPattern[Derivative[a_, b_, c_][E4][T_, Z_, U_]] :>(e4Deriv[ord, a, b, c] /. subExp[T, Z, U, 1]),
   HoldPattern[Derivative[a_, b_, c_][E6][T_, Z_, U_]] :>(e6Deriv[ord, a, b, c] /. subExp[T, Z, U, 1]),
   HoldPattern[Derivative[a_, b_, c_][E10][T_, Z_, U_]] :>(e10Deriv[ord, a, b, c] /. subExp[T, Z, U, 1]),   
   HoldPattern[Derivative[a_, b_, c_][E12][T_, Z_, U_]] :>(e12Deriv[ord, a, b, c] /. subExp[T, Z, U, 1])}
];


End[];
EndPackage[];
