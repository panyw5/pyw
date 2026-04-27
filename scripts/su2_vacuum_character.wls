#!/usr/bin/env wolframscript
(* 
   su(2) Admissible Module Vacuum Character Computation
   
   For su(2) at level k = -4n/(2n+1), the vacuum character is:
   Ch_0 = theta_1(z | (2n+1)tau) / theta_1(z | tau)
   
   We compute the q-expansion (z-expansion) for n = 1, 2, 3, 4
*)

(* Set computation order *)
order = 20;

(* Jacobi theta_1 function: theta_1(z, q) = -i sum_{n in Z} (-1)^n q^{(n+1/2)^2/2} e^{2 pi i (n+1/2) z} *)
(* Using Mathematica's built-in EllipticTheta[1, pi*z, q^(1/2)] *)

(* Define theta_1 as a q-series expansion *)
Clear[theta1Series];
theta1Series[z_, q_, maxOrder_] := Module[{result, n, nHalf},
  result = Sum[
    (-1)^n * q^((n + 1/2)^2/2) * (Exp[2 Pi I (n + 1/2) z] - Exp[-2 Pi I (n + 1/2) z])/(2 I),
    {n, -maxOrder, maxOrder}
  ];
  (* Simplify using sin *)
  result = -I * Sum[
    (-1)^n * q^((n + 1/2)^2/2) * Exp[2 Pi I (n + 1/2) z],
    {n, -maxOrder, maxOrder}
  ];
  result
];

(* Alternative: use product formula for theta_1 *)
(* theta_1(z, q) = -i q^{1/8} (y^{1/2} - y^{-1/2}) prod_{n>=1} (1-q^n)(1-y q^n)(1-y^{-1} q^n) *)
Clear[theta1Product];
theta1Product[z_, q_, maxOrder_] := Module[{y, qEighth, prod},
  y = Exp[2 Pi I z];
  qEighth = q^(1/8);
  prod = Product[(1 - q^n) * (1 - y * q^n) * (1 - y^(-1) * q^n), {n, 1, maxOrder}];
  -I * qEighth * (y^(1/2) - y^(-1/2)) * prod
];

(* Vacuum character ratio: theta_1(z | u*tau) / theta_1(z | tau) *)
(* In terms of q = e^{2 pi i tau}: theta_1(z, q^u) / theta_1(z, q) *)
Clear[vacuumCharacterRatio];
vacuumCharacterRatio[u_, z_, q_, maxOrder_] := Module[{num, den},
  num = theta1Product[z, q^u, maxOrder];
  den = theta1Product[z, q, maxOrder];
  num / den
];

(* Expand as Laurent series in y = e^{2 pi i z} and q *)
Clear[expandCharacter];
expandCharacter[u_, maxQOrder_] := Module[{ch, expanded, y, simplified},
  (* Use symbolic z, expand in q *)
  ch = vacuumCharacterRatio[u, z, q, maxQOrder + 5];
  
  (* Substitute y = e^{2 pi i z} *)
  ch = ch /. Exp[a_ z] :> y^(a/(2 Pi I)) /. Exp[a_ Pi I z] :> y^(a/(2 Pi I));
  
  (* Series expand in q *)
  expanded = Series[ch, {q, 0, maxQOrder}];
  
  (* Simplify coefficients *)
  expanded = Normal[expanded];
  expanded = Collect[expanded, q, Simplify];
  expanded
];

(* Main computation *)
Print["=== su(2) Admissible Module Vacuum Characters ==="];
Print["Formula: Ch_0 = theta_1(z | (2n+1)tau) / theta_1(z | tau)"];
Print["Level: k = -4n/(2n+1)"];
Print[""];

(* Compute for n = 1, 2, 3, 4 *)
results = {};

Do[
  u = 2 n + 1;
  k = -4 n / (2 n + 1);
  
  Print["========================================"];
  Print["n = ", n, ", u = ", u, ", k = ", k // N];
  Print["========================================"];
  
  (* Compute character ratio *)
  ch = vacuumCharacterRatio[u, z, q, order];
  
  (* Convert to y = e^{2 pi i z} *)
  y = Exp[2 Pi I z];
  
  (* Expand in q, keeping y symbolic *)
  (* First, let's compute numerically at a specific z to verify *)
  zVal = 1/10;
  yVal = Exp[2 Pi I zVal];
  
  chNumeric = ch /. z -> zVal;
  chExpanded = Series[chNumeric, {q, 0, 8}] // Normal;
  
  Print["Numeric expansion at z = 1/10:"];
  Print["Ch_0 = ", chExpanded // Simplify];
  Print[""];
  
  (* Now compute symbolic expansion *)
  (* The character has the form sum_m c_m(y) q^m where c_m(y) are Laurent polynomials in y *)
  
  (* Extract leading q-power *)
  (* For theta_1(z, q^u) / theta_1(z, q), the leading power is (u-1)/8 *)
  leadingPower = (u - 1) / 8;
  Print["Leading q-power: q^(", leadingPower, ")"];
  
  (* Compute first few terms symbolically *)
  Print[""];
  Print["First few terms (symbolic in y = e^{2 pi i z}):"];
  
  (* Use a cleaner approach: expand theta functions directly *)
  num = -I * q^(u/8) * (y^(1/2) - y^(-1/2)) * Product[(1 - q^(u*m)) * (1 - y * q^(u*m)) * (1 - y^(-1) * q^(u*m)), {m, 1, 10}];
  den = -I * q^(1/8) * (y^(1/2) - y^(-1/2)) * Product[(1 - q^m) * (1 - y * q^m) * (1 - y^(-1) * q^m), {m, 1, 10}];
  
  ratio = num / den // Simplify;
  ratio = ratio /. y^(a_) :> y^a;
  
  (* The (y^{1/2} - y^{-1/2}) terms cancel *)
  (* Leading factor: q^{(u-1)/8} *)
  
  ratioExpanded = Series[ratio, {q, 0, 6}] // Normal;
  ratioExpanded = Collect[ratioExpanded, q, Simplify];
  
  Print["Ch_0 = ", ratioExpanded];
  Print[""];
  
  AppendTo[results, {n, u, k, ratioExpanded}];
  
  , {n, 1, 4}
];

Print[""];
Print["=== Summary ==="];
Print[""];

Do[
  {nVal, uVal, kVal, chVal} = results[[i]];
  Print["n = ", nVal, ": k = ", kVal, ", u = ", uVal];
  Print["  Ch_0 = ", chVal];
  Print[""];
  , {i, Length[results]}
];

(* Verify n=1 case matches the user's example *)
Print[""];
Print["=== Verification for n=1 (k = -4/3, u = 3) ==="];
Print["User's formula: Ch_0 = q^(1/4) + (1 + 1/z + z) q^(5/4) + ..."];
Print[""];

(* For n=1, u=3, leading power is (3-1)/8 = 1/4 *)
(* Let's compute more carefully *)

u = 3;
(* theta_1(z, q^3) / theta_1(z, q) *)
(* = q^{(3-1)/8} * [product terms] *)
(* = q^{1/4} * ... *)

(* Expand to higher order *)
num3 = -I * q^(3/8) * (y^(1/2) - y^(-1/2)) * Product[(1 - q^(3*m)) * (1 - y * q^(3*m)) * (1 - y^(-1) * q^(3*m)), {m, 1, 15}];
den3 = -I * q^(1/8) * (y^(1/2) - y^(-1/2)) * Product[(1 - q^m) * (1 - y * q^m) * (1 - y^(-1) * q^m), {m, 1, 15}];

ratio3 = num3 / den3 // Simplify;
ratio3Expanded = Series[ratio3, {q, 0, 15}] // Normal;
ratio3Expanded = Collect[ratio3Expanded, q, Simplify];

Print["Detailed expansion for n=1:"];
Print["Ch_0 = ", ratio3Expanded];
