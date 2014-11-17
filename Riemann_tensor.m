(*functions and packages in Mathematica*)


metric[lineelement_, independentvars_List] := Block[
   {lenindependent, differentials, diffmatrix, metricform, varmetric, 
    gh, sum, equation, rule, varhelp, zeros, zerorule},
   lenindependent = Length[independentvars];
   differentials = Map[Dt, independentvars];
   diffmatrix = Outer[Times, differentials, differentials];
   metricform = Array[gh, {lenindependent, lenindependent}];
   varmetric = Variables[metricform];
   If[Length[metricform] == Length[diffmatrix],
    sum = 0;
    Do[
     Do[
      sum = sum + metricform[[i, j]] diffmatrix[[i, j]],
      {j, 1, lenindependent}], {i, 1, lenindependent}],
    sum = 0];
   If[sum === 0, Return[sum],
    sum = sum - lineelement;
    equation = CoefficientList[sum, differentials] == 0;
    rule = Solve[equation, varmetric];
    metricform = metricform /. rule;
    varmetric = Variables[metricform];
    varhelp = {};
    Do[
     If[Not[FreeQ[varmetric[[i]], gh]], 
      AppendTo[varhelp, varmetric[[i]]]],
     {i, 1, Length[varmetric]}];
    zeros = Table[0, {Length[varhelp]}];
    SubstRule[x_, y_] := x -> y;
    zerorule = Thread[SubstRule[varhelp, zeros]];
    metricform = Flatten[metricform /. zerorule, 1];
    metricform = Expand[(metricform + Transpose[metricform])/2]
    ];
   metricform
   ];
Off[Solve::svars];
Off1[Solve::svars];

Christoffel[m_, a_, b_, g_, ginv_] := Block[{n},
  Expand[
   Sum[ginv[[m, n]] (D[g[[n, a]], IndepVar[[b]]] +
        D[g[[n, b]], IndepVar[[a]]] -
        D[g[[a, b]], IndepVar[[n]]]),
     {n, 1, Length[g]}]/2]]    (*check this function*)

Riemann[a_, b_, c_, d_, g_, ing_] := Block[{},
  Expand[
   D[Christoffel[a, b, d, g, ing], IndepVar[[c]]] -
    D[Christoffel[a, b, c, g, ing], IndepVar[[d]]] +
    Sum[Christoffel[e, b, d, g, ing]
      Christoffel[a, e, c, g, ing],
     {e, 1, Length[g]}] -
    Sum[Christoffel[e, b, c, g, ing]
      Christoffel[a, e, d, g, ing],
     {e, 1, Length[g]}]
   ]
  ]

Ricci[m_, q_, g_, ing_] := Block[{a},
  Expand[
   Sum[Riemann[a, m, a, q, g, ing],
    {a, 1, Length[g]}]]] (*curvature scalar*)

RicciScalar[g_, ing_] := Block[{},
  Expand[Sum[ing[[a, b]] Ricci[a, b, g, ing],
    {a, 1, Length[g]}, {b, 1, Length[g]}]]]

Einstein[m_, n_, g_, ing_] :=
 Ricci[m, n, g, ing] - 
  RicciScalar[g, ing] g[[m, n]]/
    2 (*gives derived nonlinear partial differential equations of second order*)

Bianchi[a_, g_, ing_] := Block[{},
  Expand[
   Sum[D[Sum[ing[[n, m]] Einstein[m, a, g, ing],
       {m, 1, Length[g]}], IndepVar[[n]]],
     {n, 1, Length[g]}]
    + Sum[Sum[Christoffel[n, m, n, g, ing]
       Sum[ing[[m, 1]] Einstein[l, a, g, ing],
        {l, 1, Length[g]}], {m, 1, Length[g]}],
     {n, 1, Length[g]}]
    - Sum[Sum[Christoffel[n, m, a, g, ing]
       Sum[ing[[m, l]] Einstein[l, n, g, ing],
        {l, 1, Length[g]}], {m, l, Length[g]}],
     {n, 1, Length[g]}]
   ]
  ](*Bianchi identities: energy conservation*)
 
MatrixForm[metric[tx Dt[t]^2 + x Dt[x]^2, {x, t}]]
g = metric[Dt[x]^2 + Dt[y]^2 + Dt[z]^2, {x, y, z}]  (*the metric tensor on the Cartesian space*)
invg = Inverse[g]
Christoffel[1, 1, 1, g, invg] (*no nonzero Christoffel symbols of this Cartesian metric*)
Ricci[1, 2, g, invg] (*no nonzero Ricci tensor*)

h = metric[Dt[r]^2 + Dt[z]^2 + r^2 Dt[p]^2, {r, p, z}] (*in cylindrical coordinates*)
invh = Inverse[h]
Table[Christoffel[i, j, k, h, invh], {i, 1, 3}, {j, 1, 3}, {k, 1, 3}]
Table[Riemann[a, b, c, d, h, invh], {a, 1, 3}, {b, 1, 3}, {c, 1, 3}, {d, 1, 3}]
  
l = metric[Dt[r]^2 + r^2 Dt[th]^2 + r^2 Sin[th]^2 Dt[phi]^2, {r, th, phi}] (*in spherical coordinates*)
invl = Inverse[l]
Table[Christoffel[i, j, k, l, invl], {i, 1, 3}, {j, 1, 3}, {k, 1, 3}]
Table[Riemann[a, b, c, d, l, invl], {a, 1, 3}, {b, 1, 3}, {c, 1, 3}, {d, 1, 3}]

u = metric[-r^2 Dt[th]^2 - r^2 Dt[phi]^2 Sin[th]^2 - Exp[lambda[r]] Dt[r]^2 + Exp[nu[r]] Dt[t]^2, {t, r, th, phi}] (*spherically symmetric metric for a charged mass point*)
invu = Inverse[u]



  
