(*Gauss-Legendre get point*)
gauleg[x1_, x2_, n_] := Module[
    {EPS = 3*10^-16, m, xm, xl, z, z1, p1, p2, p3, pp, x, w},
    m = Quotient[n + 1, 2];
    xm = 0.5 * (x2 + x1);
    xl = 0.5 * (x2 - x1);
    x = ConstantArray[0., n];
    w = ConstantArray[0., n];
    Do[
    z = Cos[Pi * (i - 0.25) / (n + 0.5)];
    While[True,
        p1 = 1.;
        p2 = 0.;
        Do[
            p3 = p2;
            p2 = p1;
            p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j,
            {j, 1, n}
        ];
    pp = n * (z * p1 - p2) / (z^2 - 1);
    z1 = z;
    z = z1 - p1 / pp;
    If[Abs[z - z1] <= EPS, Break[] ]
    ];
    x[[i]] = xm - xl * z;
    x[[n + 1 - i]] = xm + xl * z;
    w[[i]] = 2 * xl / ((1 - z^2) * pp^2);
    w[[n + 1 - i]] = w[[i]],
    {i, 1, m}
    ];
    {x, w}
];
