files = {"Function_Cal-G_Input-k2.wls"};
paths = FileNameJoin[{"Package", "Functions",#}] & /@ files;
If[$OperatingSystem === "Windows", 
    paths = StringReplace[paths, "\\" -> "\\\\"]
];
Get /@ paths;
Clear[files, paths];

(*Input list A , list B , Gauss point x*)
CalAsBs[ a0_, b0_, kp0_,zp0_]:=Module[
    {s1=Table[0,{Length[a0]}],s2=Table[0,{Length[b0]}]
    },
    (*Functions*)
    (*Input p=p2*)
    P[p_]:={0,0,0,Sqrt[p]};
    (*Because of the symmetry of y, we do not need to calculate y*)
    (*Q[q_,y_,z_]:={0, Sqrt[q(1-z^2)(1-y^2)], y Sqrt[q(1-z^2)], Sqrt[q] z};*)
    Q[q_,z_]:={0, Sqrt[q(1-z^2)], 0, Sqrt[q] z};

    fa[p_,q_,z_,A_,B_]:=G[(P[p]-Q[q,z]).(P[p]-Q[q,z])](*/((P[p]-Q[q,z]).(P[p]-Q[q,z]))*)*
    A/(q*A^2+B^2)*
    (P[p].Q[q,z]+2((P[p]-Q[q,z]).Q[q,z])((P[p]-Q[q,z]).P[p])/((P[p]-Q[q,z]).(P[p]-Q[q,z])));

    fb[p_,q_,z_,A_,B_]:=G[(P[p]-Q[q,z]).(P[p]-Q[q,z])](*/((P[p]-Q[q,z]).(P[p]-Q[q,z]))*)*
    B/(q*A^2+B^2);
    
    CloseKernels[];
    LaunchKernels[kn];
    SetSharedVariable[s1,s2];
    ParallelDo[
        (*Because of the symmetry of y, we do not need to calculate y*)
        
        Quiet[
        s1[[i]]=1+4/(3 kp0[[1]][[i]])*(2\[Pi])^(-3)*Sum[
            kp0[[2]][[l]]*zp0[[2]][[n]]*kp0[[1]][[l]](*/2*)*Sqrt[1-zp0[[1]][[n]]^2](*E index map*)*Log[\[Lambda] ]*kp0[[1]][[l]]*
            fa[kp0[[1]][[i]],kp0[[1]][[l]],zp0[[1]][[n]],a0[[l]],b0[[l]]]
            ,{l,1,Length[kp0[[1]]]},{n,1,Length[zp0[[1]]]}
        ];
        s2[[i]]=m+4*(2\[Pi])^(-3)*Sum[
            kp0[[2]][[l]]*zp0[[2]][[n]]*kp0[[1]][[l]](*/2*)*Sqrt[1-zp0[[1]][[n]]^2](*E index map*)*Log[\[Lambda] ]*kp0[[1]][[l]]*
            fb[kp0[[1]][[i]],kp0[[1]][[l]],zp0[[1]][[n]],a0[[l]],b0[[l]]]
            ,{l,1,Length[kp0[[1]]]},{n,1,Length[zp0[[1]]]}
        ]
        , {General::munfl}]
        ,{i,1,Length[a0]}
        
    ];
    {s1,s2}
]

