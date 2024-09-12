(*Calculate the boundary line number of A and B*)
files = {"Function_Cal-G_Input-k2.wls"};
paths = FileNameJoin[{"Package", "Functions",#}] & /@ files;
If[$OperatingSystem === "Windows", 
    paths = StringReplace[paths, "\\" -> "\\\\"]
];
Get /@ paths;
Clear[files, paths];

files = {"Parameters-in-QCM.wls","kn.wls"};
paths = FileNameJoin[{"Package", "Parameters",#}] & /@ files;
If[$OperatingSystem === "Windows", 
    paths = StringReplace[paths, "\\" -> "\\\\"]
];
Get /@ paths;
Clear[files, paths];

files = {"Calculation-parameters.wls"};
paths = FileNameJoin[{"Parameters",#}] & /@ files;
If[$OperatingSystem === "Windows", 
    paths = StringReplace[paths, "\\" -> "\\\\"]
];
Get /@ paths;
Clear[files, paths];

CalAB[a0_,b0_, kp0_,zp0_,z20_,z40_]:=Module[
    {
        s1=Table[0,{Length[kp0[[1]]]}],
        s2=Table[0,{Length[kp0[[1]]]}],
        s3=Table[0,{Length[kp0[[1]]]}],
        s4=Table[0,{Length[kp0[[1]]]}],
        P={0,0,0,I M}
    },

    (*Relative momentum*)
    p[p0_,z_]:={0, 0, Sqrt[p0 (1-z^2)], Sqrt[p0] z};
    (*Integrate momentum*)
    Q[q_,z_]:={0, 0,  Sqrt[q (1-z^2)], Sqrt[q] z};

    fa[p0_,q_,z_,A_,B_,t_]:=G[(p[p0,t]-Q[q,z]).(p[p0,t]-Q[q,z])]*
    A/((Q[q,z]+\[Eta] P).(Q[q,z]+\[Eta] P)*A^2+B^2)*
    ((p[p0,t]+\[Eta] P).(Q[q,z]+\[Eta] P)+2((p[p0,t]-Q[q,z]).(Q[q,z]+\[Eta] P))((p[p0,t]-Q[q,z]).(p[p0,t]+\[Eta] P))/((p[p0,t]-Q[q,z]).(p[p0,t]-Q[q,z])));

    fb[p0_,q_,z_,A_,B_,t_]:=G[(p[p0,t]-Q[q,z]).(p[p0,t]-Q[q,z])]*
    B/((Q[q,z]+\[Eta] P).(Q[q,z]+\[Eta] P)*A^2+B^2);

    CloseKernels[];
    LaunchKernels[kn];
    SetSharedVariable[s1,s2,s3,s4];
    ParallelDo[
        Quiet[
        s1[[i]]=z20+4 z20^2/(3 (p[kp0[[1]][[i]], 1]+\[Eta] P).(p[kp0[[1]][[i]], 1]+\[Eta] P))*(2\[Pi])^(-3)*Sum[
            kp0[[2]][[l]]*zp0[[2]][[n]]*kp0[[1]][[l]]*Sqrt[1-zp0[[1]][[n]]^2]*Log[\[Lambda] ]*kp0[[1]][[l]]
            fa[kp0[[1]][[i]],kp0[[1]][[l]],zp0[[1]][[n]],a0[[n]][[l]],b0[[n]][[l]], 1]
            ,{l,1,Length[kp0[[1]] ]},{n,1,Length[zp0[[1]] ]}
        ];
        s2[[i]]=z40 m+4 z20^2*(2\[Pi])^(-3)*Sum[
            kp0[[2]][[l]]*zp0[[2]][[n]]*kp0[[1]][[l]]*Sqrt[1-zp0[[1]][[n]]^2]*Log[\[Lambda] ]*kp0[[1]][[l]]*
            fb[kp0[[1]][[i]],kp0[[1]][[l]],zp0[[1]][[n]],a0[[n]][[l]],b0[[n]][[l]], 1]
            ,{l,1,Length[kp0[[1]] ]},{n,1,Length[zp0[[1]] ]}
        ];
        s3[[i]]=z20+4 z20^2/(3 (p[kp0[[1]][[i]], -1]+\[Eta] P).(p[kp0[[1]][[i]], -1]+\[Eta] P))*(2\[Pi])^(-3)*Sum[
            kp0[[2]][[l]]*zp0[[2]][[n]]*kp0[[1]][[l]]*Sqrt[1-zp0[[1]][[n]]^2]*Log[\[Lambda] ]*kp0[[1]][[l]]
            fa[kp0[[1]][[i]],kp0[[1]][[l]],zp0[[1]][[n]],a0[[n]][[l]],b0[[n]][[l]], -1]
            ,{l,1,Length[kp0[[1]] ]},{n,1,Length[zp0[[1]] ]}
        ];
        s4[[i]]=z40 m+4 z20^2*(2\[Pi])^(-3)*Sum[
            kp0[[2]][[l]]*zp0[[2]][[n]]*kp0[[1]][[l]]*Sqrt[1-zp0[[1]][[n]]^2]*Log[\[Lambda] ]*kp0[[1]][[l]]*
            fb[kp0[[1]][[i]],kp0[[1]][[l]],zp0[[1]][[n]],a0[[n]][[l]],b0[[n]][[l]], -1]
            ,{l,1,Length[kp0[[1]] ]},{n,1,Length[zp0[[1]] ]}
        ]
        , {General::munfl}]

        ,{i,1,Length[kp0[[1]] ]} 
    ];
    {s1,s2,s3,s4} 
]
