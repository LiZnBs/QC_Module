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

CalArBr[a0_,kp0_,zp0_,z20_,z40_]:=Module[
    {
        s1=Table[0,{Length[zp0[[1]]]}],
        s2=Table[0,{Length[zp0[[1]]]}],
        P={0,0,0,I M},
        Ac=Table[0,{Length[kp0[[1]]]}],
        Bc=Table[0,{Length[kp0[[1]]]}]
    },

    (*Relative momentum*)
    p[z_]:={0, 0, Sqrt[\[Lambda] (1-z^2)], Sqrt[\[Lambda] ] z};
    (*Integrate momentum*)
    Q[q_,y_,z_]:={0, Sqrt[q(1-z^2)(1-y^2)], y Sqrt[q(1-z^2)], Sqrt[q] z};

    fa[z_,q_,y_,qz_,A_,B_]:=G[(p[z]+\[Eta] P-Q[q,y,qz]).(p[z]+\[Eta] P-Q[q,y,qz])]*
    A/(q*A^2+B^2)*
    ((p[z]+\[Eta] P).Q[q,y,qz]+2((p[z]+\[Eta] P-Q[q,y,qz]).Q[q,y,qz])((p[z]+\[Eta] P-Q[q,y,qz]).(p[z]+\[Eta] P))/((p[z]+\[Eta] P-Q[q,y,qz]).(p[z]+\[Eta] P-Q[q,y,qz])));

    fb[z_,q_,y_,qz_,A_,B_]:=G[(p[z]+\[Eta] P-Q[q,y,qz]).(p[z]+\[Eta] P-Q[q,y,qz])]*
    B/(q*A^2+B^2);
    Do[
        Which[
            kp0[[1]][[i]]<a0[[1]][[1]],
            Ac[[i]]=a0[[2]][[1]];
            Bc[[i]]=a0[[3]][[1]],
            kp0[[1]][[i]]>a0[[1]][[Length[a0[[1]] ] ]],
            Ac[[i]]=a0[[2]][[Length[a0[[1]] ] ]];
            Bc[[i]]=a0[[3]][[Length[a0[[1]] ] ]],
            kp0[[1]][[i]]>=a0[[1]][[1]]&& kp0[[1]][[i]]<=a0[[1]][[Length[a0[[1]] ] ]],
            Do[
                If[a0[[1]][[j]]<=kp0[[1]][[i]]&&a0[[1]][[j+1]]>=kp0[[1]][[i]],
                Ac[[i]]=(kp0[[1]][[i]]-a0[[1]][[j]])/(a0[[1]][[j+1]]-a0[[1]][[j]])*a0[[2]][[j+1]]+(a0[[1]][[j+1]]-kp0[[1]][[i]])/(a0[[1]][[j+1]]-a0[[1]][[j]])*a0[[2]][[j]];
                Bc[[i]]=(kp0[[1]][[i]]-a0[[1]][[j]])/(a0[[1]][[j+1]]-a0[[1]][[j]])*a0[[3]][[j+1]]+(a0[[1]][[j+1]]-kp0[[1]][[i]])/(a0[[1]][[j+1]]-a0[[1]][[j]])*a0[[3]][[j]];
                Break[];
                ]
            ,{j,1,Length[a0[[1]] ]}
            ]
        ] 
        ,{i,1,Length[kp0[[1]] ]}
    ];
    CloseKernels[];
    LaunchKernels[kn];
    SetSharedVariable[s1,s2];
    ParallelDo[
        Quiet[
        s1[[i]]=z20+4 z20^2/(3 (p[zp0[[1]][[i]] ]+\[Eta] P).(p[zp0[[1]][[i]] ]+\[Eta] P))*(2\[Pi])^(-3)*Sum[
            kp0[[2]][[l]]*zp0[[2]][[n]]*zp0[[2]][[m]]/2*kp0[[1]][[l]]*Sqrt[1-zp0[[1]][[n]]^2]*Log[\[Lambda] ]*kp0[[1]][[l]]*
            fa[zp0[[1]][[i]],kp0[[1]][[l]],zp0[[1]][[m]],zp0[[1]][[n]],Ac[[l]],Bc[[l]] ]
            ,{l,1,Length[kp0[[1]] ]},{n,1,Length[zp0[[1]] ]},{m,1,Length[zp0[[1]] ]}
        ];
        s2[[i]]=z40 m+4 z20^2*(2\[Pi])^(-3)*Sum[
            kp0[[2]][[l]]*zp0[[2]][[n]]*zp0[[2]][[m]]/2*kp0[[1]][[l]]*Sqrt[1-zp0[[1]][[n]]^2]*Log[\[Lambda] ]*kp0[[1]][[l]]*
            fb[zp0[[1]][[i]],kp0[[1]][[l]],zp0[[1]][[m]],zp0[[1]][[n]],Ac[[l]],Bc[[l]] ]
            ,{l,1,Length[kp0[[1]] ]},{n,1,Length[zp0[[1]] ]},{m,1,Length[zp0[[1]] ]}
        ]
        , {General::munfl}]

        ,{i,1,Length[zp0[[1]] ]} 
    ];
    {s1,s2}
]
