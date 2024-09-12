files = {"Renormalization-point.wls"};
paths = FileNameJoin[{"Package", "Parameters",#}] & /@ files;
If[$OperatingSystem === "Windows", 
    paths = StringReplace[paths, "\\" -> "\\\\"]
];
Get /@ paths;
Clear[files, paths];

files = {"Parameters-in-QCM.wls"};
paths = FileNameJoin[{"Package", "Parameters",#}] & /@ files;
If[$OperatingSystem === "Windows", 
    paths = StringReplace[paths, "\\" -> "\\\\"]
];
Get /@ paths;
Clear[files, paths];

CalAB[a0_, b0_, kp0_]:=Module[
    {
    (*A,B*)
    s1,s2,
    (*Z2,Z4 point*)
    r1,r2,
    (*Find \[Mu] place*)
    t=0,i=1
    },
    While[i<Length[kp0[[1]]] && kp0[[1]][[i]]<=\[Mu]^2,
            t=i;i=i+1
        ];
        r1=(\[Mu]^2-kp0[[1]][[t]])/(kp0[[1]][[t+1]]-kp0[[1]][[t]])*a0[[t+1]]+(kp0[[1]][[t+1]]-\[Mu]^2)/(kp0[[1]][[t+1]]-kp0[[1]][[t]])*a0[[t]];
        r2=(\[Mu]^2-kp0[[1]][[t]])/(kp0[[1]][[t+1]]-kp0[[1]][[t]])*b0[[t+1]]+(kp0[[1]][[t+1]]-\[Mu]^2)/(kp0[[1]][[t+1]]-kp0[[1]][[t]])*b0[[t]];
    s1=Table[1+a0[[j]]-r1,{j,1,Length[a0]}];
    s2=Table[m+b0[[j]]-r2,{j,1,Length[b0]}];
    {s1,s2,r1,r2}
]