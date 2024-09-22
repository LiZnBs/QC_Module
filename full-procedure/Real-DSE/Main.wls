(*Input parameters*)
files = {"real-DSE-parameters.mx"};
paths = FileNameJoin[{ParentDirectory[],"parameters",#}] & /@ files;
If[$OperatingSystem === "Windows", 
    paths = StringReplace[paths, "\\" -> "\\\\"]
];
Get /@ paths;
Clear[files, paths];

(*Gauss get point*)
files = {"Gauleg-method_Get-point.wls"};
paths = FileNameJoin[{"Package", "Functions",#}] & /@ files;
If[$OperatingSystem === "Windows", 
    paths = StringReplace[paths, "\\" -> "\\\\"]
];
Get /@ paths;
Clear[files, paths];

(*Calculate procedure*)
files = {"Cal-AsBs_Input-AB.wls","Cal-AB_Input-AsBs.wls"};
paths = FileNameJoin[{#}] & /@ files;
If[$OperatingSystem === "Windows", 
    paths = StringReplace[paths, "\\" -> "\\\\"]
];
Get /@ paths;
Clear[files, paths];

kp=gauleg[-1,1,realnk];zp=gauleg[-1,1,realnz];
(*E index map*)
kp=ReplacePart[kp,1->\[Lambda]^kp[[1]] ];

(*Initial value*)
A1=Table[1,{realnk}];B1=Table[0.0034,{realnk}];As=Table[0,{realnk}];Bs=Table[0,{realnk}];
A2=Table[0,{realnk}];B2=Table[0,{realnk}];(*A2s=Table[1,{100}];B2s=Table[1,{100}];*)
total=100;
Z2=0;Z4=0;
Print["Start calculating Real-DSE","\n"];
While[
    total>DSEerf
    ,
    t=CalAsBs[ A1, B1, kp, zp];
    As=t[[1]];Bs=t[[2]];
    t=CalAB[As,Bs,kp];
    A2=t[[1]];B2=t[[2]];
    total=Total[MapThread[Abs[Subtract[#1, #2] ] &, {A2, A1}] ] + Total[MapThread[Abs[Subtract[#1, #2] ] &, {B2, B1}] ];
    A1=A2;B1=B2;
    Print[total];
];
Z2=2-t[[3]];Z4=2-t[[4]]/m;
    (*Select result*) 
A=Table[0,{i,1,Length[A1]}];
B=Table[0,{i,1,Length[B1]}];
M=Table[0,{i,1,Length[A1]}];
n=1;i=1;
While[
        i<=Length[A] && kp[[1]][[i]]<prg
        ,
        A[[n]]={kp[[1]][[n]],A1[[n]]};
        B[[n]]={kp[[1]][[n]],B1[[n]]};
        M[[n]]={kp[[1]][[n]],B1[[n]]/A1[[n]]};
        n=n+1;i=i+1;
    ];
A=Drop[A, n - Length[A] - 1];
B=Drop[B, n - Length[B] - 1];  
M=Drop[M, n - Length[M] - 1];
    (*Plot*) 
pic = ListLinePlot[{A,B,M},
    Frame -> True,
    PlotStyle -> {Blue, Red,Black},
    FrameLabel -> {Style["p^2 (GeV)", Bold, 12.5], 
    Style["(GeV)", Bold, 12.5]},
    Epilog -> {Blue,Text["A", {1, 0.5}], Red, Text["B", {1, 0.4}],Black,Text["M", {1, 0.3}]},
    PlotRange -> Full
    ];

    (*Export*)
    If[Not[DirectoryQ["Result"] ], CreateDirectory["Result"] ];  
    path = FileNameJoin[{"Result", "Real-dse_" <>StringReplace[DateString[], {" " -> "-", ":" -> "."}]<> ".png"}];
    If[$OperatingSystem === "Windows", 
    path = StringReplace[path, "\\" -> "\\\\"]
    ];
    Export[path, pic];
    Clear[path];

    path = FileNameJoin[{"Result", "Real-dse_A-B_p2_All-results" <> ".wls"}];
    If[$OperatingSystem === "Windows", 
    path = StringReplace[path, "\\" -> "\\\\"]
    ];
    Export[path, {kp[[1]],A1,B1}];
    Clear[path];

    path = FileNameJoin[{"Result", "Real-dse_Z2-Z4" <> ".wls"}];
    If[$OperatingSystem === "Windows", 
    path = StringReplace[path, "\\" -> "\\\\"]
    ];
    Export[path, {Z2,Z4}];
    Clear[path];
(*Export into txt*)
path = FileNameJoin[{"Result", "Z2.txt"}];
Export[path, Z2, "Text"];