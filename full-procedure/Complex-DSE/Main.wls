(*Files*)
files = {"complex-DSE-parameters.mx"};
paths = FileNameJoin[{ParentDirectory[],"parameters",#}] & /@ files;
If[$OperatingSystem === "Windows", 
    paths = StringReplace[paths, "\\" -> "\\\\"]
];
Get /@ paths;
Clear[files, paths];
Print[M];
files = {"Gauleg-method_Get-point.wls"};
paths = FileNameJoin[{"Package", "Functions",#}] & /@ files;
If[$OperatingSystem === "Windows", 
    paths = StringReplace[paths, "\\" -> "\\\\"]
];
Get /@ paths;
Clear[files, paths];

files = {"Cal-ApBp_Input-AB.wls","Cal-AB_Input-ApBp.wls","Cal-ArBr.wls"};
paths = FileNameJoin[{#}] & /@ files;
If[$OperatingSystem === "Windows", 
    paths = StringReplace[paths, "\\" -> "\\\\"]
];
Get /@ paths;
Clear[files, paths];
(*Get the renormalization constant of Real DSE*)
files = {"Real-dse_Z2-Z4.wls"};
paths = FileNameJoin[{"Parameters",#}] & /@ files;
If[$OperatingSystem === "Windows", 
    paths = StringReplace[paths, "\\" -> "\\\\"]
];
Z= Get /@ paths;
Clear[files, paths];

(*Get A,B of Real DSE*)
files = {"Real-dse_A-B_p2_All-results.wls"};
paths = FileNameJoin[{"Parameters",#}] & /@ files;
If[$OperatingSystem === "Windows", 
    paths = StringReplace[paths, "\\" -> "\\\\"]
];
t= Get /@ paths;
Clear[files, paths];
t=Flatten[t,1];(*Flatten the list {A,B}*)
realresultA=Table[{t[[1]][[i]],t[[2]][[i]]},{i,1,Length[t[[1]] ]}];
realresultB=Table[{t[[1]][[i]],t[[3]][[i]]},{i,1,Length[t[[1]] ]}];
(*Initial*)
kp=gauleg[-1,1,complexnk];zp=gauleg[-1,1,complexnz];
(*E index map*)
kp=ReplacePart[kp,1->\[Lambda]^kp[[1]] ];
(*Initial value*)
(*Loop points*)
A1u=Table[1,{complexnk}];B1u=Table[0.0034,{complexnk}];
A1d=Table[1,{complexnk}];B1d=Table[0.0034,{complexnk}];
(*Inside point*)
(*Ap[[zp]][[kp]]*)
Ap=Table[0,{complexnz},{complexnk}];Bp=Table[0,{complexnz},{complexnk}];
(*Loop points 2*)
A2u=Table[0,{complexnk}];B2u=Table[0,{complexnk}];
A2d=Table[0,{complexnk}];B2d=Table[0,{complexnk}];
total=100;
(*Exact value*)
Z2=Z[[1]][[1]];Z4=Z[[1]][[2]];
Ar=Table[0,{complexnz}];Br=Table[0,{complexnz}];
t=CalArBr[t,kp,zp,Z2,Z4];
Ar=t[[1]];Br=t[[2]];
Clear[t];
Print["Start calculating Complex-DSE","\n"];
(*Calculate Ar,Br over, absolutely right!*)
While[
    total>DSEerf
    ,
    t=CalApBp[ A1u, B1u, A1d, B1d, Ar, Br, kp, zp];
    Ap=t[[1]];Bp=t[[2]];
    t=CalAB[ Ap, Bp, kp, zp, Z2, Z4];
    A2u=t[[1]];B2u=t[[2]];A2d=t[[3]];B2d=t[[4]];
    total=Total[MapThread[Abs[Subtract[#1, #2] ] &, {A2u, A1u}] ] + Total[MapThread[Abs[Subtract[#1, #2] ] &, {B2u, B1u}] ]+Total[MapThread[Abs[Subtract[#1, #2] ] &, {A2d, A1d}] ] + Total[MapThread[Abs[Subtract[#1, #2] ] &, {B2d, B1d}] ];
    A1u=A2u;B1u=B2u;A1d=A2d;B1d=B2d;
    Print[total];
]; 
(*Result on the real axis*)
resultra=Table[{kp[[1]][[i]]-(\[Eta] M)^2,Re[Ap[[12]][[i]] ]},{i,1,complexnk}];
resultrb=Table[{kp[[1]][[i]]-(\[Eta] M)^2,Re[Bp[[12]][[i]] ]},{i,1,complexnk}];
Ar=Table[{zp[[1]][[i]],kp[[1]][[j]],Ap[[i]][[j]]},{i,1,complexnz},{j,1,complexnk}];
Br=Table[{zp[[1]][[i]],kp[[1]][[j]],Bp[[i]][[j]]},{i,1,complexnz},{j,1,complexnk}];
(*Select results*)
n=1;i=1;
While[
        resultra[[i]][[1]]<prg,n=n+1;i=i+1;
    ];
resultra=Drop[resultra, n - Length[resultra] - 1];
resultrb=Drop[resultrb, n - Length[resultrb] - 1];

n=1;i=1;
While[
        realresultA[[i]][[1]]<prg,n=n+1;i=i+1;
    ];
realresultA=Drop[realresultA, n - Length[realresultA] - 1];
realresultB=Drop[realresultB, n - Length[realresultB] - 1];

    (*Plot*) 
picA = ListLinePlot[{resultra,realresultA},
    Frame -> True,
    PlotStyle -> {Blue, Red},
    FrameLabel -> {Style["p^2 (GeV)", Bold, 12.5], 
    Style["(GeV)", Bold, 12.5]},
    Epilog -> {Blue,Text["rA", {1, 0.5}], Red, Text["realA", {1, 0.4}]},
    PlotRange -> Full
    ];

picB = ListLinePlot[{resultrb,realresultB},
    Frame -> True,
    PlotStyle -> {Blue, Red},
    FrameLabel -> {Style["p^2 (GeV)", Bold, 12.5], 
    Style["(GeV)", Bold, 12.5]},
    Epilog -> {Blue,Text["rB", {1, 0.5}], Red, Text["realB", {1, 0.4}]},
    PlotRange -> Full
    ];

    (*Export*)
If[Not[DirectoryQ["Result"] ], CreateDirectory["Result"] ];  
path = FileNameJoin[{"Result", "Compare_real-dse_complex-dse-on-real-axis_A" <>"m-"<>ToString[m]<>"_M-"<>ToString[M]<>".png"}];
If[$OperatingSystem === "Windows", 
path = StringReplace[path, "\\" -> "\\\\"]
];
Export[path, picA];
Clear[path];
    
path = FileNameJoin[{"Result", "Compare_real-dse_complex-dse-on-real-axis_B"  <>"m-"<>ToString[m]<>"_M-"<>ToString[M]<> ".png"}];
If[$OperatingSystem === "Windows", 
path = StringReplace[path, "\\" -> "\\\\"]
];
Export[path, picB];
Clear[path];

(*Export into txt*)
path = FileNameJoin[{"Result", "Complex-dse_A.txt"}];
content=Map["("<>ToString[Re[#], CForm] <> ", " <> ToString[Im[#], CForm]<>")" &, Ap, {2}];
content=Flatten[content,1];
Export[path, content, "Text"];

path = FileNameJoin[{"Result", "Complex-dse_B.txt"}];
content=Map["("<>ToString[Re[#], CForm] <> ", " <> ToString[Im[#], CForm]<>")" &, Bp, {2}];
content=Flatten[content,1];
Export[path, content, "Text"];
Clear[path,content];