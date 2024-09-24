files = {"complex-DSE-parameters.mx"};
paths = FileNameJoin[{ParentDirectory[],"parameters",#}] & /@ files;
If[$OperatingSystem === "Windows", 
    paths = StringReplace[paths, "\\" -> "\\\\"]
];
Get /@ paths;
Clear[files, paths];
M=M*1.0001;
path = FileNameJoin[{ParentDirectory[],"parameters", "complex-DSE-parameters" <> ".mx"}];
If[$OperatingSystem === "Windows", 
path = StringReplace[path, "\\" -> "\\\\"]
];
DumpSave[path, {m, M, \[Lambda], DSEerf, prg, complexnk, complexnz, kn, \[Eta], d, \[Omega], mt, \[Tau], Nf, \[Lambda]QCD, \[Gamma]m, \[Mu]}];
Clear[path];