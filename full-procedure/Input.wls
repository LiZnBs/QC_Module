(*m(19)=0.0034*)
Print["Input m(19):"];
m=Input[];
Print["M:"];
M=Input[];
Print["Your input numbers are m=",m," and M=",M];
path = FileNameJoin[{"Real-DSE", "Parameters", "m-M" <> ".wls"}];
If[$OperatingSystem === "Windows", 
path = StringReplace[path, "\\" -> "\\\\"]
];
Export[path, {m,M}];
Clear[path];

path = FileNameJoin[{"Complex-DSE", "Parameters", "m-M" <> ".wls"}];
If[$OperatingSystem === "Windows", 
path = StringReplace[path, "\\" -> "\\\\"]
];
Export[path, {m,M}];
Clear[path];

path = FileNameJoin[{"BSE", "DSE_Results", "M" <> ".txt"}];
If[$OperatingSystem === "Windows", 
path = StringReplace[path, "\\" -> "\\\\"]
];
Export[path, M, "Text"];
Clear[path];
