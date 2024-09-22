(*Input m and M*)
m=Input["Input m(19):"];
M=Input["Input M:"];
Print["Your input numbers are m="<>ToString[m]<>" and M="<>ToString[M] ];
Pause[0.4];

(*calculation parameters*)
\[Lambda]=100^2;
prg=50;
DSEerf=0.007;
kn=16;
(*Real DSE*)
(*Because of the symmetry of y, we could not calculate y*)
realnk=150;
realnz=24;
(*Complex DSE*)
complexnk=40;
(*Because of the symmetry of y, we could not calculate y*)
complexnz=24;
(*BSE*)
BSEnk=40;
BSEnz=24;
BSEny=24;
BSEerf=Exp[-10];
(*judge procedure*)
While[True,
calpara=Input["Your calculation parameters:
    1. \[Lambda]="<>ToString[\[Lambda]]<>"
    2. DSEerf="<>ToString[DSEerf]<>"
    3. prg="<>ToString[prg]<>"
    4. realnk="<>ToString[realnk]<>"
    5. realnz="<>ToString[realnz]<>"
    6. complexnk="<>ToString[complexnk]<>"
    7. complexnz="<>ToString[complexnz]<>"
    8. kn="<>ToString[kn]<>"
    9. BSEnk="<>ToString[BSEnk]<>"
    10. BSEnz="<>ToString[BSEnz]<>"
    11. BSEny="<>ToString[BSEny]<>"
    12. BSEerf="<>ToString[N[BSEerf] ]<>"
    Do you want to change these parameters?(Answer: 'y' or 'n')"];
    Which[
        ToString[calpara]=="y",num=Input["Please input the number of the parameter you want to change:"];
            Which[
            num==1,\[Lambda]=Input["Please input the new value of \[Lambda]:"],
            num==2,DSEerf=Input["Please input the new value of DSEerf:"],
            num==3,prg=Input["Please input the new value of prg:"],
            num==4,realnk=Input["Please input the new value of realnk:"],
            num==5,realnz=Input["Please input the new value of realnz:"],
            num==6,complexnk=Input["Please input the new value of complexnk:"],
            num==7,complexnz=Input["Please input the new value of complexnz:"],
            num==8,kn=Input["Please input the new value of kn:"],
            num==9,BSEnk=Input["Please input the new value of BSEnk:"],
            num==10,BSEnz=Input["Please input the new value of BSEnz:"],
            num==11,BSEny=Input["Please input the new value of BSEny:"],
            num==12,BSEerf=Input["Please input the new value of BSEerf:"],
            True,Print["Please check your input!\n"];Pause[1]
            ];
            (*judge*)
            Which[
                complexnk>realnk,Print["The value of complexnk should be less than realnk!"];
                complexnk=realnk;Pause[1]
            ]
            ,
        ToString[calpara]=="n",Break[], 
        True,Print["Please check your input!\n"];Pause[1]
    ];
];
Print["Your calculation parameters:
   1. \[Lambda]="<>ToString[\[Lambda]]<>"
    2. DSEerf="<>ToString[DSEerf]<>"
    3. prg="<>ToString[prg]<>"
    4. realnk="<>ToString[realnk]<>"
    5. realnz="<>ToString[realnz]<>"
    6. complexnk="<>ToString[complexnk]<>"
    7. complexnz="<>ToString[complexnz]<>"
    8. kn="<>ToString[kn]<>"
    9. BSEnk="<>ToString[BSEnk]<>"
    10. BSEnz="<>ToString[BSEnz]<>"
    11. BSEny="<>ToString[BSEny]<>"
    12. BSEerf="<>ToString[N[BSEerf] ]
    <>"\n"
    ];
Pause[0.5];
(*module parameters*)
(*D=1.024 GeV^2*)
d=1.024;
(*\[Omega]=0.5 GeV*)
\[Omega]=0.5;
(*m_t=0.5 GeV*)
mt=0.5;
(*\[Tau]=e^2-1 GeV*)
\[Tau]=E^2-1;
(*Nf=4*)
Nf=4;
(*\[Lambda]QCD=0.234 GeV*)
\[Lambda]QCD=0.234;
(*\[Gamma]m=12/(33-2 Nf)*)
\[Gamma]m=12/(33-2 Nf);
(*renormalised point*)
\[Mu]=19;
(*\[Eta]=0.5*)
\[Eta]=0.5;
(*judge procedure*)
While[True,
modpara=Input["Your module parameters:
    1. d="<>ToString[d]<>"
    2. \[Omega]="<>ToString[\[Omega]]<>"
    3. mt="<>ToString[mt]<>"
    4. Nf="<>ToString[Nf]<>"
    5. \[Lambda]QCD="<>ToString[\[Lambda]QCD]<>"
    6. \[Mu]="<>ToString[\[Mu]]<>"
    7. \[Eta]="<>ToString[\[Eta]]<>"
    Do you want to change these parameters?(Answer: 'y' or 'n')"];
    Which[
        ToString[modpara]=="y",num=Input["Please input the number of the parameter you want to change:"];
            Which[
            num==1,d=Input["Please input the new value of d:"],
            num==2,\[Omega]=Input["Please input the new value of \[Omega]:"],
            num==3,mt=Input["Please input the new value of mt:"],
            num==4,Nf=Input["Please input the new value of Nf:"],
            num==5,\[Lambda]QCD=Input["Please input the new value of \[Lambda]QCD:"],
            num==6,Print["You can not change this parameter!\n"],
            num==7,\[Eta]=Input["Please input the new value of \[Eta]:"],
            True,Print["Please check your input!\n"];Pause[1]
            ],
        ToString[calpara]=="n",Break[], 
        True,Print["Please check your input!\n"];Pause[1]
    ];
];

(*Export*)
path = FileNameJoin[{"parameters", "real-DSE-parameters" <> ".mx"}];
If[$OperatingSystem === "Windows", 
path = StringReplace[path, "\\" -> "\\\\"]
];
DumpSave[path, {m, \[Lambda], DSEerf, prg, realnk, realnz, kn, d, \[Omega], mt, \[Tau], Nf, \[Lambda]QCD, \[Gamma]m, \[Mu]}];
Clear[path];

path = FileNameJoin[{"parameters", "complex-DSE-parameters" <> ".mx"}];
If[$OperatingSystem === "Windows", 
path = StringReplace[path, "\\" -> "\\\\"]
];
DumpSave[path, {m, M, \[Lambda], DSEerf, prg, complexnk, complexnz, kn, \[Eta], d, \[Omega], mt, \[Tau], Nf, \[Lambda]QCD, \[Gamma]m, \[Mu]}];
Clear[path];

path = FileNameJoin[{"parameters", "BSE-parameters" <> ".txt"}];
If[$OperatingSystem === "Windows", 
path = StringReplace[path, "\\" -> "\\\\"]
];
Export[path, 
    {N[d], N[\[Omega]], N[mt], N[\[Tau]], N[\[Lambda]QCD], N[\[Gamma]m], N[M], N[\[Lambda]],
    N[BSEerf], BSEnk, BSEnz, BSEny}];
Clear[path];