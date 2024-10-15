m={0.002,0.0034,0.004,0.005,0.006,0.007,0.008,0.01}
M={0.102,0.134,0.1454,0.163,0.1785,0.1935,0.2072,0.2325}
f={0.0925411,0.093628,0.094167,0.094902,0.0956966,0.0964696,0.0972156,0.0982149}
r={0.241155,0.246978,0.249134,0.251990,0.254830,0.257463,0.259971,0.264597}
p=Table[{m[[i]],2*r[[i]]/f[[i]]*m[[i]]},{i,1,Length[m]}]
pf=Table[{m[[i]],f[[i]]},{i,1,Length[m]}]
pr=Table[{m[[i]],r[[i]]},{i,1,Length[m]}]
picp=ListPlot[p,
Frame->True,
PlotStyle->{Blue},
PlotRange->{{0,0.01},{0,0.06}},
FrameLabel->{Style["m (GeV)",Bold,12.5],Style[Row[{Superscript["M", 2], " (",Superscript["GeV", 2], ")"}], Bold, 12.5]}]
(* Fix[pp_] := Module[
        {
        model, p1, p2, f(*parameters*)
        },

    model = p1*x + p2;
    f = FindFit[pp, model, {p1, p2}, x];
    model /. f
    ];*)
lineM=2*r[[2]]/f[[2]]*x;
picl=Plot[lineM,{x,0,0.01},PlotStyle->{Red}];
pic=Show[picp,picl];
picf=ListPlot[pf,
Frame->True,
PlotStyle->{Blue},
FrameLabel->{Style["m (GeV)",Bold,12.5],Style["f_\[Pi] (GeV)",Bold,12.5]}];

picr=ListPlot[pr,
Frame->True,
PlotStyle->{Blue},
FrameLabel->{Style["m (GeV)",Bold,12.5],Style["r_P (GeV^2)",Bold,12.5]}];
Export["pm-M^2.png",pic];
Export["pm-f.png",picf];
Export["pm-r.png",picr];
pd1=Table[{m[[i]],2*r[[i]]/f[[i]]},{i,1,Length[m]}];
pd2=Table[{m[[i]],M[[i]]^2/m[[i]]},{i,1,Length[m]}];
picd=ListPlot[{pd1,pd2},
Frame->True,
PlotStyle->{Blue,Red},
FrameLabel->{Style["m (GeV)",Bold,12.5],Style["2*r_P/f_\[Pi] (GeV)",Bold,12.5]}];
Export["pm-d.png",picd];

