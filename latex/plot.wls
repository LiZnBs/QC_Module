m={0.002,0.0034,0.004,0.005,0.006,0.007,0.008,0.01}
M={0.103,0.134,0.1454,0.163,0.1785,0.1925,0.2072,0.239}
f={0.092536,0.093628,0.094167,0.094902,0.0956966,0.0964723,0.0972156,0.0986266}
p=Table[{m[[i]],M[[i]]^2},{i,1,Length[m]}]
picp=ListLinePlot[p,
Frame->True,
PlotStyle->{Blue},
Epilog -> {Blue, PointSize[Large], Point[p]},
FrameLabel->{Style["m (GeV)",Bold,12.5],Style["M^2 (GeV^2)",Bold,12.5]}]
Fix[pp_] := Module[
        {
        model, p1, p2, f(*parameters*)
        },

    model = p1*x + p2;
    f = FindFit[pp, model, {p1, p2}, x];
    model /. f
    ];
line=Fix[p]
picl=Plot[line,{x,0,0.01},PlotStyle->{Red}]
pic=Show[picp,picl]
Export["pm-M^2.png",pic]
