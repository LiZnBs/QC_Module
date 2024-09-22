(*Necessary input: k=k^2, Output: G*)
(*Main functions*)
F[k_]:=(1- E^(-k/(4 mt^2)))/k;
Fuv[k_]:=(8\[Pi]^2 \[Gamma]m F[k])/Log[\[Tau]+(k/\[Lambda]QCD^2+1)^2];
G[k_]:=d 8\[Pi]^2/(\[Omega]^4) E^(-k/\[Omega]^2) +Fuv[k];


(*Print[G[200]];*)