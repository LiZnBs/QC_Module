CalApBp[ au0_, bu0_, ad0_, bd0_, ar0_, br0_, kp0_, zp0_]:=Module[
    {
    s1=Table[0,{Length[zp0[[1]]]}],
    s2=Table[0,{Length[zp0[[1]]]}],
    (*Total momentum*)
    P={0,0,0,I M}
    },
    (*Because of the symmetry of y, we do not need to calculate y*)
    (*Q[q_,y_,z_]:={0, Sqrt[q(1-z^2)(1-y^2)], y Sqrt[q(1-z^2)], Sqrt[q] z};*)
    Q[q_,z_]:={0, 0, Sqrt[q(1-z^2)], Sqrt[q] z};
    CloseKernels[];
    LaunchKernels[kn];
    SetSharedVariable[s1,s2];
    ParallelDo[ 
    temp1=Table[0,{Length[kp0[[1]] ]}];
    temp2=Table[0,{Length[kp0[[1]] ]}];
        Do[
            temp1[[j]]=(
                Sum[
                au0[[k]]
                *kp0[[2]][[k]]*(-1-I*\[Eta]*M/Sqrt[kp0[[1]][[k]] ])*kp0[[1]][[k]] *Log[\[Lambda] ]
                /((Q[kp0[[1]][[k]],1 ]+\[Eta] P).(Q[kp0[[1]][[k]],1 ]+\[Eta] P)
                -(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P).(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P)) 
            ,{k,1,Length[kp0[[1]] ]}
                ]+Sum[
                ad0[[k]] 
                *kp0[[2]][[k]]*(+1-I*\[Eta]*M/Sqrt[kp0[[1]][[k]] ])*kp0[[1]][[k]]*Log[\[Lambda] ]
                /((Q[kp0[[1]][[k]] ,-1]+\[Eta] P).(Q[kp0[[1]][[k]] ,-1]+\[Eta] P)
                -(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P).(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P))  
            ,{k,1,Length[kp0[[1]] ]}
                ]+Sum[  
                ar0[[k]] 
                *zp0[[2]][[k]]*I*M*Sqrt[\[Lambda] ] 
                /((Q[\[Lambda],zp0[[1]][[k]] ]+\[Eta] P).(Q[\[Lambda],zp0[[1]][[k]] ]+\[Eta] P)
                -(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P).(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P))   
                ,{k,1,Length[zp0[[1]] ]} 
                ]  
            )
            /
            (
            Sum[
                kp0[[2]][[k]]*(-1-I*\[Eta]*M/Sqrt[kp0[[1]][[k]] ])*kp0[[1]][[k]] *Log[\[Lambda] ]
                /((Q[kp0[[1]][[k]],1 ]+\[Eta] P).(Q[kp0[[1]][[k]],1 ]+\[Eta] P)  
                -(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P).(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P)) 
            ,{k,1,Length[kp0[[1]] ]}
                ]+Sum[
                kp0[[2]][[k]]*(+1-I*\[Eta]*M/Sqrt[kp0[[1]][[k]] ])*kp0[[1]][[k]]*Log[\[Lambda] ]
                /((Q[kp0[[1]][[k]] ,-1]+\[Eta] P).(Q[kp0[[1]][[k]] ,-1]+\[Eta] P)
                -(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P).(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P))  
            ,{k,1,Length[kp0[[1]] ]}
                ]+Sum[  
                zp0[[2]][[k]]*I*M*Sqrt[\[Lambda] ] 
                /((Q[\[Lambda],zp0[[1]][[k]] ]+\[Eta] P).(Q[\[Lambda],zp0[[1]][[k]] ]+\[Eta] P) 
                -(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P).(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P))   
                ,{k,1,Length[zp0[[1]] ]} 
                ]  
            ); 
            temp2[[j]]=(
                Sum[
                bu0[[k]]
                 *kp0[[2]][[k]]*(-1-I*\[Eta]*M/Sqrt[kp0[[1]][[k]] ])*kp0[[1]][[k]] *Log[\[Lambda] ]
                /((Q[kp0[[1]][[k]],1 ]+\[Eta] P).(Q[kp0[[1]][[k]],1 ]+\[Eta] P)  
                -(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P).(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P)) 
            ,{k,1,Length[kp0[[1]] ]}
                ]+Sum[
                bd0[[k]]
                *kp0[[2]][[k]]*(+1-I*\[Eta]*M/Sqrt[kp0[[1]][[k]] ])*kp0[[1]][[k]]*Log[\[Lambda] ]
                /((Q[kp0[[1]][[k]] ,-1]+\[Eta] P).(Q[kp0[[1]][[k]] ,-1]+\[Eta] P)
                -(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P).(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P))  
            ,{k,1,Length[kp0[[1]] ]}
                ]+Sum[  
                br0[[k]]
                *zp0[[2]][[k]]*I*M*Sqrt[\[Lambda] ] 
                /((Q[\[Lambda],zp0[[1]][[k]] ]+\[Eta] P).(Q[\[Lambda],zp0[[1]][[k]] ]+\[Eta] P)
                -(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P).(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P))   
                ,{k,1,Length[zp0[[1]] ]} 
                ]  
            )
            /
            (
            Sum[
                kp0[[2]][[k]]*(-1-I*\[Eta]*M/Sqrt[kp0[[1]][[k]] ])*kp0[[1]][[k]] *Log[\[Lambda] ]
                /((Q[kp0[[1]][[k]],1 ]+\[Eta] P).(Q[kp0[[1]][[k]],1 ]+\[Eta] P)  
                -(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P).(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P)) 
            ,{k,1,Length[kp0[[1]] ]}
                ]+Sum[
                kp0[[2]][[k]]*(+1-I*\[Eta]*M/Sqrt[kp0[[1]][[k]] ])*kp0[[1]][[k]]*Log[\[Lambda] ]
                /((Q[kp0[[1]][[k]] ,-1]+\[Eta] P).(Q[kp0[[1]][[k]] ,-1]+\[Eta] P)
                -(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P).(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P))  
            ,{k,1,Length[kp0[[1]] ]}
                ]+Sum[  
                zp0[[2]][[k]]*I*M*Sqrt[\[Lambda] ] 
                /((Q[\[Lambda],zp0[[1]][[k]] ]+\[Eta] P).(Q[\[Lambda],zp0[[1]][[k]] ]+\[Eta] P)
                -(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P).(Q[kp0[[1]][[j]], zp0[[1]][[i]] ]+\[Eta] P))   
                ,{k,1,Length[zp0[[1]] ]} 
                ]  
            )
            ,{j,1,Length[kp0[[1]] ]}
        ];
        s1[[i]]=temp1;
        s2[[i]]=temp2;
        ,{i,1,Length[zp0[[1]] ]}
        
    ];
(*{{{zp1,kp1},{zp1,kp2}},{{zp2,kp1},{zp2,kp2}}}*)
{s1,s2}
]