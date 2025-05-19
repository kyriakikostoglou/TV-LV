function kernels=CREATEKERN(Q,lag,c1_1,c2_1)
%1st-order kernel
k1_1= c1_1*lag{1};
kernels.k1_1=k1_1;
kernels.k2_1=[];

if(Q==2)    
    %2nd-order kernel
    k2_1=c2_1*lag{2}';    
    kernels.k2_1=k2_1;    
end






