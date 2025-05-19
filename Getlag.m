function lagMat=Getlag(L,Q,a,Kerlen,T)

beta=1-a;
salpha=sqrt(a);
lagMat=cell(2,1);
lag1=zeros(L,Kerlen);

for n=1:Kerlen
    lag1(1,n)=T*sqrt((a^(n-1))*beta);
    for k=2:L
        lag1(k,1)=salpha*lag1(k-1,1);
        if (n>1)
            lag1(k,n)=salpha*lag1(k,n-1)+salpha*lag1(k-1,n)-lag1(k-1,n-1);       
        end
    end
end

lagMat{1}=lag1;

if(Q==2)
    lagMat{2} = kron(lag1',lag1'); 
end

end