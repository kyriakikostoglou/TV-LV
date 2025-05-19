function CreateSimulation(inp1,Q,amp1,amp2,w1,w2,N,SNR,Kerlen)

%--------------Create sinusoidally varying SME kernels--------------------
k1_1ct=exp(-(1:Kerlen)/3).*sin(pi*(1:Kerlen)/5); %The basis SME kernel
[X,Y]=meshgrid(1:Kerlen);
M1=X(:);
M2=Y(:);
%The amplitude of the basis SME kernel is modulated by a sinusoidal signal 
for n=1:round(N/2)
    Realkernels.k1_1(n,:)=(amp1*sin(2*pi*(w1/N)*n)+1)*k1_1ct;
    if(Q==2)
        Realkernels.k2_1(n,:)=Realkernels.k1_1(n,M1).*Realkernels.k1_1(n,M2);
    end
end
for n=round(N/2)+1:N
    Realkernels.k1_1(n,:)=(amp2*sin(2*pi*(w2/N)*n)-1)*k1_1ct;
    if(Q==2)
        Realkernels.k2_1(t,:)=Realkernels.k1_1(n,M1).*Realkernels.k1_1(n,M2);
    end
end        

INP1=flipud(buffer(inp1,Kerlen,Kerlen-1));
v1=zeros(N,1);
out_noisefree=zeros(N,1);

for n=1:N
    v1(n)=INP1(:,n)'*Realkernels.k1_1(n,:)';
    if(Q==1)
        out_noisefree(n,1)=v1(n);
    elseif(Q==2)
        out_noisefree(n,1)=v1(n)+(v1(n))^2;
        %out_noisefree(n)=v1(n)+INP1(:,n)'*(reshape((Realkernels.k2_1(n,:))',Kerlen,Kerlen))*INP1(:,n);
    end

end   
%-------------Add noise to system output-----------------------------------
rng('shuffle')
out=awgn(out_noisefree,SNR,'measured');
%-------------Save---------------------------------------------------------
ff=sprintf('SIM_w%d_w%d_amp%1.2f_amp%1.2f.mat',w1,w2,amp1,amp2);
save(ff,'Realkernels','inp1','out_noisefree','out','SNR','amp1','amp2','w1','w2','Q')
end
