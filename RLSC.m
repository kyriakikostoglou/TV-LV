function J=RLSC(out,inp1,Q,X,metric,ignore,T)
L=X(1);
a=X(2);
lam=X(3);
pin=X(4);

if(Q==1)
    L_len=0;
    totPar = L;
else    
    L_len = sum(L:-1:1);
    totPar = L+L_len;
end

%----Built V matrix--------------------------------------------------------
beta        = 1-a;
sbeta       = sqrt(beta);
salpha      = sqrt(a); 
N=length(inp1);
v1=zeros(L,N);
v1(1,1)=T*sbeta*inp1(1);

for i=1:N
    if(i>1)
        v1(1,i)=salpha*v1(1,i-1)+T*sbeta*inp1(i);
    end
    for k=2:L
        v1(k,1)=salpha*v1(k-1,1);
        if(i>1)
            v1(k,i)=salpha*(v1(k,i-1)+v1(k-1,i))-v1(k-1,i-1);
        end
    end
end

v_ev_1=v1;
Vmat=[];
Vmat=[Vmat v_ev_1'];
Vmat2_1=[];

if (Q>=2)
    for k1=1:L
        for k2=k1:L
            v_ev_2_1=(v_ev_1(k1,:)).*v_ev_1(k2,:);
            Vmat2_1=[Vmat2_1 v_ev_2_1'];
        end
    end
    Vmat=[Vmat Vmat2_1];
end
V=Vmat;

%----RLSC recursive estimator----------------------------------------------
yhat = zeros(N,1);
th = zeros(totPar,1);        
p = pin*eye(totPar);
for k=1:N
    phit= V(k,:);
    phi=phit';
    yh=phit*th;
    epsi=out(k,1)-yh;
    pphi=p*phi;
    rt=phit*pphi;
    K=pphi/(lam+rt);
    W=(p-(K*(phit*p)));
    p=(1/lam)*W;
    th=(th+(K*epsi));    
    yhat(k,1)=yh;          
end

e = out-yhat;
if(metric==1)     %BIC
    n=length(e(ignore:end));
    R=(norm(e(ignore:end)))^2;
    J=0.5*n*log(R/n)+0.5*(totPar)*log(n);
elseif(metric==2) %AIC
    n=length(e(ignore:end));
    R=(norm(e(ignore:end)))^2;
    J=(0.5*n*log(R/n))+totPar;
end

end