function results=TEST(out,inp1,Q,X,ignore,T,method,Kerlen,smooth)
L=X(1);
a=X(2);
pin=X(end);
if(Q==1)
    L_len=0;
    totPar = L;
    totPar2=0;
else    
    L_len = sum(L:-1:1);
    totPar = L+L_len;
    totPar2=L+L*L;
end
N=length(inp1);

%---Build V matrix and DLF matrix------------------------------------------
V=GetV(inp1,L,Q,a,T);
lag= Getlag(L,Q,a,Kerlen,T);

%---Recursive Estimators--------------------------------------------------
yhat = zeros(N,1);
thm= zeros(N,totPar);
th=zeros(totPar,1);
if(Q==2)
    thm2= zeros(N,totPar2);
end

p = pin*eye(totPar);


if(method==1)
%----RLSC recursive estimator----------------------------------------------
    
    lam=X(3);
    lambda=zeros(N,1);
    lambda(1:end)=lam;
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
        thm(k,:)=th;
        yhat(k,1)=yh; 
        
    end


%----RLSA recursive estimator----------------------------------------------
elseif(method==2)
    
    A=X(3);
    B=X(4);
    lame=X(5);    
    st=0;
    lambda=zeros(N,1);
    for k=1:N
        phit= V(k,:);
        phi=phit';
        yh=phit*th;
        epsi=out(k)-yh;
        ep2=epsi*epsi;  
        pphi=p*phi;
        rt=(phit*(pphi));
        Ct=((rt*ep2)/(st*(1+rt)));
        Ct=max(Ct,0);
        S=gammainc(Ct/2,totPar/2);  
        S=1-S;
        lam1=real(min(max(A,S),B));
        lambda(k)=lam1;
        K=(pphi/(real(lam1)+rt));
        W=(p-(K*(phit*p)));
        p=(1/lam1)*W;
        th=(th+(K*epsi));
        thm(k,:)=th;
        yhat(k,1)=yh;        
        st=lame*st+(1-lame)*ep2;            
    end

%----KF recursive estimator----------------------------------------------
elseif(method==3)

    r1=X(3);
    r2=X(4);
    R1=r1*eye(totPar);
    R2=r2;
    for k=1:N
        phit= V(k,:);
        phi=phit';
        yh=phit*th;
        epsi=out(k,1)-yh;
        pphi=p*phi;
        rt=phit*pphi;
        K=pphi/(R2+rt);
        PPo{k}=p+R1;
        p=p+R1-(K*(phit*p)); 
        PP{k}=p;
        th=(th+(K*epsi));    
        thm(k,:)=th;
        yhat(k,1)=yh;   
        
    end

%----KFA recursive estimator----------------------------------------------
else
    
    lam=X(3);
    r2=X(4);
    r1=0;
    R2=r2;
    for k=1:N
        phit= V(k,:);
        phi=phit';
        yh=phit*th;
        epsi=out(k,1)-yh;
        r1=lam*r1+(1-lam)*epsi^2;
        R1=r1*eye(totPar);    
        pphi=p*phi;
        rt=phit*pphi;
        K=pphi/(R2+rt);
        PPo{k}=p+R1;
        p=p+R1-(K*(phit*p)); 
        PP{k}=p;
        th=(th+(K*epsi));    
        thm(k,:)=th;
        yhat(k,1)=yh;         
    end


end
e = out-yhat;
J=1-var(e(ignore:end))/var(out(ignore:end));

%---Apply smoothing if smooth is set to 1----------------------------------
if((smooth==1)&(method>=3))
        [thm,yhat]=kalsmooth(thm',V,PP,PPo,yhat);
elseif((smooth==1)&(method<3))
        [thm,yhat]=rlsmooth(thm',V,yhat,lambda);
end

if(Q==2)
    thm2=makesymcoef(L,thm,thm2); 
end 

%-----Extract the TV kernels-----------------------------------------------
if(Q==1)
    c1_1=thm(:,1:L);
    C1_1=thm(ignore:end,1:L);
    C2_1=[];
    Coef=c1_1(ignore:end,:);
    k1_1=zeros(size(c1_1,1),size(lag{1},2));
    k1_1= c1_1*lag{1};
    K1_1=k1_1(ignore:end,:);
    K2_1=[];
else
    c1_1=thm2(:,1:L);
    c2_1=thm2(:,L+1:end);    
    C1_1=thm(:,1:L);
    C2_1=thm(:,L+1:end);
    Coef=[C1_1(ignore:end,:) C2_1(ignore:end,:)];
    k1_1=zeros(size(c1_1,1),size(lag{1},2));
    k2_1=zeros(size(c2_1,1),size(lag{2},1));
    k1_1= c1_1*lag{1};
    K1_1=k1_1(ignore:end,:);
    k2_1= c2_1*lag{2}';
    K2_1=k2_1(ignore:end,:);
end
results.thm=thm(ignore:end,:);
results.yhat=yhat(ignore:end);
results.e=e(ignore:end);
results.C1_1=C1_1;
results.C2_1=C2_1;
results.K1_1=K1_1; %Estimated 1st-order Kernel
results.K2_1=K2_1; %Estimated 2nd-order Kernel
results.Coef=Coef; %Estimated Coefficients
results.J=J;       %VAF(%) 











