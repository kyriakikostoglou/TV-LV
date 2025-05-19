function V=GetV(inp1,L,Q,a,T)

beta=1-a;
sbeta=sqrt(beta);
salpha=sqrt(a); 
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

% Vmat=ones(N,1);
Vmat=[];
Vmat=[Vmat v_ev_1'];
Vmat2_1=[];

if (Q==2)
    for k1=1:L
        for k2=k1:L
            v_ev_2_1=(v_ev_1(k1,:)).*v_ev_1(k2,:);
            Vmat2_1=[Vmat2_1 v_ev_2_1'];
        end
    end
    Vmat=[Vmat Vmat2_1];
end
V=Vmat;
 
  


end
