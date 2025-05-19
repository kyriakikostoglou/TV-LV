function thm2=makesymcoef(L,thm,thm2)
%------Make symmetric coefficients if Q=2----------------------------------
for k=1:size(thm,1)
    c2_1=zeros(L);
    ka=0;     
    for ind1=1:L
        for ind2=ind1:L
            ka=ka+1;
            c2_1(ind1,ind2)=thm(k,L+ka);            
        end
    end
    for ind1=1:L
        for ind2=1:ind1-1
            c2_1(ind1,ind2)=c2_1(ind2,ind1);           
        end
    end
    c2_1 = 0.5*c2_1 + 0.5*tril(triu(c2_1));
    thm2(k,:)= [thm(k,1:L)';c2_1(:)]'; 
end
