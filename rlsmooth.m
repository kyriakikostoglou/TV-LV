function [thms,yhat]=rlsmooth(thm,V,yhat,lam)
%----RLS Smoothing---------------------------------------------------------
N=length(yhat);
thms=thm;
for k=N-1:-1:1
        th=thm(:,k)+lam(k)*(thms(:,k+1)-thm(:,k));
        thms(:,k)=th;
        phit= V(k,:);
        yh=phit*th;
        yhat(k,1)=yh;
end
thms=thms';
end


