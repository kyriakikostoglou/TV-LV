function [thms,yhat]=kalsmooth(thm,V,PP,PPo,yhat)
%------Kalman Smoothing----------------------------------------------------
N=length(yhat);
thms=thm;
for k=N-1:-1:1
        A=PP{k}/PPo{k+1};
        th=thm(:,k)+A*(thms(:,k+1)-thm(:,k));
        thms(:,k)=th;
        phit= V(k,:);
        yh=phit*th;
        yhat(k,1)=yh;
end
thms=thms';
end


