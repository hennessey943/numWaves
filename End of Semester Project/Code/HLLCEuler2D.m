function [u]=HLLCEuler2D(wL,wR,gamma,x,t)
%define nice variable names
rhoL=wL(1);
rhoR=wR(1);
tvL=wL(2);
tvR=wR(2);
nvL=wL(3);
nvR=wR(3);
pL=wL(4);
pR=wR(4);
aL=sqrt(gamma*pL/rhoL);
aR=sqrt(gamma*pR/rhoR);
EL=.5*rhoL*(tvL^2+nvL^2)+pL/(gamma-1);
ER=.5*rhoR*(tvR^2+nvR^2)+pR/(gamma-1);
uL=[rhoL;rhoL*tvL;rhoL*nvL;EL];
uR=[rhoR;rhoR*tvR;rhoR*nvR;ER];

%Compute Pressure Estimate
rhob=.5*(rhoL+rhoR);
ab=.5*(aL+aR);
ppvrs=.5*(pL+pR)-.5*(tvR-tvL)*rhob*ab;
pS=max(0,ppvrs);

%Compute Wave Speed Estimates
qL=findq(pS,pL,gamma);
qR=findq(pS,pR,gamma);
SL=tvL-aL*qL;
SR=tvR+aR*qR;

SS=(pR-pL+rhoL*tvL*(SL-tvL)-rhoR*tvR*(SR-tvR))/...
    (rhoL*(SL-tvL)-rhoR*(SR-tvR));

%Construct Solution
u=zeros(4,numel(x));
for i=1:numel(x)
    xi=x(i)/t;
    if xi<SL
        u(:,i)=uL;
    elseif SL<=xi && xi<= SS
        u(:,i)=rhoL*(SL-tvL)/(SL-SS)*[1;SS;nvL;...
            EL/rhoL+(SS-tvL)*(SS+pL/(rhoL*(SL-tvL)))];
    elseif SS<= xi && xi<= SR
        u(:,i)=rhoR*(SR-tvR)/(SR-SS)*[1;SS;nvR;...
            ER/rhoR+(SS-tvR)*(SS+pR/(rhoR*(SR-tvR)))];
    elseif xi>=SR
        u(:,i)=uR;
    end
end
end

function qK=findq(pS,pK,gamma)
if pS<=pK
    qK=1;
else
    qK=sqrt(1+(gamma+1)/(2*gamma)*(pS/pK-1));
end
end
