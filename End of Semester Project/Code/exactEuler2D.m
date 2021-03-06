function [w]=exactEuler2D(wL,wR,gamma,x,t,Qtol)
%exactEuler2D solves the Riemann problem exact for the 2D
%Euler equations in primitive variables with ideal gas
%constituive law
%   Input:
%       wL=left initial data
%       wR=right initial data
%       gamma=ideal gas constant
%       x=domain
%       t=time
%       Qtol=user specified tolerance
%   Output:
%       w=solution at time t

%% Define primitive variable initial states
rhoL=wL(1);
rhoR=wR(1);
tvL=wL(2);
tvR=wR(2);
nvL=wL(3);
nvR=wR(3);
pL=wL(4);
pR=wR(4);

%Sound Speeds
if pL/rhoL<0
    fprintf('Imaginary Sound Speed\n');
    fprintf('p_l=%d,rho_l=%d\n',pL,rhoL);
    return
elseif  pR/rhoR<0
    fprintf('Imaginary Sound Speed\n');
    fprintf('p_r=%d,rho_r=%d\n',pR,rhoR)
    return
end
cL=sqrt(gamma*pL/rhoL);
cR=sqrt(gamma*pR/rhoR);

%Check Pressure positivity:
if 2*(cL+cR)/(gamma-1)<(tvR-tvL)
    fprintf('Vacuum is generated by data\n');
    return
end

%% Initial Guess
lmax=10;
tol=1e-6;
%Linearized Guess
p0Lin=.5*(pL+pR)-(tvR-tvL)*(rhoL+rhoR)*(cL+cR)/8;
p0Lin=max([p0Lin,tol]);
%Double Rarefaction guess:
pRare=((cL+cR-.5*(gamma-1)*(tvR-tvL))/(cL/pL^((gamma-1)/...
    (2*gamma))+cR/pR^((gamma-1)/(2*gamma))))^(2*gamma/(gamma-1));
%Double Shock guess:
gL=(2/((gamma+1)*rhoL)/(p0Lin+(gamma-1)/(gamma+1)*pL))^.5;
gR=(2/((gamma+1)*rhoR)/(p0Lin+(gamma-1)/(gamma+1)*pR))^.5;

p0Shock=(pL*gL+pR*gR-(tvR-tvL))/(gL+gR);
p0Shock=max([p0Shock,tol]);
% Determine Best Initial guess:
pmax=max([pRare,p0Shock,p0Lin]);
pmin=min([pRare,p0Shock,p0Lin]);
Q=pmax/pmin;
pStar=p0Lin;
if Q<Qtol && pmin<pStar&& pStar<pmax
    p=p0Lin;
else if pStar<=pmin
        p=pRare;
    else
        p=p0Shock;
    end
end

%% Newton Iteration to find middle (star) state:

for l=1:lmax
    [lwave,locusL,derivL]=fk(p,pL,rhoL,gamma);
    [rwave,locusR,derivR]=fk(p,pR,rhoR,gamma);
    pNew=p-(locusL+locusR+tvR-tvL)/(derivL+derivR);
    CHA=2*abs((pNew-p)/(pNew+p));
    p=pNew;
    if CHA<Qtol
        fprintf('Newton Iteration Converged with %f iterations \n', l);
        break
    end
end

%% Construct rest of primitive variables and sound speeds:
pM = p;
[rhoML,rhoMR]=rhofinder(rhoL,rhoR,pM,pL,pR,gamma,lwave,rwave);
tvM= .5*(tvL+tvR)+.5*(locusR-locusL);
cML=sqrt(gamma*pM/rhoML);
cMR=sqrt(gamma*pM/rhoMR);
%% Construct Solution
w=zeros(4,numel(x));
for j=1:numel(x);
    S=x(j)/t;
if S<tvM
    %left
    if pL<pM
        %Shock
        Sl=tvL-cL*((gamma+1)/...
            (2*gamma)*pM/pL+(gamma-1)/(2*gamma))^.5;
        if S<Sl
            w(:,j)=[rhoL;tvL;nvL;pL];
        else
            w(:,j)=[rhoML;tvM;nvL;pM];
        end
    else
        %Rarefaction
        Shl=tvL-cL;
        Stl=tvM-cML;
        if S<Shl
            w(:,j)=[rhoL;tvL;nvL;pL];
        elseif S>Stl
            w(:,j)=[rhoML;tvM;nvL;pM];
        else
            [rho,v,p]=rareleft(rhoL,tvL,pL,cL,gamma,S);
            w(:,j)=[rho;v;nvL;p];
        end
    end
else
    if pR<pM
        %Shock
        Sr=tvR+cR*((gamma+1)/...
            (2*gamma)*pM/pR+(gamma-1)/(2*gamma))^.5;
        if S>Sr
            w(:,j)=[rhoR;tvR;nvR;pR];
        else
            w(:,j)=[rhoMR;tvM;nvR;pM];
        end
    else
        %Rarefaction
        Shr=tvR+cR;
        Str=tvM+cMR;
        if S>Shr
            w(:,j)=[rhoR;tvR;nvR;pR];
        elseif S<Str
            w(:,j)=[rhoMR;tvM;nvR;pM];
        else
            [rho,v,p]=rareright(rhoR,tvR,pR,cR,gamma,S);
            w(:,j)=[rho;v;nvR;p];
        end
    end
end
end


%% Auxillary Functions for Riemann solver
%Determines expressions for the newton iteration
    function [kWave,locusK,derivK]=fk(p,pK,rhoK,gamma)
        if pK>=p %Rarefaction
            locusK=2*sqrt(gamma*pK/rhoK)/(gamma-1)*...
                ((p/pK)^((gamma-1)/(2*gamma))-1);
            derivK=((p/pK)^((gamma-1)/(2*gamma))...
                *sqrt(gamma*pK/rhoK))/(p*gamma);
            kWave=0;
        else %Shock
            locusK=(p-pK)*(2/((gamma+1)*rhoK)/...
                (p+(gamma-1)/(gamma+1)*pK))^.5;
            derivK=sqrt(1/(2*rhoK))*(p*(1+gamma)+pK*(3*gamma-1))...
                /((p*(1+gamma)+pK*(gamma-1))^1.5);
            kWave=1;
        end
    end

%Finds middle rho states
    function [rhoML,rhoMR]=rhofinder(rhoL,rhoR,p,pL,pR,gamma,lWave,rWave)
        if lWave==0
            rhoML=rhoL*(p/pL)^(1/gamma);
            if rWave==0
                rhoMR=rhoR*(p/pR)^(1/gamma);
            else
                rhoMR=rhoR*(p/pR+(gamma-1)/(gamma+1))/...
                    ((gamma-1)/(gamma+1)*(p/pR)+1);
            end
        else
            rhoML=rhoL*(p/pL+(gamma-1)/(gamma+1))/...
                ((gamma-1)/(gamma+1)*(p/pL)+1);
            if rWave==0
                rhoMR=rhoR*(p/pR)^(1/gamma);
            else
                rhoMR=rhoR*(p/pR+(gamma-1)/(gamma+1))/...
                    ((gamma-1)/(gamma+1)*(p/pR)+1);
            end
        end
    end
% Defines solution in left rarefaction
    function [rho,v,p]=rareleft(rhoL,tvL,pL,cL,gamma,S)
        rho=rhoL*(2/(gamma+1)+(gamma-1)/((gamma+1)*cL)...
            *(tvL-S))^(2/(gamma-1));
        v=2/(gamma+1)*(cL+(gamma-1)/2*tvL+S);
        p=pL*(2/(gamma+1)+(gamma-1)/((gamma+1)*cL)...
            *(tvL-S))^(2*gamma/(gamma-1));
    end
%Defines solution in right rarefaction
    function [rho,v,p]=rareright(rhoR,tvR,pR,cR,gamma,S)
        rho=rhoR*(2/(gamma+1)-(gamma-1)/((gamma+1)*cR)...
            *(tvR-S))^(2/(gamma-1));
        v=2/(gamma+1)*(-cR+(gamma-1)/2*tvR+S);
        p=pR*(2/(gamma+1)-(gamma-1)/((gamma+1)*cR)...
            *(tvR-S))^(2*gamma/(gamma-1));
    end

end

