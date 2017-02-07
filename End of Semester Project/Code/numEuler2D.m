function [err,data,pdata]=numEuler2D(cuInit,gamma,...
    dx,dy,cfl,tF,Case,Order,iBound)
%numEuler2D solves the 2D Euler equations, with ideal gas
%constitutive law, numerically given some initial data
%   Inputs:
%       uInit=4xM matrix of initial data
%       gamma=ideal gas constant
%       cfl=cfl number
%       tF=final time
%       Case=type of numerical method
%       Order=Order of numerical method
%   Outputs:
%       data=primitive data saved at requested snapshot times
%       time=time(s) at which data was saved

%% Initialize Variables

%Give large bound on time for loop to end
dtMax=1e-1;
dtMin=1e-6;
nMax=tF/dtMin;

% Find total number of grid points including ghost cells
M=size(cuInit,1);
N=size(cuInit,2);

ng=Order;
Mtot=M+2*ng;
Ntot=N+2*ng;
ja=ng+1;
jb=Mtot-1;
ka=ng+1;
kb=Ntot-1;

%Initialize vectors
data=zeros(M,N,4);
xFluxMat=zeros(Mtot-1,N,4);
yFluxMat=zeros(M,Ntot-1,4);
timeNow=0;
cu=zeros(Mtot,Ntot,4);

%For twilight zone, define x and y vectors
x = linspace( -2+dx/2-ng*dx, 2-dx/2+ng*dx, Mtot );
y = linspace( -2+dy/2-ng*dy, 2-dy/2+ng*dy, Ntot );
%% Set Initial Conditions
cu(ja:jb,ka:kb,:)=cuInit;
%Set BCs on initial conditions
[cu]=setBCs(cu,ja,jb,ka,kb,Order,iBound,0,x,y,gamma);
%Calculate initial dt
[prim]=convertVar(cu,gamma,Mtot,Ntot);
xLamMax=1e-6;
yLamMax=1e-6;
for j=(ja-1):jb
    for k=(ka-1):kb
        [~,xLam]=xnumFlux(Case,prim(j,k,:),...
            prim(j+1,k,:),gamma);
        %fprintf('%d,%d\n',j,k)
        [~,yLam]=ynumFlux(Case,prim(j,k,:),...
            prim(j,k+1,:),gamma);
        if abs(xLam)>xLamMax
            xLamMax=abs(xLam);
        end
        if abs(yLam)>yLamMax
            yLamMax=abs(yLam);
        end
    end
end
%% Main Loop
for n=1:nMax
    %Convert conserved variables to primitive variables
    [prim]=convertVar(cu,gamma,Mtot,Ntot);
    
    %set dt
    dt=dtMax;
    dt=min(dt,cfl*dx*dy/(xLamMax*dy+yLamMax*dx));
    %Adjust time step to end at tF
    if timeNow+dt>tF
        dt=tF-timeNow;
    end
    %Compute Fluxes
    xLamMax=1e-6;
    yLamMax=1e-6;
    for j=(ja-1):jb
        for k=(ka-1):kb
            [fluxOut,xLam]=xnumFlux(Case,prim(j,k,:),...
                prim(j+1,k,:),gamma);
            %fprintf('%d,%d\n',j,k);
            xFluxMat(j,k,:)=fluxOut;
            [fluxOut,yLam]=ynumFlux(Case,prim(j,k,:),...
                prim(j,k+1,:),gamma);
            yFluxMat(j,k,:)=fluxOut;
            if abs(xLam)>xLamMax
                xLamMax=abs(xLam);
            end
            if abs(yLam)>yLamMax
                yLamMax=abs(yLam);
            end
        end
    end
    %Update Solution
    for j=ja:jb
        xH=(x(j)+x(j+1))/2;
        for k=ka:kb
            
            yH=(y(k)+y(k+1))/2;
            H=findH(xH,yH,timeNow,gamma);
            cu(j,k,:)=squeeze(cu(j,k,:)-...
                (dt/dx)*(xFluxMat(j,k,:)-xFluxMat(j-1,k,:))-...
                (dt/dy)*(yFluxMat(j,k,:)-yFluxMat(j,k-1,:)))+...
                dt*H';
        end
    end
    
    %Calculate time
    timeNow=timeNow+dt;
    %Set BCs
    cu=setBCs(cu,ja,jb,ka,kb,Order,iBound,timeNow,x,y,gamma);
    %calculate exact solution
    exact=zeros(M,N,4);
    for j=1:Mtot
        for k=1:Ntot
            exact(j,k,:)=getExCons(x(j),y(k),timeNow,gamma);
        end
    end
    
%     err=cu-exact;
%     figure(1)
%     surf(x,y,err(:,:,1))
%     
%     figure(2)
%     surf(x,y,err(:,:,2))
%     
%     figure(3)
%     surf(x,y,err(:,:,3))
%     
%     figure(4)
%     surf(x,y,err(:,:,4))
%     
%     drawnow
%     pause
    %Check if run has completed and save solution
    if timeNow>tF-eps
        [prim]=convertVar(cu,gamma,Mtot,Ntot);
        pdata=prim(ja:jb,ka:kb,:);
        data=cu(ja:jb,ka:kb,:);
        err=computeErr(pdata,M,N,x,y,timeNow,dx,dy);
        break
    end
end
%% Auxillary Functions

%% Bounday Conditions
    function [u]=setBCs(u,ja,jb,ka,kb,Order,iBound,timeNow,x,y,gamma)
        if Order==1
            if iBound==1 %Transmissive everywhere
                for k=ka:kb %#ok<*FXUP>
                    u(ja-1,k,:)=u(ja+1,k,:);
                    u(jb+1,k,:)=u(jb-1,k,:);
                end
                for j=ja:jb
                    u(j,ka-1,:)=u(j,ka+1,:);
                    u(j,kb+1,:)=u(j,kb-1,:);
                end
                u(ja-1,ka-1,:)=u(ja+1,ka+1,:);
                u(jb+1,kb+1,:)=u(jb-1,kb-1,:);
                u(ja-1,kb+1,:)=u(ja+1,kb-1,:);
                u(jb+1,ka-1,:)=u(jb-1,ka+1,:);
            elseif iBound==2 %Reflective everywhere
                for k=ka:kb
                    u(ja-1,k,:)=u(ja,k,:);
                    u(ja-1,k,2)=-u(ja,k,2);
                    u(ja-1,k,3)=-u(ja,k,3);
                    u(jb+1,k,:)=u(jb,k,:);
                    u(jb+1,k,2)=-u(jb,k,2);
                    u(jb+1,k,3)=-u(jb,k,3);
                end
                for j=ja:jb
                    u(j,ka-1,:)=u(j,ka,:);
                    u(j,ka-1,2)=-u(j,ka,2);
                    u(j,ka-1,3)=-u(j,ka,3);
                    u(j,kb+1,:)=u(j,kb,:);
                    u(j,kb+1,2)=-u(j,kb,2);
                    u(j,kb+1,3)=-u(j,kb,3);
                end
                u(ja-1,ka-1,:)=u(ja,ka,:);
                u(ja-1,ka-1,2)=-u(ja,ka,2);
                u(ja-1,ka-1,3)=-u(ja,ka,3);
                
                u(jb+1,ka-1,:)=u(jb,ka,:);
                u(jb+1,ka-1,2)=-u(jb,ka,2);
                u(jb+1,ka-1,3)=-u(jb,ka,3);
                
                u(ja-1,kb+1,:)=u(ja,kb,:);
                u(ja-1,kb+1,2)=-u(ja,kb,2);
                u(ja-1,kb+1,3)=-u(ja,kb,3);
                
                u(jb+1,kb+1,:)=u(jb,kb,:);
                u(jb+1,kb+1,2)=-u(jb,kb,2);
                u(jb+1,kb+1,3)=-u(jb,kb,3);
                
            elseif iBound==3 %Transmissive horiz,
                %Reflective vert
                for k=ka:kb %#ok<*FXUP>
                    u(ja-1,k,:)=u(ja+1,k,:);
                    u(jb+1,k,:)=u(jb-1,k,:);
                end
                for j=ja:jb
                    u(j,ka-1,:)=u(j,ka,:);
                    u(j,ka-1,2)=-u(j,ka,2);
                    u(j,kb+1,:)=u(j,kb,:);
                    u(j,kb+1,2)=-u(j,kb,2);
                end
            elseif iBound==4 %Reflective left,top,
                %Transmissive bot,right
                for k=ka:kb
                    u(ja-1,k,:)=u(ja,k,:);
                    u(ja-1,k,2)=-u(ja,k,2);
                    u(jb+1,k,:)=u(jb-1,k,:);
                end
                for j=ja:jb
                    u(j,ka-1,:)=u(j,ka+1,:);
                    u(j,kb+1,:)=u(j,kb,:);
                    u(j,kb+1,2)=-u(j,kb,2);
                end
            elseif iBound==5 %Reflective left
                %Transmissive else
                for k=ka:kb
                    u(ja-1,k,:)=u(ja,k,:);
                    u(ja-1,k,2)=-u(ja,k,2);
                    u(jb+1,k,:)=u(jb-1,k,:);
                end
                for j=ja:jb
                    u(j,ka-1,:)=u(j,ka+1,:);
                    u(j,kb+1,:)=u(j,kb-1,:);
                end
                
            elseif iBound==6 %Transmissive left
                %Reflective else
                for k=ka:kb
                    u(ja-1,k,:)=u(ja+1,k,:);
                    u(jb+1,k,:)=u(jb-1,k,:);
                    u(jb+1,k,2)=-u(jb-1,k,2);
                end
                for j=ja:jb
                    u(j,ka-1,:)=u(j,ka,:);
                    u(j,ka-1,2)=-u(j,ka,2);
                    u(j,kb+1,:)=u(j,kb,:);
                    u(j,kb+1,2)=-u(j,kb,2);
                end
                u(ja-1,ka-1,:)=u(ja,ka,:);
                u(ja-1,ka-1,2)=u(ja,ka,2);
                u(ja-1,ka-1,3)=-u(ja,ka,3);
                
                u(jb+1,ka-1,:)=u(jb,ka,:);
                u(jb+1,ka-1,2)=-u(jb,ka,2);
                u(jb+1,ka-1,3)=-u(jb,ka,3);
                
                u(ja-1,kb+1,:)=u(ja,kb,:);
                u(ja-1,kb+1,2)=u(ja,kb,2);
                u(ja-1,kb+1,3)=-u(ja,kb,3);
                
                u(jb+1,kb+1,:)=u(jb,kb,:);
                u(jb+1,kb+1,2)=-u(jb,kb,2);
                u(jb+1,kb+1,3)=-u(jb,kb,3);
            elseif iBound==7 %Twilight Zone
                for k=ka:kb %#ok<*FXUP>
                    exactL=getExCons(x(ja-1),y(k),timeNow,gamma);
                    exactR=getExCons(x(jb+1),y(k),timeNow,gamma);
                    u(ja-1,k,:)=exactL;
                    u(jb+1,k,:)=exactR;
                    u(ja,k,:)=getExCons(x(ja),y(k),timeNow,gamma);
                    u(jb,k,:)=getExCons(x(jb),y(k),timeNow,gamma);
                end
                for j=ja:jb
                    exactB=getExCons(x(j),y(ka-1),timeNow,gamma);
                    exactT=getExCons(x(j),y(kb+1),timeNow,gamma);
                    u(j,ka-1,:)=exactB;
                    u(j,kb+1,:)=exactT;
                    u(j,ka,:)=getExCons(x(j),y(ka),timeNow,gamma);
                    u(j,kb,:)=getExCons(x(j),y(kb),timeNow,gamma);
                    
                end
                u(ja-1,ka-1,:)=getExCons(x(ja-1),y(ka-1),timeNow,gamma);
                u(jb+1,kb+1,:)=getExCons(x(jb+1),y(kb+1),timeNow,gamma);
                u(ja-1,kb+1,:)=getExCons(x(ja-1),y(kb+1),timeNow,gamma);
                u(jb+1,ka-1,:)=getExCons(x(jb+1),y(ka-1),timeNow,gamma);
            end
        end
    end
%% Primitive Variable Converter
    function [w]=convertVar(u,gamma,Mtot,Ntot)
        %density
        rho=u(:,:,1);
        %momenta
        mx=u(:,:,2);
        my=u(:,:,3);
        %energy
        E=u(:,:,4);
        w(:,:,1)=rho;
        for i=1:Mtot
            for j=1:Ntot
                %Velocities
                w(i,j,2)=mx(i,j)/rho(i,j);
                w(i,j,3)=my(i,j)/rho(i,j);
                %Pressure
                w(i,j,4)=(gamma-1)*(E(i,j)-...
                    .5*(mx(i,j)^2+my(i,j)^2)/rho(i,j));
            end
        end
    end

%% Twilight Zone
%Source Term
    function H=findH(x,y,t,gamma)
        H(1)=(cos(x)*cos(y)*sin(t))^7*(2*cos(x)*cos(y)*cos(t)-...
            sin(t)*(2*cos(y)*sin(x)+cos(x)*sin(y)));
        H(2)=cos(x)^8*(cos(y)*sin(t))^7*(2*cos(y)*cos(t)-sin(t)*sin(y));
        H(3)=.5*(cos(x)*cos(y)*sin(t))^7*(2*cos(x)*cos(y)*cos(t)+...
            sin(t)*(-2*cos(y)*sin(x)+3*cos(x)*sin(y)));
        H(4)=1.2*cos(t)*cos(x)^8*cos(y)^8*sin(t)^7-...
            (2*cos(t)*cos(x)^8*cos(y)^8*sin(t)^7)/(gamma-1)+...
            .5*(1-.25*cos(x)^8*cos(y)^8*sin(t)^8+...
            (1-.25*cos(x)^8*cos(y)^8*sin(t)^8)/(gamma-1)+...
            5/8*(1+.25*cos(x)^8*cos(y)^8*sin(t)^8))+...
            .75*cos(x)^7*cos(y)^8*sin(t)^8*sin(x)+...
            (2*cos(x)^7*cos(y)^8*sin(t)^8*sin(x))/(gamma-1);
    end
%Exact Solution (primitive variables)
    function exact=getEx(x,y,t)
        exact(1)=1+.25*(cos(x)*cos(y)*sin(t))^8;
        exact(2)=1;
        exact(3)=.5;
        exact(4)=1-.25*(cos(x)*cos(y)*sin(t))^8;
    end

    function exact=getExCons(x,y,t,gamma)
        exact(1)=1+.25*(cos(x)*cos(y)*sin(t))^8;
        exact(2)=exact(1);
        exact(3)=.5*exact(1);
        exact(4)=.5*exact(1)*(1+.25)+...
            (1-.25*(cos(x)*cos(y)*sin(t))^8)/(gamma-1);
    end

%Error Calculator
    function err=computeErr(pdata,M,N,x,y,t,dx,dy)
        err=0;
        for j=1:M
            for k=1:N
                ex=getEx(x(j),y(k),t);
                err=abs(squeeze(pdata(j,k,:))-ex')*dy*dx+err;
                %err=max(abs(squeeze(pdata(j,k,:))-ex'),err);
            end
        end
    end

%% Numerical Flux Functions
    function [fluxOut,xLam]=xnumFlux(Case,A,B,gamma)
        if Case==1
            %Exact
            [z,xLam]=eulerExact(A,B,gamma,1e-8);
            fluxOut=xFlux(z,gamma);
        elseif Case==2
            %HLLC
            [fluxOut,xLam]=eulerHLLC(A,B,gamma);
        end
    end
    function [fluxOut,yLam]=ynumFlux(Case,A,B,gamma)
        
        C(:,:,1)=A(:,:,1);
        C(:,:,2)=A(:,:,3);
        C(:,:,3)=A(:,:,2);
        C(:,:,4)=A(:,:,4);
        D(:,:,1)=B(:,:,1);
        D(:,:,2)=B(:,:,3);
        D(:,:,3)=B(:,:,2);
        D(:,:,4)=B(:,:,4);
        if Case==1
            [z,yLam]=eulerExact(C,D,gamma,1e-8);
            q=z;
            q(2)=z(3);
            q(3)=z(2);
            fluxOut=yFlux(q,gamma);
            
        elseif Case==2
            [z,yLam]=eulerHLLC(C,D,gamma);
            fluxOut=z;
            fluxOut(2)=z(3);
            fluxOut(3)=z(2);
        end
    end

%% Flux Calculators
    function z=xFlux(A,gamma)
        rho=A(1);
        u=A(2);
        v=A(3);
        p=A(4);
        E=p/(gamma-1)+.5*rho*(u^2+v^2);
        z=[rho*u;rho*u^2+p;rho*u*v;u*(E+p)];
        
    end
    function z=yFlux(A,gamma)
        rho=A(1);
        u=A(2);
        v=A(3);
        p=A(4);
        E=p/(gamma-1)+.5*rho*(u^2+v^2);
        z=[rho*v;rho*u*v;rho*v^2+p;v*(E+p)];
    end


%% HLLC Approximate Riemann Solver
    function [fluxOut,lam]=eulerHLLC(wL,wR,gamma)
        rhoL=wL(1);
        rhoR=wR(1);
        tvL=wL(2);
        tvR=wR(2);
        nvL=wL(3);
        nvR=wR(3);
        pL=wL(4);
        pR=wR(4);
        aL=sqrt(gamma*rhoL/pL);
        aR=sqrt(gamma*rhoR/pR);
        EL=.5*rhoL*(tvL^2+nvL^2)+pL/(gamma-1);
        ER=.5*rhoR*(tvR^2+nvR^2)+pR/(gamma-1);
        uL=[rhoL;rhoL*tvL;rhoL*nvL;EL];
        uR=[rhoR;rhoR*tvR;rhoR*nvR;ER];
        
        rhob=.5*(rhoL+rhoR);
        ab=.5*(aL+aR);
        ppvrs=.5*(pL+pR)-.5*(tvR-tvL)*rhob*ab;
        ps=max(0,ppvrs);
        qL=qFind(ps,pL,gamma);
        qR=qFind(ps,pR,gamma);
        SL=tvL-aL*qL;
        SR=tvR+aR*qR;
        SS=(pR-pL+rhoL*tvL*(SL-tvL)-rhoR*tvR*(SR-tvR))/...
            (rhoL*(SL-tvL)-rhoR*(SR-tvR));
        if 0 <= SL
            lam=abs(tvL)+aL;
            fluxOut=xFlux(wL,gamma);
        elseif SL<=0 && 0<=SS
            lam=abs(tvL)+aL;
            uSL=rhoL*(SL-tvL)/(SL-SS)*[1;SS;nvL;...
                EL/rhoL+(SS-tvL)*(SS+pL/(rhoL*(SL-tvL)))];
            FL=xFlux(wL,gamma);
            fluxOut=FL+SL*(uSL-uL);
        elseif SS<=0 && 0<=SR
            lam=abs(tvR)+aR;
            uSR=rhoR*(SR-tvR)/(SR-SS)*[1;SS;nvR;...
                ER/rhoR+(SS-tvR)*(SS+pR/(rhoR*(SR-tvR)))];
            FR=xFlux(wR,gamma);
            fluxOut=FR+SR*(uSR-uR);
        elseif 0>=SR
            lam=abs(tvR)+aR;
            fluxOut=xFlux(wR,gamma);
        end
        
    end

%% HLLC Aux Functions
    function q=qFind(pS,pK,gamma)
        if pS<=pK
            q=1;
        else
            q=sqrt(1+(gamma+1)/(2*gamma)*(pS/pK-1));
        end
    end

%% Exact Riemann Solver
    function [w,lam]=eulerExact(wL,wR,gamma,Qtol)
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
            fprintf('time=%d,j=%d,k=%d\n',timeNow,j,k);
            return
        elseif  pR/rhoR<0
            fprintf('Imaginary Sound Speed\n');
            fprintf('p_r=%d,rho_r=%d\n',pR,rhoR)
            fprintf('time=%d,j=%d,k=%d\n',timeNow,j,k);
            return
        end
        cL=sqrt(gamma*pL/rhoL);
        cR=sqrt(gamma*pR/rhoR);
        
        %Check Pressure positivity:
        if 2*(cL+cR)/(gamma-1)<abs(tvR-tvL)
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
        pRare=((cL+cR-.5*(gamma-1)*(tvR-tvL))/(cL/pL^((gamma-1)/(2*gamma))+...
            cR/pR^((gamma-1)/(2*gamma))))^(2*gamma/(gamma-1));
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
        else if pStar<pmin
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
                %fprintf('Newton Iteration Converged with %f iterations \n', n);
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
        S=0;
        if S<tvM
            %left
            if pL<pM
                %Shock
                Sl=tvL-cL*((gamma+1)/...
                    (2*gamma)*pM/pL+(gamma-1)/(2*gamma))^.5;
                if S<Sl
                    w=[rhoL;tvL;nvL;pL];
                    lam=abs(tvL)+cL;
                else
                    w=[rhoML;tvM;nvL;pM];
                    lam=abs(tvM)+cML;
                end
            else
                %Rarefaction
                Shl=tvL-cL;
                Stl=tvM-cML;
                if S<Shl
                    w=[rhoL;tvL;nvL;pL];
                    lam=abs(tvL)+cL;
                elseif S>Stl
                    w=[rhoML;tvM;nvL;pM];
                    lam=abs(tvM)+cML;
                else
                    [rho,v,p]=rareleft(rhoL,tvL,pL,cL,gamma,S);
                    w=[rho;v;nvL;p];
                    lam=abs(v)+sqrt(gamma*p/rho);
                end
            end
        else
            if pR<pM
                %Shock
                Sr=tvR+cR*((gamma+1)/...
                    (2*gamma)*pM/pR+(gamma-1)/(2*gamma))^.5;
                if S>Sr
                    w=[rhoR;tvR;nvR;pR];
                    lam=abs(tvR)+cR;
                else
                    w=[rhoMR;tvM;nvR;pM];
                    lam=abs(tvM)+cMR;
                end
            else
                %Rarefaction
                Shr=tvR+cR;
                Str=tvM+cMR;
                if S>Shr
                    w=[rhoR;tvR;nvR;pR];
                    lam=abs(tvR)+cR;
                elseif S<Str
                    w=[rhoMR;tvM;nvR;pM];
                    lam=abs(tvM)+cMR;
                else
                    [rho,v,p]=rareright(rhoR,tvR,pR,cR,gamma,S);
                    w=[rho;v;nvR;p];
                    lam=abs(v)+sqrt(gamma*p/rho);
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



