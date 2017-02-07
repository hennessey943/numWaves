function uInit=eulerInit(Case,xmin,xmax,ymin,ymax,M,N,gamma)
one=ones(M,N);
uInit=zeros(M,N,4);
x=linspace(xmin,xmax,M);
y=linspace(ymin,ymax,N);

if Case ==1
    %Twilight Zone
    rho=1;
    u=1;
    v=.5;
    p=1;
    E=.5*rho*(u^2+v^2)+p/(gamma-1);
    uInit(:,:,1)=rho*one;
    uInit(:,:,2)=u*one;
    uInit(:,:,3)=v*one;
    uInit(:,:,4)=E*one;
elseif Case==2
    %Constant x, 0 y
    rho=1;
    u=2;
    v=0;
    p=.7;
    E=.5*rho*(u^2+v^2)+p/(gamma-1);
    uInit(:,:,1)=rho*one;
    uInit(:,:,2)=rho*u*one;
    uInit(:,:,3)=rho*v*one;
    uInit(:,:,4)=E*one;
    
    
elseif Case==3
    %Constant y, 0 x
    rho=1;
    u=0;
    v=2;
    p=.7;
    E=.5*rho*(u^2+v^2)+p/(gamma-1);
    uInit(:,:,1)=rho*one;
    uInit(:,:,2)=rho*u*one;
    uInit(:,:,3)=rho*v*one;
    uInit(:,:,4)=E*one;
    
elseif Case==4
    %Constant x,y
    rho=1;
    u=1;
    v=2;
    p=.7;
    E=.5*rho*(u^2+v^2)+p/(gamma-1);
    
    uInit(:,:,1)=rho*one;
    uInit(:,:,2)=rho*u*one;
    uInit(:,:,3)=rho*v*one;
    uInit(:,:,4)=E*one;
    
elseif Case==5
    %Density advection 1D -x
    rho=tanh(x)+1;
    u=-1.2;
    v=0;
    p=.7;
    E=.5*rho*(u^2+v^2)+p/(gamma-1)*ones(1,M);
    for k=1:N
        uInit(:,k,1)=rho;
        uInit(:,k,2)=rho*u;
        uInit(:,k,3)=rho*v;
        uInit(:,k,4)=E;
    end
elseif Case==6
    %Density advection 1D -y
    rho=tanh(y)+1;
    u=0;
    v=1.1;
    p=.7;
    E=.5*rho*(u^2+v^2)+p/(gamma-1);
    for j=1:M
        uInit(j,:,1)=rho;
        uInit(j,:,2)=rho*u;
        uInit(j,:,3)=rho*v;
        uInit(j,:,4)=E;
    end
    
elseif Case==7
    %Density advection 2D
    u=1.2;
    v=-1.4;
    p=.7;
    rho=zeros(M,N);
    E=zeros(M,N);
    for j=1:M
        for k=1:N
            rho(j,k)=tanh(x(j))*tanh(y(k))+1;
            E(j,k)=.5*rho(j,k)*(u^2+v^2)+p/(gamma-1);
            uInit(j,k,1)=rho(j,k);
            uInit(j,k,2)=rho(j,k)*u;
            uInit(j,k,3)=rho(j,k)*v;
            uInit(j,k,4)=E(j,k);
        end
    end
    
elseif Case==8
    %x Riemann
    rhoL=1;
    rhoR=.125;
    uL=0;
    uR=0;
    vL=0;
    vR=1;
    pL=1;
    pR=.1;
    EL=.5*rhoL*(uL^2+vL^2)+pL/(gamma-1);
    ER=.5*rhoR*(uR^2+vR^2)+pR/(gamma-1);
    for i=1:(M-1)/2
        uInit(i,:,1)=rhoL;
        uInit(i,:,2)=rhoL*uL;
        uInit(i,:,3)=rhoL*vL;
        uInit(i,:,4)=EL;
    end
    for i=(M-1)/2+1:M
        uInit(i,:,1)=rhoR;
        uInit(i,:,2)=rhoR*uR;
        uInit(i,:,3)=rhoR*vR;
        uInit(i,:,4)=ER;
    end
    
elseif Case==9
    %y Riemann
    rhoT=1;
    rhoB=1;
    uT=2;
    uB=-2;
    vT=3;
    vB=-3;
    pT=.4;
    pB=.4;
    ET=.5*rhoT*(uT^2+vT^2)+pT/(gamma-1);
    EB=.5*rhoB*(uB^2+vB^2)+pB/(gamma-1);
    for i=1:(N)/2
        uInit(:,i,1)=rhoB;
        uInit(:,i,2)=rhoB*uB;
        uInit(:,i,3)=rhoB*vB;
        uInit(:,i,4)=EB;
    end
    for i=(N)/2+1:N
        uInit(:,i,1)=rhoT;
        uInit(:,i,2)=rhoT*uT;
        uInit(:,i,3)=rhoT*vT;
        uInit(:,i,4)=ET;
    end
    
elseif Case==10
    %Continuous x, constant y
    rho=tanh(x)+1;
    u=sin(x);
    v=1;
    p=.7;
    E=zeros(1,M);
    for j=1:M
        E(j)=.5*rho(j)*(u(j)^2+v^2)+p/(gamma-1);
    end
    for k=1:N
        uInit(:,k,1)=rho;
        uInit(:,k,2)=rho.*u;
        uInit(:,k,3)=rho*v;
        uInit(:,k,4)=E;
    end
    
elseif Case==11
    %Continuous y, constant x
    rho=tanh(y)+1;
    v=sin(y);
    u=1;
    p=.7;
    E=zeros(1,M);
    for k=1:N
        E(k)=.5*rho(k)*(u^2+v(k)^2)+p/(gamma-1);
    end
    for j=1:M
        uInit(j,:,1)=rho;
        uInit(j,:,2)=rho.*u;
        uInit(j,:,3)=rho.*v;
        uInit(j,:,4)=E;
    end
elseif Case==12
    %Continuous both
    u=.3;
    v=.2;
    
    for i=1:M
        for j=1:N
            rho=2*tanh(x(i))*tanh(y(j))+2;
            uInit(i,j,1)=rho;
            uInit(i,j,2)=rho*u;
            uInit(i,j,3)=rho*v;
            uInit(i,j,4)=.5*rho*(u^2+v^2)+...
                (tanh(x(i))*tanh(y(j))+1)/(gamma-1);
        end
    end
    
elseif Case==13
    %'Riemann' Data in x&y
    %x Riemann
    rhoL=1;
    rhoR=.125;
    uL=0;
    uR=0;
    vL=0;
    vR=1;
    pL=1;
    pR=.1;
    EL=.5*rhoL*(uL^2+vL^2)+pL/(gamma-1);
    ER=.5*rhoR*(uR^2+vR^2)+pR/(gamma-1);
    %y Riemann
    rhoT=1;
    rhoB=1;
    uT=2;
    uB=-2;
    vT=3;
    vB=-3;
    pT=.4;
    pB=.4;
    ET=.5*rhoT*(uT^2+vT^2)+pT/(gamma-1);
    EB=.5*rhoB*(uB^2+vB^2)+pB/(gamma-1);
    for j=1:M/2
        for k=1:N/2
            uInit(j,k,1)=rhoL;
            uInit(j,k,2)=uL*rhoL;
            uInit(j,k,3)=vL*rhoL;
            uInit(j,k,4)=EL;
        end
        for k=N/2+1:N
            uInit(j,k,1)=rhoR;
            uInit(j,k,2)=uR*rhoR;
            uInit(j,k,3)=vR*rhoR;
            uInit(j,k,4)=ER;
        end
    end
    for j=M/2+1:M
        for k=1:N/2
            uInit(j,k,1)=rhoB;
            uInit(j,k,2)=uB*rhoB;
            uInit(j,k,3)=vB*rhoB;
            uInit(j,k,4)=EB;
        end
        for k=N/2+1:N
            uInit(j,k,1)=rhoT;
            uInit(j,k,2)=uT*rhoT;
            uInit(j,k,3)=vT*rhoT;
            uInit(j,k,4)=ET;
        end
    end
end