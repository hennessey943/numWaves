function [dx,err]=consAdvSolv(a,N,xa,xb,cfl,tf,order,iPlot)
%Determines Numerical Solution to the Conservation Law u_t+f(u)_x=0

%% Initialize Variables
%Number of ghost points
if order==1
    ng=1;
elseif order==2
    ng=2;
elseif order==3
    ng=3;
elseif order==4
    ng=4;
end

%% Grid size and other grid quantities
dx=(xb-xa)/N;
dt=cfl*dx/a;
Nt=ceil(tf/dt);
dt=tf/Nt;

Ntot=N+1+ng;

%% First and last interior grid points
ia=ng+1;
ib=Ntot;

x = linspace(xa-ng*dx,xb,Ntot);

%% Set initial conditions
t=0;
u=setICs(x,ia,ib);
%% Main loop
for n=1:Nt-1
    %% Set BCs on current time step
    u=setBCs(u,ia,order);
    %% Compute RK-4 constants
    [k1]=computeK(u,a,ia,ib,dx,dt,order);
    %ua=u+dt*.5*k1;
    [ua]=setBCs(u+dt*.5*k1,ia,order);
    [k2]=computeK(ua,a,ia,ib,dx,dt,order);
    %ub=u+dt*.5*k2;
    [ub]=setBCs(u+dt*.5*k2,ia,order);
    [k3]=computeK(ub,a,ia,ib,dx,dt,order);
    %uc=u+dt*k3;
    [uc]=setBCs(u+dt*k3,ia,order);
    [k4]=computeK(uc,a,ia,ib,dx,dt,order);
    %Update solution
    for i=ia:ib
        u(i)=u(i)+dt/6*(k1(i)+2*k2(i)+2*k3(i)+k4(i));
    end
    
    %Update time
    t=t+dt;
    
    %% Check error
    [ue]=uex(x-a*t);
    ind=ia:ib;
    err=max(abs(u(ind)-ue(ind)));
    %% Plot Solution
    if iPlot==1
        plotStuff(x,u,ue,ind);
    end
%     if t>=m*dtplot
%         time(m)=t;
%         saved_data(m,1:N)=u;
%         m=m+1;
%     end
end

%% Initial Condition Functions
    function z=f(x)
        if x<-pi*.25
            z=0;
        else
            z=(cos(2*x)).^5;
        end
%         if x>-.75&&x<=.25
%             z=1;
%         else
%             z=0;
%         end
    end

%Sets condition at t=0
    function u=setICs(x,ia,ib)
        u=zeros(size(x));
        for j=ia:ib
            u(j)=f(x(j));
        end
    end

%% Exact Solution
    function z=uex(x)
        z=zeros(1,numel(x));
        for j=1:numel(x)
            z(j)=f(x(j));
        end
    end

%% Boundary Condition Function
    function u=setBCs(u,ia,order)
        %u(ia)=0;
        if order==1
            u(ia-1)=u(ia);
        elseif order==2
            u(ia-2)=u(ia);
            u(ia-1)=u(ia);
        elseif order==3
            u(ia-3)=u(ia);
            u(ia-2)=u(ia);
            u(ia-1)=u(ia);
        elseif order==4
            u(ia-4)=u(ia);
            u(ia-3)=u(ia);
            u(ia-2)=u(ia);
            u(ia-1)=u(ia);
        end
    end
%% Compute RK-4 Coefficient Function
    function K=computeK(u,a,ia,ib,dx,dt,order)
        K=zeros(size(u));
        for j=ia:ib
            if order==1
                flp=u(j);
                flm=u(j-1);
                flpxx=0;
                flmxx=0;
            elseif order==2
                uxp=(u(j)-u(j-1))/dx;
                uxm=(u(j-1)-u(j-2))/dx;
                flp=u(j)+dx*(1-a*dt/dx)*uxp/2;
                flm=u(j-1)+dx*(1-a*dt/dx)*uxm/2;
                flpxx=0;
                flmxx=0;
            elseif order==3
                %uxp=(u(j)-u(j-1))/dx;
                uxp=(3*u(j)-4*u(j-1)+u(j-2))/(2*dx);
                uxxp=(u(j)-2*u(j-1)+u(j-2))/(dx^2);
                %uxm=(u(j-1)-u(j-2))/dx;
                uxm=(3*u(j-1)-4*u(j-2)+u(j-3))/(2*dx);
                uxxm=(u(j-1)-2*u(j-2)+u(j-3))/(dx^2);
                %flp=u(j)+dx*(1-a*dt/dx)*uxp/2+dx^2*(1-2*a*dt/dx+a^2*dt^2/dx^2)*uxxp/8;
                flp=u(j)+dx*uxp/2+dx^2*uxxp/8;
                %flm=u(j-1)+dx*(1-a*dt/dx)*uxm/2+dx^2*(1-2*a*dt/dx+a^2*dt^2/dx^2)*uxxm/8;
                flm=u(j-1)+dx*uxm/2+dx^2*uxxm/8;
                flpxx=uxxp;
                flmxx=uxxm;
            elseif order==4
                uxp=(11*u(j)-18*u(j-1)+9*u(j-2)-2*u(j-3))/(6*dx);
                uxm=(11*u(j-1)-18*u(j-2)+9*u(j-3)-2*u(j-4))/(6*dx);
                uxxp=(2*u(j)-5*u(j-1)+4*u(j-2)-u(j-3))/(dx^2);
                uxxm=(2*u(j-1)-5*u(j-2)+4*u(j-3)-u(j-4))/(dx^2);
                uxxxp=(u(j)-3*u(j-1)+3*u(j-2)-u(j-3))/(dx^3);
                uxxxm=(u(j-1)-3*u(j-2)+3*u(j-3)-u(j-4))/(dx^3);
                flp=u(j)+dx*uxp/2+dx^2*uxxp/8+dx^3*uxxxp/48;
                flm=u(j-1)+dx*uxm/2+dx^2*uxxm/8+dx^3*uxxxm/48;
                flpxx=uxxp+dx/2*uxxxp;
                flmxx=uxxm+dx/2*uxxxm;
            end
            %K(j)=-a*(flp-flm)/dx;
            K(j)=-a*(flp-dx^2*flpxx/24-(flm-dx^2*flmxx/24))/dx;
            %K(j)=-a*(11*u(j)-18*u(j-1)+9*u(j-2)-2*u(j-3))/(6*dx);
        end
    end
%% Plot Function
    function plotStuff(x,u,ue,ind)
        subplot(2,1,1)
        plot(x(ind),u(ind),'rx',x(ind),ue(ind),'k-')
        xlabel('x')
        ylabel('u')
        legend('num u','exact u')
        subplot(2,1,2)
        plot(x(ind),u(ind)-ue(ind),'rx')
        xlabel('x')
        ylabel('error')
        drawnow
        pause
    end    
end