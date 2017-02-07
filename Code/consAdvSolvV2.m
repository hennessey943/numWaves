function [dx,err,err1,data,x]=consAdvSolvV2(a,N,xa,xb,cfl,tf,order,iPlot)
%Determines Numerical Solution to the Conservation Law u_t+f(u)_x=0
err=0;
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

Ntot=N+1+2*(ng);

%% First and last interior grid points
ia=ng+1;
ib=Ntot-ng;

x = linspace(xa-ng*dx,xb+(ng)*dx,Ntot);

%% Set initial conditions
t=0;
u=setICs(x,ia,ib);
%% Main loop
for n=1:Nt
    %% Set BCs on current time step
    u=setBCs(u,ia,ib,order);
    %% Compute RK-4 constants
    [k1]=computeK(u,a,ia,ib,dx,order);
    %ua=u+dt*.5*k1;
    [ua]=setBCs(u+dt*.5*k1,ia,ib,order);
    [k2]=computeK(ua,a,ia,ib,dx,order);
    %ub=u+dt*.5*k2;
    [ub]=setBCs(u+dt*.5*k2,ia,ib,order);
    [k3]=computeK(ub,a,ia,ib,dx,order);
    %uc=u+dt*k3;
    [uc]=setBCs(u+dt*k3,ia,ib,order);
    [k4]=computeK(uc,a,ia,ib,dx,order);
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
end
err1=0;
for i=1:N
    err1=dx*abs(u(i)-ue(i))+err1;
end

data=u;
%% Initial Condition Functions
    function z=f(x)
%         if x<-pi*.25
%             z=0;
%         else
%             z=(cos(2*x)).^5;
%         end
        if x>-.75&&x<=.25
            z=1;
        else
            z=0;
        end
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
    function u=setBCs(u,ia,ib,order)
        %u(ia)=0;
        if order==1
            u(ia-1)=u(ia);
            %u(ib+1)=-u(ib-1);
        elseif order==2
            u(ia-2)=u(ia);
            u(ia-1)=u(ia);
            u(ib+1)=3*u(ib)-3*u(ib-1)+u(ib-2);
            u(ib+2)=3*u(ib+1)-3*u(ib)+u(ib-1);
        elseif order==3
            u(ia-3)=u(ia);
            u(ia-2)=u(ia);
            u(ia-1)=u(ia);
            u(ib+1)=4*u(ib)-6*u(ib-1)+4*u(ib-2)-u(ib-3);
            u(ib+2)=4*u(ib+1)-6*u(ib)+4*u(ib-1)-u(ib-2);
            u(ib+3)=4*u(ib+2)-6*u(ib+1)+4*u(ib)-u(ib-1);
        elseif order==4
            u(ia-4)=u(ia);
            u(ia-3)=u(ia);
            u(ia-2)=u(ia);
            u(ia-1)=u(ia);
            u(ib+1)=4*u(ib)-6*u(ib-1)+4*u(ib-2)-u(ib-3);
            u(ib+2)=4*u(ib+1)-6*u(ib)+4*u(ib-1)-u(ib-2);
            u(ib+3)=4*u(ib+2)-6*u(ib+1)+4*u(ib)-u(ib-1);
            u(ib+4)=4*u(ib+3)-6*u(ib+2)+4*u(ib+1)-u(ib);
        end
    end
%% Compute RK-4 Coefficient Function
    function K=computeK(u,a,ia,ib,dx,order)
        K=zeros(size(u));
        for j=ia:ib
            if order==1
                flp=u(j);
                flm=u(j-1);
                flpxx=0;
                flmxx=0;
            elseif order==2
                uxp=(u(j+1)-u(j-1))/(2*dx);
                uxm=(u(j)-u(j-2))/(2*dx);
                flp=u(j)+dx*uxp/2;
                flm=u(j-1)+dx*uxm/2;
                flpxx=0;
                flmxx=0;
            elseif order==3
                uxp=(u(j+1)-u(j-1))/(2*dx);
                uxm=(u(j)-u(j-2))/(2*dx);
                uxxp=(u(j+1)-2*u(j)+u(j-1))/(dx^2);
                uxxm=(u(j)-2*u(j-1)+u(j-2))/(dx^2);
                flp=u(j)+dx*uxp/2+dx^2*uxxp/8;
                flm=u(j-1)+dx*uxm/2+dx^2*uxxm/8;
                flpxx=uxxp;
                flmxx=uxxm;
            elseif order==4
                uxxp=(u(j+1)-2*u(j)+u(j-1))/(dx^2);
                uxxm=(u(j)-2*u(j-1)+u(j-2))/(dx^2);
                uxp=(-u(j+2)+8*u(j+1)-8*u(j-1)+u(j-2))/(12*dx);
                uxm=(-u(j+1)+8*u(j)-8*u(j-2)+u(j-3))/(12*dx);
                flp=u(j)+dx*uxp/2+dx^2*uxxp/8;
                flm=u(j-1)+dx*uxm/2+dx^2*uxxm/8;
                flpxx=uxxp;
                flmxx=uxxm;
            end
            K(j)=-a*(flp-dx^2*flpxx/24-(flm-dx^2*flmxx/24))/dx;
        end
    end
%% Plot Function
    function plotStuff(x,u,ue,ind)
        subplot(2,1,1)
        plot(x(ind),u(ind),'r-',x(ind),ue(ind),'k-')
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