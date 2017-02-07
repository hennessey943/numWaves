function [dx,v_err,sig_err]=colFOSsolver(N,xa,xb,cfl,tf,iPlot,Case)
% Solves the first order system for the wave equation on a collocated grid

%% Number of Ghost points
ng=1;
dx=(xb-xa)/N;
Ntot=N+1+2*ng;

%% first and last interior grid points
ia=ng+1;
ib=Ntot-ng;
x=linspace(xa-ng*dx,xb+ng*dx,Ntot);
x_ghost=linspace(xa-2*dx,xb+2*dx,Ntot+2);

%% grid size and other quantities

dt=cfl*dx;
Nt=ceil(tf/dt);
dt=tf/Nt;

[un]=setICs(x,x_ghost,dx,ia,ib,Case);

%% allocate space for u
u=zeros(size(un));

%% Loop through time
t=0;
for n=1:Nt-1
    %% Set BCs on u
    un=setBCs(un,ia,ib);
    %% Compute Runge-Kutta coefficient matrices
    [k1]=computeK(un,ia,ib,dx);
    [ua]=setBCs(un+dt*.5*k1,ia,ib);
    [k2]=computeK(ua,ia,ib,dx);
    [ub]=setBCs(un+dt*.5*k2,ia,ib);
    [k3]=computeK(ub,ia,ib,dx);
    [uc]=setBCs(un+dt*k3,ia,ib);
    [k4]=computeK(uc,ia,ib,dx);
    %% Compute update over domain interior
    for i=ia:ib
        u(:,i)=un(:,i)+dt/6*(k1(:,i)+2*k2(:,i)+2*k3(:,i)+k4(:,i));
    end

    %% update solution histories
    un=u;
    t=t+dt;
    
    %% Plot solution
    if(iPlot==1)
        plotStuff(x,u,ia,ib,t,Case);
    end
    
    %% Check error
    if Case==1
    [v,sigma]=ue(x,t);
    ind=ia:ib;
    v_err=max(abs(u(1,ind)-v(ind)));
    sig_err=max(abs(u(2,ind)-sigma(ind)));
    end
end

%% Compute k for RK4
    function K=computeK(un,ia,ib,dx)
        K=zeros(size(un));
        A=[0,1;1,0];
        for j=ia:ib
            K(:,j)=.5*A*(un(:,j+1)-un(:,j-1))/dx;
        end
        return
    end
%% Boundary condition setter
function un=setBCs(un,ia,ib)
    un(1,ia)=0;
    un(1,ia-1)=2*un(1,ia)-un(1,ia+1);
    un(2,ia-1)=un(2,ia+1);
    un(1,ib+1)=un(1,ib-1);
    un(2,ib)=0;
    un(2,ib+1)=2*un(2,ib)-un(2,ib-1);
    return
end

%% Exact Solution
    function [v,sigma]=ue(x,t)
        v=-pi*.5*sin(pi*t*.5)*sin(pi*x*.5);
        sigma=pi*.5*cos(pi*t*.5)*cos(pi*x*.5);
        
        return
    end

%% Solution at t=0
    function z=f(x,Case)
        if Case==1
            z=sin(pi*x*.5);
        elseif Case==2
            z=exp(-100*(x-.5)^2);
        elseif Case==3
            z=zeros(size(x));
            for m=1:numel(x)
                if x(m)>.25&&x(m)<=.5
                    z(m)=4*x(m)-1;
                elseif x(m)>.5&&x(m)<.75
                    z(m)=-4*x(m)+3;
                else
                    z(m)=0;
                end
            end
        elseif Case==4
            z=zeros(size(x));
            for m=1:numel(x)
                if x(m)>.25&&x(m)<.75
                    z=1;
                else
                    z=0;
                end
            end
        end
        return
    end

    function z=fx(m,x,dx,Case)
        if Case==1
            z=pi*.5*cos(pi*x(m)*.5);
        elseif Case==2
            z=-200*(x(m)-.5)*exp(-100*(x(m)-.5)^2);
        elseif Case==3 || Case==4
            % 4th order Dx of f(x)
            z=(-f(x(m+2),Case)+8*f(x(m+1),Case)-...
                8*f(x(m-1),Case)+f(x(m-2),Case))/(12*dx);
%         elseif Case==4
%             z=(f(x(m+1))-f(x(m-1)))/(2*dx);
        end
        return
    end

    function z=g(~)
        z=0;
        return
    end
%% Plot function
    function plotStuff( x,u,ia,ib,t,Case )
        if Case==1
        [v,sigma] = ue( x,t );
        
        ind = ia:ib;
        subplot( 2,1,1 )
        plot( x(ind),u(1,ind),'rx', x(ind),v(ind),'k-' );
        hold on
        plot(x(ind),u(2,ind),'bx',x(ind),sigma(ind),'g-');
        hold off
        xlabel( 'x' );
        ylabel( 'u' );
        legend('num v','exact v','num sigma','exact sigma')
        
        subplot(2,1,2)
        plot( x(ind),u(1,ind)-v(ind),'rx' );
        hold on
        plot(x(ind),u(2,ind)-sigma(ind),'bx');
        hold off
        xlabel( 'x' );
        ylabel( 'error' );
        drawnow;
        pause
        else
            ind=ia:ib;
            plot(x(ind),u(1,ind));
            hold on
            plot(x(ind),u(2,ind));
            hold off
            title(['Case ' num2str(Case) 'at t=' num2str(t)])
            %axis([0,1,-10,10])
            xlabel('x');
            ylabel('soln');
            legend('v','\sigma')
            drawnow
            pause
        end
        
        return
    end

%% Initial condition setter
    function [un]=setICs(x,x_ghost,dx,ia,ib,Case)
        %% set initial condition at t=0(then BCs on IC)
        un=zeros(2,numel(x));
        for j=ia:ib
            un(1,j)=0;
            un(2,j)=fx(j+1,x_ghost,dx,Case);
        end
        un=setBCs(un,ia,ib);
    end
end