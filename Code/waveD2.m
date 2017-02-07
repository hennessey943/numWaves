function [err,dx,dy]=waveD2(M,N,xa,xb,ya,yb,c,tf,iPlot)

%% Set Global variables
global dx dy dt ng ja jb ka kb x y Mtot Ntot

%% Initialize grid
dx=(xb-xa)/M;
dy=(yb-ya)/N;
dt=.9*(dx*dy)/(c*sqrt(dx^2+dy^2));
Nt=ceil(tf/dt);
dt=tf/Nt;
%assign ghost points
ng=1;

% data set for grid
Mtot=M+1+2*ng;
ja=ng+1;
jb=Mtot-ng;

Ntot=N+1+2*ng;
ka=ng+1;
kb=Ntot-ng;

x=linspace(xa-ng*dx,xb+ng*dx,Mtot);
y=linspace(ya-ng*dx,yb+ng*dx,Ntot);

%% set Initial conditions
[unm1,un]=setICs;

%% allocate space for u

u=zeros(size(un));
if iPlot==1
    figure
end



%% Loop through time
t=dt;
for n=1:Nt-1
    
    %% Set boundary conditions on current step
    un=setBCs(un,t);
    %% Compute update over domain interior
    for j=ja:jb
        for k=ka:kb
            %uxx=-(sin(x(j)-sqrt(2)*c*t)*cos(y(k))+sin(x(j)+sqrt(2)*c*t)*cos(y(k)));
            %uyy=-uxx;
            uxx=(un(j+1,k)-2*un(j,k)+un(j-1,k))/dx^2;
            uyy=(un(j,k+1)-2*un(j,k)+un(j,k-1))/dy^2;
            utt=uxx+uyy;
            u(j,k)=2*un(j,k)-unm1(j,k)+c^2*dt^2*utt;
            %u(j,k)=(sin(x(j)-sqrt(2)*c*t)*cos(y(k))+sin(x(j)+sqrt(2)*c*t)*cos(y(k)));
        end
    end

    %% Update solution histories
    unm1=un;
    un=u;
    t=t+dt;
    %% Plot stuff
    if iPlot==1
        plotStuff(u,t)
    end
    
    
    %% Calculate exact solution and check error
    ue=uex(t);
    xind=ja:jb;
    yind=ka:kb;
    
    err=max(max(abs(u(xind,yind)-ue(xind,yind))));
    %err=norm(u(xind,yind)-ue(xind,yind),1);
    

end




%% Define boundary functions and derivatives
    function z=gl(~,~)
        z=0;
    end
    function z=gr(~,~)
        z=0;
    end
    function z=gltt(~,~)
        z=0;
    end
    function z=grtt(~,~)
        z=0;
    end
    function z=bl(~,~)
        z=0;
    end
    function z=br(~,~)
        z=0;
    end
    function z=f(m,n)
        z=2*sin(m)*cos(n);
    end
    function z=g(a,b)
        z=-sqrt(2)*c*cos(a)*cos(b)+sqrt(2)*c*cos(a)*cos(b);
    end
    function z=fxx(a,b)
        z=-2*sin(a)*cos(b);
    end
    function z=fyy(a,b)
        z=-2*sin(a)*cos(b);
    end

    function z=uex(t)
        z=zeros(M,N);
        for a=ja:jb
            for b=ka:kb
                z(a,b)=sin(x(a)-sqrt(2)*c*t)*cos(y(b))+sin(x(a)+sqrt(2)*c*t)*cos(y(b));
            end
        end
    end
%% Initial condition function
    function [unm1,un]=setICs
        
        %Set intitial condition at t=0:
        unm1=zeros(Mtot,Ntot);
        for a=ja:jb
            for b=ka:kb
                unm1(a,b)=f(x(a),y(b));
            end
        end
        
        unm1=setBCs(unm1,0);
        z=uex(0);
        
        %% Initialize at t=dt
        %% we use a Taylor expansion in time
        un=zeros(Mtot,Ntot);
        for a=ja:jb
            for b=ka:kb
                un(a,b)=f(x(a),y(b))+dt*g(x(a),y(b))+...
                    c^2*dt^2/2*(fxx(x(a),y(b))+fyy(x(a),y(b)));
                %un(a,b)=sin(x(a)-sqrt(2)*c*dt)*cos(y(b))+sin(x(a)+sqrt(2)*c*dt)*cos(y(b));
            end
        end
        un=setBCs(un,dt);
        
        return
    end
%% Boundary condition function
    function un=setBCs(un,t)
        for a=ka:kb
            uyyl=(un(ja,a+1)-2*un(ja,a)+un(ja,a-1))/dy^2;
            %uyyl=-(sin(x(ja)-sqrt(2)*c*t)*cos(y(a))+sin(x(ja)+sqrt(2)*c*t)*cos(y(a)));
            un(ja-1,a)=2*un(ja,a)-un(ja+1,a)+dx^2/c^2*gltt(y(a),t)-dx^2*uyyl;
            
            uyyr=(un(jb,a+1)-2*un(jb,a)+un(jb,a-1))/dy^2;
            %uyyr=-(sin(x(jb)-sqrt(2)*c*t)*cos(y(a))+sin(x(jb)+sqrt(2)*c*t)*cos(y(a)));
            un(jb+1,a)=2*un(jb,a)-un(jb-1,a)+dx^2/c^2*grtt(y(a),t)-dx^2*uyyr;
        end
        for b=ja:jb
            un(b,ka-1)=un(b,ka+1)-2*dy*bl(x(b),t);
            un(b,kb+1)=un(b,kb-1)+2*dy*br(x(b),t);
        end
    end

%% Plot function
    function plotStuff(u,t)
        ue=uex(t);
        xind=ja:jb;
        yind=ka:kb;
        subplot(2,1,1)
        surf(x(xind),y(yind),u(xind,yind))
        hold on
        surf(x(xind),y(yind),ue(xind,yind))
        xlabel('x');
        ylabel('y');
        zlabel('u');
        hold off
        
        subplot(2,1,2)
        surf(u(xind,yind)-ue(xind,yind))
        xlabel('x')
        ylabel('y')
        zlabel('error');
        drawnow;
        pause
        return
    end
end




