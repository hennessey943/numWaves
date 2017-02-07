function [data,dx,x,y]=waveD2RichTest(M,N,xa,xb,ya,yb,c,cfl,tf,iPlot,iOrder)

%% Set Global variables
global dx dy dt ng ja jb ka kb x y Mtot Ntot

%% Initialize grid
if iOrder~=2&&iOrder~=4
    fprintf('%f is not a valid Order for this Solver',iOrder)
    err=NaN;
    dx=NaN;
    dy=NaN;
    return
else

dx=(xb-xa)/M;
dy=(yb-ya)/N;
dt=cfl*(dx*dy)/(c*sqrt(dx^2+dy^2));
Nt=ceil(tf/dt);
dt=tf/Nt;
%assign ghost points
if iOrder==2
    ng=1;
elseif iOrder==4
    ng=2;
end

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
    for k=ka:kb
        for j=ja:jb
            if iOrder==2
                %uxx=-(sin(x(j)-sqrt(2)*c*t)*cos(y(k))+sin(x(j)+sqrt(2)*c*t)*cos(y(k)));
                %uyy=-uxx;
                uxx=(un(j+1,k)-2*un(j,k)+un(j-1,k))/dx^2;
                uyy=(un(j,k+1)-2*un(j,k)+un(j,k-1))/dy^2;
                utt=uxx+uyy;
                u(j,k)=2*un(j,k)-unm1(j,k)+c^2*dt^2*utt;
                %u(j,k)=(sin(x(j)-sqrt(2)*c*t)*cos(y(k))+sin(x(j)+sqrt(2)*c*t)*cos(y(k)));
            elseif iOrder==4
                uxx=(-un(j-2,k)+16*un(j-1,k)-30*un(j,k)+16*un(j+1,k)-un(j+2,k))/(12*dx^2);
                uyy=(-un(j,k-2)+16*un(j,k-1)-30*un(j,k)+16*un(j,k+1)-un(j,k+2))/(12*dx^2);
                utt=c^2*(uxx+uyy);
                uxxxx=(un(j-2,k)-4*un(j-1,k)+6*un(j,k)-4*un(j+1,k)+un(j+2,k))/dx^4;
                uxxyy=(un(j-1,k-1)-2*un(j-1,k)+un(j-1,k+1)-2*(un(j,k-1)-2*un(j,k)+un(j,k+1))+...
                        un(j+1,k-1)-2*un(j+1,k)+un(j+1,k+1))/(dx^2*dy^2);
                uyyyy=(un(j,k-2)-4*un(j,k-1)+6*un(j,k)-4*un(j,k+1)+un(j,k+2))/dx^4;
                utttt=c^4*(uxxxx+2*uxxyy+uyyyy);
                u(j,k)=2*un(j,k)-unm1(j,k)+dt^2*utt+dt^4*utttt/12;
            end
        end
    end

    %% Update solution histories
    unm1=un;
    un=u;
    t=t+dt;
    %% Plot stuff
    if iPlot==1
        plotStuff(u)
    end
    
end
end
data=u;
t


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

    function z=g(~,~)
        z=0;
    end
    function z=gxx(~,~)
        z=0;
    end
    function z=gyy(~,~)
        z=0;
    end
    function z=f(m,n)
%         if pi/4<=m&&m<=3*pi/4&&pi/4<=n&&n<=3*pi/4
%             z=(cos(2*m)*cos(2*n))^5;
%         else
%             z=0;
%         end
%         if m>pi/4 && m<=pi/2 && n>pi/4&& n<=pi/2
%             z=(4*m-pi)*(4*n-pi);
%         elseif m>pi/2&&m<=3*pi/4 && n>pi/2&&n<=3*pi/4
%             z=(-4*m+3*pi)*(-4*n+3*pi);
%         else
%             z=0;
%         end
        if m>pi/3&&m<=2*pi/3 && n>pi/3&& n<=2*pi/3
            z=1;
        else
            z=0;
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
        
        %% Initialize at t=dt
        %% we use a Taylor expansion in time
        un=zeros(Mtot,Ntot);
        
        for a=ja:jb
            for b=ka:kb
                if iOrder==2
                    fxx=(f(x(a+1),y(b))-2*f(x(a),y(b))+...
                        f(x(a-1),y(b)))/dx^2;
                    fyy=(f(x(a),y(b+1))-2*f(x(a),y(b))+...
                        f(x(a),y(b-1)))/dy^2;
                    un(a,b)=f(x(a),y(b))+dt*g(x(a),y(b))+...
                            c^2*dt^2/2*(fxx+fyy);
                elseif iOrder==4
                    fxx=(-f(x(a-2),y(b))+16*f(x(a-1),y(b))-...
                        30*f(x(a),y(b))+16*f(x(a+1),y(b))-...
                        f(x(a+2),y(b)))/(12*dx^2);
                    fyy=(-f(x(a),y(b-2))+16*f(x(a),y(b-1))-...
                        30*f(x(a),y(b))+16*f(x(a),y(b+1))-...
                        f(x(a),y(b+2)))/(12*dy^2);
                    fxxxx=(f(x(a-2),y(b))-4*f(x(a-1),y(b))+...
                        6*f(x(a),y(b))-4*f(x(a+1),y(b))+...
                        f(x(a+2),y(b)))/dx^4;
                    fxxyy=(f(x(a-1),y(b-1))-2*f(x(a-1),y(b))+...
                        f(x(a-1),y(b+1))-2*(f(x(a),y(b-1))-2*f(x(a),y(b))+...
                        f(x(a),y(b+1)))+f(x(a+1),y(b-1))-2*f(x(a+1),y(b))+...
                        f(x(a+1),y(b+1)))/(dx^2*dy^2);
                    fyyyy=(f(x(a),y(b-2))-4*f(x(a),y(b-1))+...
                        6*f(x(a),y(b))-4*f(x(a),y(b+1))+...
                        f(x(a),y(b+2)))/dy^4;
                    un(a,b)=f(x(a),y(b))+dt*g(x(a),y(b))+...
                            c^2*dt^2/2*(fxx+fyy)+...
                            c^2*dt^3/6*(gxx(x(a),y(b))+gyy(x(a),y(b)))+...
                            c^4*dt^4/24*(fxxxx+2*fxxyy+fyyyy);
                end      
            end
        end
        un=setBCs(un,dt);
        
        return
    end
%% Boundary condition function
    function un=setBCs(un,t)
        for a=ka:kb
            if iOrder==2
                uyyl=(un(ja,a+1)-2*un(ja,a)+un(ja,a-1))/dy^2;
                %uyyl=-(sin(x(ja)-sqrt(2)*c*t)*cos(y(a))+sin(x(ja)+sqrt(2)*c*t)*cos(y(a)));
                un(ja-1,a)=2*un(ja,a)-un(ja+1,a)+dx^2/c^2*gltt(y(a),t)-dx^2*uyyl;
                
                uyyr=(un(jb,a+1)-2*un(jb,a)+un(jb,a-1))/dy^2;
                %uyyr=-(sin(x(jb)-sqrt(2)*c*t)*cos(y(a))+sin(x(jb)+sqrt(2)*c*t)*cos(y(a)));
                un(jb+1,a)=2*un(jb,a)-un(jb-1,a)+dx^2/c^2*grtt(y(a),t)-dx^2*uyyr;
            elseif iOrder==4
                un(ja-2,a)=-un(ja+2,a);
                un(ja-1,a)=-un(ja+1,a);
                un(jb+1,a)=-un(jb-1,a);
                un(jb+2,a)=-un(jb-2,a);
            end
                
        end
        for b=ja:jb
            if iOrder==2
                un(b,ka-1)=un(b,ka+1)-2*dy*bl(x(b),t);
                un(b,kb+1)=un(b,kb-1)+2*dy*br(x(b),t);
            elseif iOrder==4
                un(b,ka-2)=un(b,ka+2);
                un(b,ka-1)=un(b,ka+1);
                un(b,kb+1)=un(b,kb-1);
                un(b,kb+2)=un(b,kb-2);
            end
        end
        %Corner Points
        un(ja-1,ka-1)=-un(ja+1,ka+1);
        un(jb+1,ka-1)=-un(jb-1,ka+1);
        un(ja-1,kb+1)=-un(ja+1,kb-1);
        un(jb+1,kb+1)=-un(jb-1,kb-1);
    end

%% Plot function
    function plotStuff(u)
        xind=ja:jb;
        yind=ka:kb;
        surf(x(xind),y(yind),u(xind,yind))
        xlabel('x');
        ylabel('y');
        zlabel('u');
        drawnow;
        pause
        return
    end
end




