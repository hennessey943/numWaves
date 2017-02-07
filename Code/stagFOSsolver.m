function [dx,v_err,sig_err]=stagFOSsolver(N,xa,xb,cfl,tf,iPlot,Case)
% Solves the first order system for the wave equation on a staggered grid
%treat v=u_t and sigma=u_x separately
%% Number of Ghost points
ng=0;
dx=(xb-xa)/(N-.5);
Ntot=N+2*ng; %For each staggered grid

%% First and last interior grid points
% Staggered Left
ia=ng+1;
ib=Ntot-ng;
x_left=linspace(xa-ng*dx,xb-dx/2+ng*dx,Ntot);
%Staggered Right note ja=ia and jb=ib
ja=ng+1;
jb=Ntot-ng;
x_right=linspace(xa+dx/2-ng*dx,xb+ng*dx,Ntot);

%% Set up ghost grid for calculating fx and fxxx for the discontinuous cases
x_ghost=linspace(xa+dx/2-2*dx,xb+2*dx,N+4);
%% Grid size and other quantities
dt=cfl*dx;
Nt=ceil(tf/dt);
dt=tf/Nt;
[vn,sigman]=setICs(x_left,x_right,x_ghost,dx,dt,ja,jb,Case);
%% Allocate space for solution
v=zeros(1,Ntot);
sig=zeros(1,Ntot);
%% Loop through time
tv=0;
ts=dt/2;
for n=1:Nt-1
    %% Set BCs on sig and v
    [vn,sigman]=setBCs(vn,sigman,ia,jb);
    %% Compute update over domain interior

    for i=ia+1:ib
        v(i)=vn(i)+dt/dx*(sigman(i)-sigman(i-1));
    end
    for i=ja:jb-1
        sig(i)=sigman(i)+dt/dx*(v(i+1)-v(i));
    end
    
    %% Update solution histories
    sigman=sig;
    vn=v;
    tv=tv+dt;
    ts=ts+dt;
    
    %% Plot solution
    if iPlot==1
        plotStuff(x_left,v,x_right,sig,ia,ib,ja,jb,tv,ts,Case);
    end
    
    %% Check error
    if Case==1
    [ve,sige]=ue(x_left,x_right,tv,ts);
    ind=ia:ib;
    jnd=ja:jb;
    v_err=max(abs(v(ind)-ve(ind)));
    sig_err=max(abs(sig(jnd)-sige(jnd)));
    end
end

%% Boundary condition setter
    function [vn,sigman]=setBCs(vn,sigman,ia,jb)
        vn(ia)=0;
        sigman(jb)=0;
        return
    end
%% Initial condition setter
    function [vn,sigman]=setICs(x_left,x_right,x_ghost,dx,dt,ja,jb,Case)
        vn=zeros(1,numel(x_left));
        sigman=zeros(1,numel(x_right));
        for j=ja:jb-1
            sigman(j)=fx(j+2,x_ghost,dx,Case)+dt^2/8*fxxx(j+2,x_ghost,dx,Case);
        end
        [vn,sigman]=setBCs(vn,sigman,1,jb);
        return
    end
        
%% Exact Solution
    function [ve,sige]=ue(x_left,x_right,tv,ts)
        
        ve=-pi*.5*sin(pi*tv*.5)*sin(pi*x_left*.5);
        sige=pi*.5*cos(pi*ts*.5)*cos(pi*x_right*.5);
        return
    end

%% Solution at t=0
    function z=f(x,Case)
        %z=sin(pi*x*.5);
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
%     function z=fx(x)
%         z=pi*.5*cos(pi*x*.5);
%         return
%     end
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
    function z=fxxx(m,x,dx,Case)
        if Case==1
            z=-(pi*.5)^3*cos(pi*x(m)*.5);
        else
            %2nd order Dxxx of f
            z=(f(x(m+2),Case)-2*f(x(m+1),Case)+2*f(x(m-1),Case)-...
                f(x(m-2),Case))/(2*dx^3);
        end
        return
    end
%% Plot function
    function plotStuff(x_left,v,x_right,sig,ia,ib,ja,jb,tv,ts,Case)
        if Case==1
        [ve,sige]=ue(x_left,x_right,tv,ts);
        ind=ia:ib;
        jnd=ja:jb;
        subplot(2,1,1)
        plot(x_left(ind),v(ind),'rx',x_left(ind),ve(ind),'k-');
        hold on
        plot(x_right(jnd),sig(jnd),'bx',x_right(jnd),sige(jnd),'g-');
        hold off
        xlabel('x')
        ylabel('u')
        legend('num v', 'exact v', 'num sig','exact sig')
        
        subplot(2,1,2)
        plot(x_left(ind),v(ind)-ve(ind),'rx')
        hold on
        plot(x_right(jnd),sig(jnd)-sige(jnd),'bx')
        hold off
        xlabel('x')
        ylabel('error')
        drawnow;
        pause
        else
            ind=ia:ib;
            plot(x_left(ind),v(ind));
            hold on
            plot(x_right(ind),sig(ind));
            hold off
            xlabel('x')
            ylabel('u')
            legend('num v','num sigma')
            title(['Case ' num2str(Case) 'at t=' num2str(tv)])
            drawnow
            pause
        end
    end
end
