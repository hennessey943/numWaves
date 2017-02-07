function [dx,err] = nonConstWave( N,xa,xb,tf,iPlot )
%xa=0;
%xb=1;

%% number of ghost points
ng = 1;
dx = (xb-xa)/N;
Ntot = N+2*ng;
%% first and last interior grid points
ia = ng+1;
ib = Ntot-ng;
x = linspace( xa+dx/2-ng*dx, xb-dx/2+ng*dx, Ntot );

%% grid size and other grid quantities

a=max(cos(x));
dt = 0.8*dx/a;
Nt = ceil( tf/dt );
dt = tf/Nt;


[unm1,un] = setICs( x,ia,ib,dx,dt);

%% allocate space for u
u = zeros(size(un));


%% Loop through time
t = dt;
for n = 1:Nt-1
    
    %% set BCs on un
    un = setBCs( un,ia,ib);
    
    %% compute update over domain interior
    for i = ia:ib
        i_plus=(i-1)*dx;
        i_minus=(i-2)*dx;
        uxx = ((c(i_plus))^2*(un(i+1)-un(i))-...
            (c(i_minus))^2*(un(i)-un(i-1)))/(dx^2);%+...
            %h(i_plus,t);
        u(i) = 2.*un(i)-unm1(i)+dt^2*(uxx+h(x(i),t));
    end   
        %% update solution histories
        unm1 = un;
        un   = u;
        t = t+dt;
        
        %% Plot Solution
        if( iPlot == 1 )
            plotStuff( x,u,ia,ib,t );
        end
       
    %% Check error
    u_exact = ue( x,t );
    
    ind = ia:ib;
    err = max(abs(u(ind)-u_exact(ind)));
%    return
end

%% Boundary condition setter
    function un = setBCs( un,ia,ib)
        un(ia-1)=un(ia);
        un(ib+1)=un(ib);
        return
    end

%% Exact Solution
    function z = ue( x,t )
        z = cos(t)*cos(pi*x);
        return
    end

    function z = uex ( x,t )
        z=-pi*cos(t)*sin(pi*x);
        return
    end

    function z= uett(x,t)
        z=-cos(t)*cos(pi*x);
        return
    end
    function z=uexx(x,t)
        z=-pi^2*cos(t)*cos(pi*x);
        return
    end
%% Source Term
    function z=h(x,t)
        z=uett(x,t)-2*c(x)*cx(x)*uex(x,t)-(c(x))^2*uexx(x,t);
        %z=(pi^2-1)*cos(pi*x)*cos(t);
        return
    end
%% Solution at t=0
    function z = f( x )
        z = cos(pi*x);
        return
    end
    function z = fp( x )
        z=-pi*sin(pi*x);
        return
    end

    function z=fpp( x )
        z=-pi^2*cos(pi*x);
        return
    end

    function z = g( x )
        % time derivative of solution at t=0
        z = 0;
        return
    end

%% Wave Speed
    function z=c(x)
        z=cos(x);
        %z=1;
        return
    end

    function z=cx(x)
        z=-sin(x);
        return
    end

%% Plot function
    function plotStuff( x,u,ia,ib,t )
        
        u_exact = ue( x,t );
        
        ind = ia:ib;
        subplot( 2,1,1 )
        plot( x(ind),u(ind),'rx', x(ind),u_exact(ind),'k-' );
        xlabel( 'x' );
        ylabel( 'u' );
        legend('Numerical Solution','Exact Solution')
        
        subplot(2,1,2)
        plot( x(ind),u(ind)-u_exact(ind),'rx' );
        xlabel( 'x' );
        ylabel( 'error' );
        drawnow;
        pause
        
        return
    end

%% Initial Condition setter
    function [unm1,un] = setICs( x,ia,ib,dx,dt )
        %% set intial condition at t=0 (and then BCs on initial condition)
        unm1 = zeros(size(x));
        for i = ia:ib
            unm1(i) = f( x(i) );
        end
        %plot(ia:ib,unm1)
        unm1 = setBCs(  unm1,ia,ib );
        
        %% now set initial condition at t=dt (and then BCs)
        %% here we use a Taylor expansion in time
        un = zeros(size(x));
        for i = ia:ib
            i_plus=(i-1)*dx;
            i_minus=(i-2)*dx;
            
            f0  = f(x(i));
            
            g0  = g(x(i));
            h0  = h(x(i),0);
            uxx = 2*c(x(i))*cx(x(i))*fp(x(i))+(c(x(i)))^2*fpp(x(i));
            %uxx = ((c(i_plus))^2*(unm1(i+1)-unm1(i))-...
            %(c(i_minus))^2*(unm1(i)-unm1(i-1)))/(dx^2);
            utt = uxx+h0;
            
            un(i) = f0+dt*g0+0.5*dt^2*(utt);
            
        end
        un = setBCs(  un,ia,ib);
        return
    end
end