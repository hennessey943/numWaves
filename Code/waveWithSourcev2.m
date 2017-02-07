function [dx,err] = waveWithSourcev2( N,xa,xb,c,tf,iPlot )
  %N = 50;
  %xa = 0;
  %xb = 1;
  %c = 1;
  %tf = .7;

  %% number of ghost points
  ng = 1;

  %% grid size and other grid quantities
  dx = (xb-xa)/N;
  dt = 0.9*dx/c;
  Nt = ceil( tf/dt );
  dt = tf/Nt;

  Ntot = N+1+2*ng;

  %% first and last interior grid points
  ia = ng+1;
  ib = Ntot-ng;

  x = linspace( xa-dx, xb+dx, Ntot );

  [unm1,un] = setICs( x,ia,ib,dx,dt,c);
  
  %% allocate space for u
  u = zeros(size(un));

  if( iPlot == 1 )
    figure
  end

  %% Loop through time
  t = dt;
  for n = 1:Nt-1

    %% set BCs on un
    un = setBCs( un,x,ia,ib,dx,dt,t );

    %% compute update over domain interior
    for i = ia:ib
      uxx = (...
            un(i+1)...
        -2.*un(i)...
           +un(i-1))/(dx^2);
       u(i) = 2.*un(i)-unm1(i)+(c*dt)^2*uxx+dt^2*h(x(i),t);
    end

    %% update solution histories
    unm1 = un;
    un   = u;
    t = t+dt;
    
    if( iPlot == 1 )
      plotStuff( x,u,ia,ib,t,c );
    end
  end
    
  ue = uex( x,t,c );

  ind = ia:ib;
  err = max(abs(u(ind)-ue(ind)));
  return
end

function un = setBCs( un,x,ia,ib,dx,dt,t)
  un(ia-1) = 2*un(ia)-un(ia+1)...
            +dx^2/dt^2*(gl(t-dt)-2*gl(t)+gl(t+dt))-dx^2*h(x(ia),t);
  un(ib+1) = un(ib-1)+2*dx*gr(t);
  return
end

function z = uex( x,t,c )
  z = cos(2*x)*cos(t);
  return
end

function z = f( x )
  %% solution at t=0
  z = cos(2*x);
  return
end

function z = g( x )
  %% time derivative of solution at t=0
  z = 0;
  return
end

function z = h( x,t )
    %% Source term
    z = 3*cos(2*x)*cos(t);
    return
end

function z=gl(t)
z=cos(t);
return
end

function z=gr(t)
z=-2*sin(2)*cos(t);
return
end

function plotStuff( x,u,ia,ib,t,c )

  ue = uex( x,t,c );

  ind = ia:ib;
  subplot( 2,1,1 )
  plot( x(ind),u(ind),'rx', x(ind),ue(ind),'k-' );
  xlabel( 'x' );
  ylabel( 'u' );

  subplot(2,1,2)
  plot( x(ind),u(ind)-ue(ind),'rx' );
  xlabel( 'x' );
  ylabel( 'error' );
  drawnow;
  pause
  
  return
end

function [unm1,un] = setICs( x,ia,ib,dx,dt,c )
  %% set intial condition at t=0 (and then BCs on initial condition)
  unm1 = zeros(size(x));
  for i = ia:ib
    unm1(i) = f( x(i) );
  end
  figure
  %plot(ia:ib,unm1)
  unm1 = setBCs( unm1,x,ia,ib,dx,dt,0 );
  
  %% now set initial condition at t=dt (and then BCs)
  %% here we use a Taylor expansion in time
  un = zeros(size(x));
  for i = ia:ib
    fm1 = f(x(i-1));
    f0  = f(x(i));
    fp1 = f(x(i+1));
    
    g0  = g(x(i));
    h0  = h(x(i),0);
    
    uxx = (...
       1.*fp1...
      -2.*f0...
      +1.*fm1)/(dx^2);
    
    utt = c^2*uxx+h0;
    
    un(i) = f0+dt*g0+0.5*dt^2*(utt);

  end
  un = setBCs( un,x,ia,ib,dx,dt,dt );
  return
end