function [dx,err] = waveWithSourceO4( N,xa,xb,c,tf,iPlot )
  %N = 50;
  %xa = 0;
  %xb = 1;
  %c = 1;
  %tf = .7;

  %% number of ghost points
  ng = 2;

  %% grid size and other grid quantities
  dx = (xb-xa)/N;
  dt = 0.9*dx/c;
  Nt = ceil( tf/dt );
  dt = tf/Nt;

  Ntot = N+1+2*ng;

  %% first and last interior grid points
  ia = ng+1;
  ib = Ntot-ng;

  x = linspace( xa-ng*dx, xb+ng*dx, Ntot );

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
    un = setBCs( un,x,ia,ib,dx,dt,t,c );

    %% compute update over domain interior
    for i = ia:ib
        
        uxx = (-un(i-2)+16*un(i-1)-30*un(i)+16*un(i+1)-un(i+2))/(12*dx^2);
        utt=c^2*uxx+h(x(i),t);
%         
        uxxxx=(un(i-2)-4*un(i-1)+6*un(i)-4*un(i+1)+un(i+2))/(dx^4);
        %Note for O(dx^4) accuracy, htt must be approximated to O(dx^6)
        %hxx=-4*h(x(i),t);
        htt=-h(x(i),t);
        %hxx=(h(x(i-1),t)-2*h(x(i),t)+h(x(i+1),t))/dx^2;
        hxx=(-h(x(i-2),t)+16*h(x(i-1),t)-30*h(x(i),t)+16*h(x(i+1),t)-h(x(i+2),t))/(12*dx^2);
        %htt=(h(x(i),t-1)-2*h(x(i),t)+h(x(i),t+1))/dt^2;
        %htt=(-h(x(i),t-2*dt)+16*h(x(i),t-1*dt)-30*h(x(i),t)+16*h(x(i),t+1*dt)-h(x(i),t+2*dt))/(12*dx^2);
        
        utttt=uxxxx+hxx+htt;
  %utt=-uex(x(i),t,c);
  %utttt=uex(x(i),t,c);
      
        u(i) = 2*un(i)-unm1(i)+dt^2*utt+dt^4*utttt/12;

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

function un = setBCs( un,x,ia,ib,dx,dt,t,c)
  

%     un(ia-1)=cos(2*x(ia-1))*cos(t);
%     un(ia-2)=cos(2*x(ia-2))*cos(t);
%     un(ib+1)=cos(2*x(ib+1))*cos(t);
%     un(ib+2)=cos(2*x(ib+2))*cos(t);
    gltt=(-gl(t-2*dt)+16*gl(t-dt)-30*gl(t)+16*gl(t+dt)-gl(t+2*dt))/(12*dt^2);
    gltttt=(gl(t-2*dt)-4*gl(t-dt)+6*gl(t)-4*gl(t+dt)+gl(t+2*dt))/dt^4;
    htt=(h(x(ia),t-dt)-2*h(x(ia),t)+h(x(ia),t+dt))/dt^2;
    hxx=(h(x(ia-1),t)-2*h(x(ia),t)+h(x(ia+1),t))/dx^2;
    un(ia-1) = 2*un(ia)-un(ia+1)...
            +dx^2*(gltt-h(x(ia),t))...
            +dx^4/(12)*(gltttt-htt-hxx);
            
    un(ia-2)=4*un(ia-1)-6*un(ia)+4*un(ia+1)-un(ia+2)...
            +dx^4*(gltttt-htt-hxx);
        
    grtt=(-gr(t-2*dt)+16*gr(t-dt)-30*gr(t)+16*gr(t+dt)-gr(t+2*dt))/(12*dt^2);
    hx=(h(x(ib+1),t)-h(x(ib-1),t))/(2*dx);
    
    un(ib+1) = un(ib-1)+dx^3/3*(grtt-hx)+2*dx*gr(t);
    un(ib+2) = 2*un(ib+1)-2*un(ib-1)+un(ib-2)+2*dx^3*(grtt-hx);
%     un(ib)=gr(t);
%     un(ia)=gl(t);
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
  unm1 = setBCs( unm1,x,ia,ib,dx,dt,0,c );
  
  %% now set initial condition at t=dt (and then BCs)
  %% here we use a Taylor expansion in time
  un = zeros(size(x));
  for i = ia:ib
    fm2 = f(x(i-2));
    fm1 = f(x(i-1));
    f0  = f(x(i));
    fp1 = f(x(i+1));
    fp2 = f(x(i+2));
    
    g0  = g(x(i));
    hm1=h(x(i-1),0);
    htm1=h(x(i),-dt);
    h0  = h(x(i),0);
    hp1=h(x(i+1),0);
    htp1=h(x(i),dt);
    
    
    uxx = (-fm2+16*fm1-30*f0+16*fp1-fp2)/(12*dx^2);
    
    utt = c^2*uxx+h0;
    
    
    uttt=(htp1-htm1)/(2*dt);
    
    uxxxx= (fm2-4*fm1+6*f0-4*fp1+fp2)/dx^4;
    hxx= (hm1-2*h0+hp1)/dx^2;
    htt= (htm1-2*h0+htp1)/dt^2;
    
    utttt=uxxxx+hxx+htt;

    %utt=-uex(x(i),dt,c);
    %utttt=uex(x(i),dt,c);
    %uttt=-cos(2*x(i))*sin(dt);
    
    un(i) = f0+0.5*dt^2*(utt)+dt^3/6*uttt+dt^4/24*utttt;
%un(i)=uex(x(i),dt,c);

  end
  un = setBCs( un,x,ia,ib,dx,dt,dt );
  return
end