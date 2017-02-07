% set number of grid points
N  = 7; % note this is off by one from the notes

% c*dt/dx
CFL = 0.8;

% flag to indicate BC type
iBC = 1;

Ntot = N+2;

A = zeros(Ntot,Ntot);
B = zeros(Ntot,Ntot);

% interior discretization
for i = 2:Ntot-1
  A(i,i+1) = 1;
  A(i,i)   = -2;
  A(i,i-1) = 1;
  
  B(i,i) = 1;
end

% BCa
if( iBC == 1 )
  B(1,1) = -1;
  B(1,3) = 1;
  
  B(Ntot,Ntot)   = 1;
  B(Ntot,Ntot-2) = -1;
  
  A(1,1) = -1;
  A(1,3) = 1;
  
  A(Ntot,Ntot)   = 1;
  A(Ntot,Ntot-2) = -1;
elseif( iBC == 2 )
  B(1,1) = 1;
  B(1,2) = 3/2;
  B(1,3) = -6/2;
  B(1,4) = 1/2;
  
  B(Ntot,Ntot)   = 1;
  B(Ntot,Ntot-1) = 3/2;
  B(Ntot,Ntot-2) = -6/2;
  B(Ntot,Ntot-3) = 1/2;
  
  A(1,1) = 1;
  A(1,2) = 3/2;
  A(1,3) = -6/2;
  A(1,4) = 1/2;
  
  A(Ntot,Ntot)   = 1;
  A(Ntot,Ntot-1) = 3/2;
  A(Ntot,Ntot-2) = -6/2;
  A(Ntot,Ntot-3) = 1/2;
 
elseif( iBC == 3 )
  B(1,1) = 1;
  B(1,3) = 1;
  
  B(Ntot,Ntot)   = 1;
  B(Ntot,Ntot-2) = 1;
  
  A(1,1) = 1;
  A(1,3) = 1;
  
  A(Ntot,Ntot)   = 1;
  A(Ntot,Ntot-2) = 1;
elseif( iBC == 4 )
  B(1,1) = 1;
  B(1,2) = 3;
  B(1,3) = -1;
  
  B(N,N)   = -1;
  B(N,N-1) = 3;
  B(N,N-2) = 1;
  
  A(1,1) = 1;
  A(1,2) = 3;
  A(1,3) = -1;
  
  A(N,N)   = -1;
  A(N,N-1) = 3;
  A(N,N-2) = 1;
end
A
B
pause

[V,D] = eig(A,B);
for iv = 1:Ntot
  VV = V(:,iv);
  plot( [1:Ntot]',real(VV),'r-x', [1:Ntot]',imag(VV),'b-s', [1:Ntot]',abs(VV),'k-d' );
  legend( 'Re', 'Imag', 'abs' );
  sigma = D(iv,iv);
  b = (2+CFL^2*sigma)/2;
  ap = b+sqrt(b^2-1);
  am = b-sqrt(b^2-1);
  fprintf( '%i: ap=%s, am=%s,  |a|=%e\n', iv, num2str(ap), num2str(am), norm(ap) );
  
  pause
end

