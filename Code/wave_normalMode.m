% set number of grid points
N  = 7; % note this is off by one from the notes

% c*dt/dx
CFL = 0.8;

% flag to indicate BC type
iBC = 1;

A = zeros(N,N);

% interior discretization
for i = 2:N-1
  A(i,i+1) = 1;
  A(i,i)   = -2;
  A(i,i-1) = 1;
end

% discretization plus BC
if( iBC == 1 )
  A(1,1) = -2+0;
  A(1,2) = 1+1;
  
  A(N,N)   = -2+0;
  A(N,N-1) = 1+1;
elseif( iBC == 2 )
  A(1,1) = -2-3/2;
  A(1,2) = 1+6/2;
  A(1,3) = 0-1/2;
  
  A(N,N)   = -2-3/2;
  A(N,N-1) = 1+6/2;
  A(N,N-2) = 0-1/2;
elseif( iBC == 3 )
  A(1,1) = -2;
  A(1,2) = 1-1;
  
  A(N,N)   = -2;
  A(N,N-1) = 1-1;
elseif( iBC == 4 )
  A(1,1) = -2;
  A(1,2) = 1-3;
  A(1,3) = 0+1;
  
  A(N,N)   = -2;
  A(N,N-1) = 1-3;
  A(N,N-2) = 0+1;
end
A
pause

[V,D] = eig(A);
for iv = 1:N
  VV = V(:,iv);
  plot( [1:N]',real(VV),'r-x', [1:N]',imag(VV),'b-s', [1:N]',abs(VV),'k-d' );
  legend( 'Re', 'Imag', 'abs' );
  sigma = D(iv,iv);
  b = (2+CFL^2*sigma)/2;
  ap = b+sqrt(b^2-1);
  am = b-sqrt(b^2-1);
  fprintf( '%i: ap=%s, am=%s,  |a|=%e\n', iv, num2str(ap), num2str(am), norm(ap) );
  
  pause
end

