%run nonconstant wave solver
xa=0;
xb=1;
N0=20;
tf=1;
iPlot=0;
[dx,err]=nonConstWave( N0,xa,xb,tf,iPlot );
m=4;
hPlot   = zeros(m,1);
errPlot = zeros(m,1);

for k = 0:m-1
  N = N0*2^k;  
  [hPlot(k+1),errPlot(k+1)] = nonConstWave( N,xa,xb,tf,0);
end

errPlot
loglog( hPlot,errPlot,'x', hPlot,(hPlot).^2, 'k-' );
legend('Numerical Result','Exact Result')