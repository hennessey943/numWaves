% run FOS solvers
N0=200;
xa=0;
xb=1;
cfl=.7;
tf=1.2;
iPlot=1;
Case=4;
%[dx,v_err,sig_err]=colFOSsolver(N0,xa,xb,cfl,tf,iPlot,Case);
[dx,v_err,sig_err]=stagFOSsolver(N0,xa,xb,cfl,tf,iPlot,Case);

if 0==9
m=4;
hPlot   = zeros(m,1);
v_errPlot = zeros(m,1);
sig_errPlot=zeros(m,1);

for k = 0:m-1
  N = N0*2^k;  
  %[hPlot(k+1),v_errPlot(k+1),sig_errPlot(k+1)] = colFOSsolver(N,xa,xb,cfl,tf,0,1);
  [hPlot(k+1),v_errPlot(k+1),sig_errPlot(k+1)] = stagFOSsolver(N,xa,xb,cfl,tf,0,1);
end

loglog( hPlot,v_errPlot,'x', hPlot,(hPlot).^2, 'k-' );
legend('Numerical Result','Exact Result')
title('v Convergence Test')
figure
loglog(hPlot,sig_errPlot,'x',hPlot,(hPlot).^2,'k-');
legend('Numerical Result','Exact Result')
title('\sigma Convergence Test')
end