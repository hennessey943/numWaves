N0 = 20;

xa = 0;
xb = pi;
ya=0;
yb=pi;
c = 1;
tf = 1;
m = 4;
iOrder=4;
hPlot   = zeros(m,1);
errPlot = zeros(m,1);
kPlot=zeros(m,1);
cfl=.8;

for k = 0:m-1
    N = N0*2^k;
    M=(N0)*2^k;
   [errPlot(k+1),hPlot(k+1),kPlot(k+1)]=waveD2wO4(M,N,xa,xb,ya,yb,c,cfl,tf,0,iOrder);

end

loglog( hPlot,errPlot,'x', hPlot,(hPlot.^4)/30, 'g-');
xlabel('dx')
ylabel('error')
legend('expected error','actual error')
figure
loglog(kPlot,errPlot,'x',kPlot,(kPlot.^4)/40,'k-');
xlabel('dy')
ylabel('error')
legend('expected error','actual error')