N0 = 20;

xa = 0;
xb = 1;
c = 1;
tf = 1;
m = 4;

hPlot   = zeros(m,1);
errPlot = zeros(m,1);

for k = 0:m-1
  N = N0*2^k;  
  [hPlot(k+1),errPlot(k+1)] = waveWithSourceO4( N,xa,xb,c,tf,0);
end

loglog( hPlot,errPlot,'x', hPlot,(hPlot.^4)/2, 'k-' );