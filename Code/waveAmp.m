%% script for the amplification factor of the 2nd order scheme
m = 501;

lambda = linspace(-1.1,1.1,m);
sigma=linspace(-1.1,1.1,m);
xi     = linspace(-pi,pi,m);
rhop   = zeros(m,m);
rhom   = zeros(m,m);

for j = 1:m
  for k = 1:m
    lamL = .75;
    sigL=.7;
    xiL  = xi(j);
    etaL=eta(k);
    b = 1-sigL^2*(1-cos(xiL))-lamL^2*(1-cos(etaL));
    
    rhop(j,k) = abs(b+sqrt(b^2-1));
    rhom(j,k) = abs(b-sqrt(b^2-1));
    
  end
end

figure
fs = 16;
lineWidth = 2;
ms = 16;
set(gca,'FontSize',fs);
surf( xi,lambda,rhop' );
shading interp
xlabel( '\xi' );
ylabel( '\eta' );
zlabel( 'a_{+}' );
title('\sigma=.7, and \lambda=.75');
% plotName = sprintf('images/wave_rhop.eps');
% fprintf('Saving file=[%s]\n',plotName);
% print('-depsc2',plotName);

figure
fs = 16;
lineWidth = 2;
ms = 16;
set(gca,'FontSize',fs);
surf( xi,lambda,rhom' );
shading interp
xlabel( '\xi' );
ylabel( '\eta' );
zlabel( 'a_{-}' );
title('\sigma=.7, and \lambda=.75');
% plotName = sprintf('images/wave_rhom.eps');
% fprintf('Saving file=[%s]\n',plotName);
% print('-depsc2',plotName);