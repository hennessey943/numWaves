function testExact2D(Case,iPlot)
M=1001;
if Case==1
    %Shock Right, Rare Left
    wL=[1;0;0;1];
    wR=[.125;0;1;.1];
    t=.5;
    x=linspace(-1,1,M);
elseif Case==2
    %Double rare
    wL=[1;-2;-3;.4];
    wR=[1;2;3;.4];
    t=.15;
    x=linspace(-2,2,M);
elseif Case==3
    %rare left shock right
    wL=[1;0;0;1000];
    wR=[1;0;2;.01];
    t=1;
    x=linspace(-100,100,M);
elseif Case==4
    %Shock left rare right
    wL=[1;0;2;.01];
    wR=[1;0;1;100];
    t=1.5;
    x=linspace(-30,30,M);
elseif Case==5
    %Double Shock
    wL=[5.99924;19.5975;0;460.894];
    wR=[5.99242;-6.19633;0;46.0950];
    t=2;
    x=linspace(-30,30,M);
elseif Case==6
    %shock left, rare right
    wL=[1;0;1;1];
    wR=[1.5;0;0;2];
    t=.15;
    x=linspace(-1,1,M);
end
gamma=1.4;
Qtol=1e-8;
[w]=exactEuler2D(wL,wR,gamma,x,t,Qtol);
[uApprox]=HLLCEuler2D(wL,wR,gamma,x,t);
wApprox=convertVar(uApprox,gamma,M);
if iPlot==1
    figure
    plot(x,w(1,:),x,wApprox(1,:))
    title('Density')
    xlabel('x')
    ylabel('\rho')
    legend('Exact','HLLC')
    figure
    plot(x,w(2,:),x,wApprox(2,:))
    title('Normal Velocity')
    xlabel('x')
    ylabel('u')
    legend('Exact','HLLC')
    figure
    plot(x,w(3,:),x,wApprox(3,:))
    title('Tangential Velocity')
    xlabel('x')
    ylabel('v')
    legend('Exact','HLLC')
    figure
    plot(x,w(4,:),x,wApprox(4,:))
    title('Pressure')
    xlabel('x')
    ylabel('p')
    legend('Exact','HLLC')
end
end
%% Primitive Variable Converter
function [w]=convertVar(u,gamma,Mtot)
%density
rho=u(1,:);
%momenta
mx=u(2,:);
my=u(3,:);
%energy
E=u(4,:);
w(1,:)=rho;
for i=1:Mtot
    %Velocities
    w(2,i)=mx(i)/rho(i);
    w(3,i)=my(i)/rho(i);
    %Pressure
    w(4,i)=(gamma-1)*(E(i)-...
        .5*(mx(i)^2+my(i)^2)/rho(i));
end
end