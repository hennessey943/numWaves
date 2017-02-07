rho0=1;
u0=0;
p0=1;
gamma=1.4;
M=500;
p=zeros(1,M);
u1=zeros(1,M);
u2=zeros(1,M);
ur1=zeros(1,M);
ur2=zeros(1,M);
for i=1:M
    p(i)=i/50;
    u1(i)=u0-(p(i)-p0)*...
        ((2/(gamma+1)*rho0)/(p(i)+(gamma-1)/(gamma+1)*p0))^.5;
    u2(i)=u0+(p(i)-p0)*...
        ((2/(gamma+1)*rho0)/(p(i)+(gamma-1)/(gamma+1)*p0))^.5;
    ur1(i)=u0-2/(gamma-1)*sqrt(gamma*p0/rho0)*...
        ((p(i)/p0)^((gamma-1)/(2*gamma))-1);
    ur2(i)=u0+2/(gamma-1)*sqrt(gamma*p0/rho0)*...
        ((p(i)/p0)^((gamma-1)/(2*gamma))-1);
end
hold on
p1=plot(p,u1,'k-');
p2=plot(p,u2,'k-');
p3=plot(p,ur1,'r--');
p4=plot(p,ur2,'r--');
title('Integral Curve and Hugoniot Locus')
xlabel('p_*')
ylabel('u')
legend([p1 p3],'Shock','Rarefaction')
