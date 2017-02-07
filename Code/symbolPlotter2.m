%Initialize

y=linspace(-pi,pi);
u1=zeros(1,numel(y));
u2=zeros(1,numel(y));
u3=zeros(1,numel(y));
u4=zeros(1,numel(y));

v2=zeros(1,numel(y));
v4=zeros(1,numel(y));
v6=zeros(1,numel(y));
v8=zeros(1,numel(y));
v10=zeros(1,numel(y));

%Calculate numerical values for symbols
for j=1:numel(y)
u1(j)=-(1-exp(-1i*y(j)));
u2(j)=-.25*(exp(1i*y(j))+3-5*exp(-1i*y(j))+exp(-2*1i*y(j)));
u3(j)=-1/6*(2*exp(1i*y(j))+3-6*exp(-1i*y(j))+exp(-2*1i*y(j)));
u4(j)=-1/24*(-exp(2*1i*y(j))+11*exp(1i*y(j))+10-26*exp(-1i*y(j))+7*exp(-2*1i*y(j))-exp(-3*1i*y(j)));
end


%Calculate error
real_err1=abs(real(u1)-real(-1i*y));
real_err2=abs(real(u2)-real(-1i*y));
real_err3=abs(real(u3)-real(-1i*y));
real_err4=abs(real(u4)-real(-1i*y));
imag_err1=abs(imag(u1)-imag(-1i*y));
imag_err2=abs(imag(u2)-imag(-1i*y));
imag_err3=abs(imag(u3)-imag(-1i*y));
imag_err4=abs(imag(u4)-imag(-1i*y));

%Plot
plot(y,real(u1),y,real(u2),y,real(u3),y,real(u4),y,real(-1i*y))
title('Real Symbol')
xlabel('\xi')
ylabel('Symbol')
legend('O1','O2','O3','O4','Exact')

figure
plot(y,imag(u1),y,imag(u2),y,imag(u3),y,imag(u4),y,imag(-1i*y))
title('Imaginary Symbol')
xlabel('\xi')
ylabel('Symbol')
legend('O1','O2','O3','O4','Exact')
figure
hold on
plot(u1)
plot(u2)
plot(u3)
plot(u4)
xlabel('Real')
ylabel('Imaginary')
legend('O1','O2','O3','O4')