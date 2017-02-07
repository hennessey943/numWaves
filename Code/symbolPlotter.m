%Initialize

y=linspace(-pi/2,pi/2);
u2=zeros(1,numel(y));
u4=zeros(1,numel(y));
u6=zeros(1,numel(y));
u8=zeros(1,numel(y));
u10=zeros(1,numel(y));

v2=zeros(1,numel(y));
v4=zeros(1,numel(y));
v6=zeros(1,numel(y));
v8=zeros(1,numel(y));
v10=zeros(1,numel(y));

%Calculate numerical values for symbols
for j=1:numel(y)
u2(j)=1i*sin(2*y(j));
u4(j)=1i*sin(2*y(j))*(1+2*(sin(y(j)))^2/3);
u6(j)=1i*sin(2*y(j))*(1+2*(sin(y(j)))^2/3+8*(sin(y(j)))^4/15);
u8(j)=1i*sin(2*y(j))*(1+2*(sin(y(j)))^2/3+8*(sin(y(j)))^4/15+...
        16*(sin(y(j)))^6/35);
u10(j)=1i*sin(2*y(j))*(1+2*(sin(y(j)))^2/3+8*(sin(y(j)))^4/15+...
        16*(sin(y(j)))^6/35+128*(sin(y(j)))^8/315);
    
v2(j)=-4*(sin(y(j)))^2;
v4(j)=-4*(sin(y(j)))^2*(1+(sin(y(j)))^2/3);
v6(j)=-4*(sin(y(j)))^2*(1+(sin(y(j)))^2/3+8*(sin(y(j)))^4/45);
v8(j)=-4*(sin(y(j)))^2*(1+(sin(y(j)))^2/3+8*(sin(y(j)))^4/45+4*(sin(y(j)))^6/35);
v10(j)=-4*(sin(y(j)))^2*(1+(sin(y(j)))^2/3+8*(sin(y(j)))^4/45+4*(sin(y(j)))^6/35+128*(sin(y(j)))^8/1575);

end


%Calculate error
real_err2=abs(real(u2)-real(1i*2*y));
real_err4=abs(real(u4)-real(1i*2*y));
real_err6=abs(real(u6)-real(1i*2*y));
real_err8=abs(real(u8)-real(1i*2*y));
real_err10=abs(real(u8)-real(1i*2*y));
imag_err2=abs(imag(u2)-imag(1i*2*y));
imag_err4=abs(imag(u4)-imag(1i*2*y));
imag_err6=abs(imag(u6)-imag(1i*2*y));
imag_err8=abs(imag(u8)-imag(1i*2*y));
imag_err10=abs(imag(u10)-imag(1i*2*y));

vreal_err2=abs(real(v2)-real(-4*y.^2));
vreal_err4=abs(real(v4)-real(-4*y.^2));
vreal_err6=abs(real(v6)-real(-4*y.^2));
vreal_err8=abs(real(v8)-real(-4*y.^2));
vreal_err10=abs(real(v10)-real(-4*y.^2));
vimag_err2=abs(imag(v2)-imag(-4*y.^2));
vimag_err4=abs(imag(v4)-imag(-4*y.^2));
vimag_err6=abs(imag(v6)-imag(-4*y.^2));
vimag_err8=abs(imag(v8)-imag(-4*y.^2));
vimag_err10=abs(imag(v10)-imag(-4*y.^2));


%Plot
plot(y,real(u2),y,real(u4),y,real(u6),y,real(u8),y,real(u10),y,real(1i*2*y))
title('Real Symbol')
xlabel('\eta')
ylabel('Symbol')

figure
plot(y,imag(u2),y,imag(u4),y,imag(u6),y,imag(u8),y,imag(u10),y,imag(1i*2*y))
title('Imaginary Symbol')
xlabel('\eta')
ylabel('Symbol')
figure

plot(y,real_err2,y,real_err4,y,real_err6,y,real_err8,y,real_err10)
title('Real Error')
xlabel('\eta')
ylabel('error')
figure
plot(y,imag_err2,y,imag_err4,y,imag_err6,y,imag_err8,y,imag_err10)
title('Imaginary Error')
xlabel('\eta')
ylabel('error')

figure
plot(y,real(v2),y,real(v4),y,real(v6),y,real(v8),y,real(v10),y,-4*y.^2)
title('Real Symbol')
xlabel('\eta')
ylabel('Symbol')
figure
plot(y,imag(v2),y,imag(v4),y,imag(v6),y,imag(v8),y,imag(v10),y,imag(-4*y.^2))
title('Imaginary Symbol')
xlabel('\eta')
ylabel('Symbol')
figure
plot(y,vreal_err2,y,vreal_err4,y,vreal_err6,y,vreal_err8,y,vreal_err10)
title('Real Error')
xlabel('\eta')
ylabel('error')
figure
plot(y,vimag_err2,y,vimag_err4,y,vimag_err6,y,vimag_err8,y,vimag_err10)
title('Imaginary Error')
xlabel('\eta')
ylabel('error')


