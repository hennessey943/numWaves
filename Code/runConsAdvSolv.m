% Runs consAdvSolv
a=2;
N=100;
xa=-1;
xb=1;
cfl=.8;
tf=.3;
order=4;
iPlot=1;
[~]=consAdvSolvV2(a,N,xa,xb,cfl,tf,order,iPlot);

%% Convergence Test
if 1==12
m=3;
N0=100;
dxplot=zeros(1,m+1);
errplot=zeros(1,m+1);
errplot1=zeros(1,m+1);
for m=0:3
    [dxplot(m+1),errplot(m+1),errplot1(m+1),~,~]=...
        consAdvSolvV2(a,N0*2^m,xa,xb,cfl,tf,order,0);
end

loglog(dxplot,dxplot.^order,'b-',dxplot,errplot,'rx')
legend('4th Order','Numerical')
figure
loglog(dxplot,dxplot.^(order/(order+1)),'b-',dxplot,errplot1,'rx')
legend('4/5 Order','Numerical')
end

%% Richardson Test
if 2==21
    order=4
    N0=350;
    alpha=2;
    N1=N0;
    N2=N0*alpha;
    N3=N0*alpha^2;
    [dx1,~,~,data1,x1]=consAdvSolvV2(a,N1,-1,1,.8,.2,order,0);
    [dx2,~,~,data2,x2]=consAdvSolvV2(a,N2,-1,1,.8,.2,order,0);
    [dx3,~,~,data3,x3]=consAdvSolvV2(a,N3,-1,1,.8,.2,order,0);
    if abs(dx1-dx2*alpha)>1e-13
        fprintf('dx1 and dx2 wrong')
        return
    end
    if abs(dx2-dx3*alpha)>1e-13
        fprinf('dx2 and dx3 wrong')
        return
    end
        
    norm1=0;
    norm2=0;
    mnorm1=0;
    mnorm2=0;
    for i=1:N1
        for j=1:N2
            if abs(x1(i)-x2(j))<1e-13
                norm1=dx1*abs(data1(i)-data2(j))+norm1;
                mnorm1=max(abs(data1(i)-data2(j)),mnorm1);
            end
        end
    end
    
    for i=1:N1
        for j=1:N2
            for k=1:N3
                if abs(x2(j)-x3(k))<1e-13 && abs(x1(i)-x3(k))<1e-13
                    norm2=dx1*abs(data2(j)-data3(k))+norm2;
                    mnorm2=max(abs(data2(j)-data3(k)),mnorm2);
                end
            end
        end
    end
    p=log(norm2/norm1)/log(1/alpha)
    pmax=log(mnorm2/mnorm1)/log(1/alpha)
end