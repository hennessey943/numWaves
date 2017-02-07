%Runs Richardson Extrapolation test on 2D wave solver
N0 = 200;
M=N0;
N=N0;
xa = 0;
xb = pi;
ya=0;
yb=pi;
c = 1;
tf = 3;
iOrder=2;
iPlot=1;
cfl=.8;

%[data,dx,x]=waveD2RichTest(M,N,xa,xb,ya,yb,c,cfl,tf,iPlot,iOrder);

%% Richardson Test
if 1==1
    N0=30;
    alpha=2;
    N1=N0;
    N2=N0*alpha;
    N3=N0*alpha^2;
    [data1,dx1,x1,y1]=waveD2RichTest(N1,N1,xa,xb,ya,yb,c,cfl,tf,0,iOrder);
    [data2,dx2,x2,y2]=waveD2RichTest(N2,N2,xa,xb,ya,yb,c,cfl,tf,0,iOrder);
    [data3,dx3,x3,y3]=waveD2RichTest(N3,N3,xa,xb,ya,yb,c,cfl,tf,0,iOrder);
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
                norm1=dx1*abs(data1(i,i)-data2(j,j))+norm1;
                mnorm1=max(abs(data1(i,i)-data2(j,j)),mnorm1);
            end
        end
    end
    
    for i=1:N1
        for j=1:N2
            for k=1:N3
                if abs(x2(j)-x3(k))<1e-13 && abs(x1(i)-x3(k))<1e-13
                    norm2=dx1*abs(data2(j,j)-data3(k,k))+norm2;
                    mnorm2=max(abs(data2(j,j)-data3(k,k)),mnorm2);
                end
            end
        end
    end
    p=log(norm2/norm1)/log(1/alpha)
    pmax=log(mnorm2/mnorm1)/log(1/alpha)
end