function convExact2D(Conv,iPlot)
if Conv==1
    M0=10;
    N=2;
    ymin=-1;
    ymax=1;
    gamma=1.4;
    Qtol=1e-8;
    wL=[1;0;0;1];
    wR=[.125;0;1;.1];
    t=.5;
    xmin=-1;
    xmax=1;
    m=6;
    dx=zeros(1,m);
    errDen1=zeros(1,m);
    errDen2=zeros(1,m);
    errU1=zeros(1,m);
    errU2=zeros(1,m);
    errV1=zeros(1,m);
    errV2=zeros(1,m);
    errP1=zeros(1,m);
    errP2=zeros(1,m);
    for k=0:m-1
        M=M0*2^k+1;
        uInit=eulerInit(8,xmin,xmax,ymin,ymax,M,N,gamma);
        x=linspace(xmin,xmax,M);
        dx(k+1)=(xmax-xmin)/(M-1);
        dy=(ymax-ymin)/(N-1);
        [~,pdata1]=numEuler2D(uInit,gamma,...
            dx(k+1),dy,.8,t,1,1,1);
        [~,pdata2]=numEuler2D(uInit,gamma,...
            dx(k+1),dy,.8,t,2,1,1);
        [w]=exactEuler2D(wL,wR,gamma,x,t,Qtol);
        for i=1:M
            errDen1(k+1)=abs(pdata1(i,1,1)-w(1,i))*dx(k+1)+...
                errDen1(k+1);
            errDen2(k+1)=abs(pdata2(i,1,1)-w(1,i))*dx(k+1)+...
                errDen2(k+1);
            errU1(k+1)=abs(pdata1(i,1,2)-w(2,i))*dx(k+1)+...
                errU1(k+1);
            errU2(k+1)=abs(pdata2(i,1,2)-w(2,i))*dx(k+1)+...
                errU2(k+1);
            errV1(k+1)=abs(pdata1(i,1,3)-w(3,i))*dx(k+1)+...
                errV1(k+1);
            errV2(k+1)=abs(pdata2(i,1,3)-w(3,i))*dx(k+1)+...
                errV2(k+1);
            errP1(k+1)=abs(pdata1(i,1,4)-w(4,i))*dx(k+1)+...
                errP1(k+1);
            errP2(k+1)=abs(pdata2(i,1,4)-w(4,i))*dx(k+1)+...
                errP2(k+1);
        end
        if iPlot==1
            figure
            plot(x,w(1,:),x,pdata1(:,1,1),x,pdata2(:,1,1))
            title('Density')
            xlabel('x')
            ylabel('\rho')
            legend('Exact','Exact Flux','HLLC Flux')
            figure
            plot(x,w(2,:),x,pdata1(:,1,2),x,pdata2(:,1,2))
            title('Transverse Velocity')
            xlabel('x')
            ylabel('u')
            legend('Exact','Exact Flux','HLLC Flux')
            figure
            plot(x,w(3,:),x,pdata1(:,1,3),x,pdata2(:,1,3))
            title('Normal Velocity')
            xlabel('x')
            ylabel('v')
            legend('Exact','Exact Flux','HLLC Flux')
            figure
            plot(x,w(4,:),x,pdata1(:,1,4),x,pdata2(:,1,4))
            title('Pressure')
            xlabel('x')
            ylabel('p')
            legend('Exact','Exact Flux','HLLC Flux')
        end
    end
    figure
    loglog(dx,errDen1,'ro',dx,errDen2,'bx',dx,dx.^.5/10,'g-')
    title('Density Convergence')
    xlabel('log(dx)')
    ylabel('log(error)')
    legend('Exact Flux Error','HLLC Flux Error','1/2 Order Convergence Reference')
elseif Conv==2
    gamma=1.4;
    M0=40;
    N0=40;
    M1=80;
    N1=80;
    M2=160;
    N2=160;
    xmin=-5;
    xmax=5;
    ymin=-5;
    ymax=5;
    tF=.5;
    x0=linspace(xmin,xmax,M0);
    x1=linspace(xmin,xmax,M1);
    x2=linspace(xmin,xmax,M2);
    uInit0=eulerInit(7,xmin,xmax,ymin,ymax,M0,N0,gamma);
    uInit1=eulerInit(7,xmin,xmax,ymin,ymax,M1,N1,gamma);
    uInit2=eulerInit(7,xmin,xmax,ymin,ymax,M2,N2,gamma);
    dx0=(xmax-xmin)/(M0-1);
    dy0=(ymax-ymin)/(N0-1);
    dx1=(xmax-xmin)/(M1-1);
    dy1=(ymax-ymin)/(N1-1);
    dx2=(xmax-xmin)/(M2-1);
    dy2=(ymax-ymin)/(N2-1);
    [~,pdata0]=numEuler2D(uInit0,gamma,...
        dx0,dy0,.8,tF,2,1,1);
    [~,pdata1]=numEuler2D(uInit1,gamma,...
        dx1,dy1,.8,tF,2,1,1);
    [~,pdata2]=numEuler2D(uInit2,gamma,...
        dx2,dy2,.8,tF,2,1,1);
    m2=zeros(1,4);
    m1=zeros(1,4);
    for j=1:M2
        for i=1:M1
            if x2(j)==x1(i)
                for k=1:4
                    m2(k)=max(abs(pdata1(i,i,k)-pdata2(i,i,k)),m2(k));
                end
            end
            for l=1:M0
                if x1(i)==x0(l)
                    for k=1:4
                        m1(k)=max(abs(pdata0(l,l,k)-pdata1(l,l,k)),m1(k));
                    end
                end
            end
        end
    end
    order=zeros(1,4);
    for i=1:4
        order(i)=log(m2(i)/m1(i))/log(1/2);
    end
    order
elseif Conv==3 %Twilight zone
    M0=5;
    ymin=-2;
    ymax=2;
    gamma=1.4;
    t=.15;
    xmin=-2;
    xmax=2;
    m=4;
    dx=zeros(1,m);
    err1=zeros(4,m);
    err2=zeros(4,m);
    for k=0:m-1
        M=M0*2^k;
        N=M;
        uInit=eulerInit(1,xmin,xmax,ymin,ymax,M,N,gamma);
        dx(k+1)=(xmax-xmin)/(M-1);
        dy=(ymax-ymin)/(N-1);
        [err1(:,k+1),~,~]=numEuler2D(uInit,gamma,...
           dx(k+1),dy,.8,t,1,1,7);
        [err2(:,k+1),~,~]=numEuler2D(uInit,gamma,...
           dx(k+1),dy,.8,t,2,1,7);
        
        k
    end
    if iPlot==1
        figure
        loglog(dx,err1(1,:),'ro',dx,err2(1,:),'bx',dx,dx,'g-')
        title('Density Convergence')
        xlabel('log(dx)')
        ylabel('log(error)')
        legend('Exact Flux Error','HLLC Flux Error','1st Order Convergence Reference')
        
        figure
        loglog(dx,err1(2,:),'ro',dx,err2(2,:),'bx',dx,dx,'g-')
        title('x-Velocity Convergence')
        xlabel('log(dx)')
        ylabel('log(error)')
        legend('Exact Flux Error','HLLC Flux Error','1st Order Convergence Reference')
        figure
        loglog(dx,err1(3,:),'ro',dx,err2(3,:),'bx',dx,dx,'g-')
        title('y-Velocity Convergence')
        xlabel('log(dx)')
        ylabel('log(error)')
        legend('Exact Flux Error','HLLC Flux Error','1st Order Convergence Reference')
        figure
        loglog(dx,err1(4,:),'ro',dx,err2(4,:),'bx',dx,dx,'g-')
        title('Pressure Convergence')
        xlabel('log(dx)')
        ylabel('log(error)')
        legend('Exact Flux Error','HLLC Flux Error','1st Order Convergence Reference')
    end
end