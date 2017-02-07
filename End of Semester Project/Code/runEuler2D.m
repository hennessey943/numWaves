function runEuler2D(inCase,Case,iBound,iPlot,ex)

xmin=-2;
xmax=2;
ymin=-2;
ymax=2;
M=31;
N=31;
x=linspace(xmin,xmax,M);
y=linspace(ymin,ymax,N);
gamma=1.4;

dx=(xmax-xmin)/(M-1);
dy=(ymax-ymin)/(N-1);
cfl=.8;
tF=1;
nPlot=1;
figure
uInit=eulerInit(inCase,xmin,xmax,ymin,ymax,M,N,gamma);
data=uInit;
%v=VideoWriter('HLLCBoxOpen.avi');
%open(v);


for i=1:nPlot
    tInc=tF/nPlot
    if ex==1
        exact=zeros(M,N,4);
        for j=1:M
            for k=1:N
                exact(j,k,:)=getEx(x(j),y(k),tInc*i);
            end
        end
    end
    [~,data,pdata]=numEuler2D(data,...
        gamma,dx,dy,cfl,tInc,Case,1,iBound);
    err=pdata-exact;
    if iPlot==1
        %figure
        if M~=N
            plot(x,pdata(:,1,1))
            xlabel('x')
            ylabel('\rho')
        else
            %surf(x,y,exdata(:,:,4))
            figure(1)
            surf(x,y,err(:,:,1))
            xlabel('y')
            ylabel('x')
            zlabel('\rho')
%             axis([-2 2 -2 2 0 2])
            
            figure(2)
            surf(x,y,err(:,:,2))
            xlabel('y')
            ylabel('x')
            zlabel('u')
%             axis([-2 2 -2 2 0 2])
            
            figure(3)
            surf(x,y,err(:,:,3))
            xlabel('y')
            ylabel('x')
            zlabel('v')
%             axis ([-2 2 -2 2 0 2])
            
            figure(4)
            surf(x,y,err(:,:,4))
            xlabel('y')
            ylabel('x')
            zlabel('p')
%             axis([-2 2 -2 2 0 2])
            
            drawnow
            pause
        end
        %title('Pressure')
        %F=getframe;
        %writeVideo(v,F)
        %     figure
        %     surf(x,y,pdata(:,:,2))
        %     title('x-Velocity')
        %     xlabel('x')
        %     ylabel('y')
        %     zlabel('u')
        %     figure
        %     surf(x,y,pdata(:,:,3))
        %     title('y-Velocity')
        %     xlabel('x')
        %     ylabel('y')
        %     zlabel('v')
        %     figure
        %     surf(x,y,pdata(:,:,4))
        %     title('Pressure')
        %     xlabel('x')
        %     ylabel('y')
        %     zlabel('p')
        i
    end
end
%close(v);

end

function exact=getEx(x,y,t)
exact(1)=1+.25*sin(pi*x)*sin(pi*y)*sin(t);
exact(2)=1;
exact(3)=.5;
exact(4)=1-.25*cos(pi*x)*cos(pi*y)*sin(t);
end
