clc 
clear
close all
ti=0; %time interval
tf=1.5;
L=1;   %length of row
N=200;  %number of points on the rod in both x and y directions
M=2000;  %number of time points
k=.01; %diffusibity constant

x=linspace(0,L,N);
dx=x(2)-x(1);
t=linspace(ti,tf,M);
dt=t(2)-t(1);
u=zeros(M,N); %rows over time, columns over space
u(1,:)=sin(pi*x)+.2*sin(10*pi*x); %initial conditions
%v=zeros(N,N);
%u=zeros(N,N,t);
%u(1,1,:)=0;
%plot3(u);

%u(:,1) = sin((2*pi/tf)*t); %boundary conditions 
%u(:,N)= cos((2*pi/tf)*t);
u(:,1) = 0;
u(:,N)= 0;
r=k*dt/(dx^2);  %r needs to be less than .5

tic;
for j=1:M-1
    
    for n=2:N-1
        u(j+1,n)= u(j,n) + r * (u(j,n+1) - 2*u(j,n) + u(j,n-1));
    end 
    %u(t+1,1)= u(t,1) + r * (u(t,2) - 2*u(t,1) + u(t,N));
    %u(t+1,N)= u(t,N) + r * (u(t,1) - 2*u(t,N) + u(t,N-1));
    %u(t+1,N)=u(t+1,1);
end
toc
%for k=[1 20 150 400 1000 1500 M]
 %   plot(X,u(k,:)) 
 %   ylim ([-5 5])
 %   title (sprintf('Plot of heat distribution at time %3.0f.',k))
 %   pause(2)
%end
%plot(X,u(1,:),'b',X,u(20,:),'r',X,u(150,:),'o',X,u(400,:),'y',X,u(1000,:),'g',X,u(1500,:),'c',X,u(M,:),'k')

%legend('t=1','t=20','t=150','t=400','t=1000','t=1500','t=2000')



figure(2)
%surf(u); shading interp; colorbar;

[X,T]=meshgrid(x,t);
pcolor(X,T,u); shading interp; colorbar;
colormap jet
set(gca,'FontSize',18)
ylabel('time','interpreter','latex','FontSize',32)
xlabel('space','interpreter','latex','FontSize',32)
title(sprintf('k = %0.3f',k),'interpreter','latex','fontsize',26)
set(gcf,'position',[403 210 415 440])
