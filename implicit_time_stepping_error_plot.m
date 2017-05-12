clc
clear
close all

L=1; %rod
N=100; %point on rod
T=1.5; %time interval
M=2000; %points in time
k=0.01; %diffusivity constant

x=linspace(0,L,N);
dx = x(2) - x(1);
t=linspace(0,T,M);
dt=t(2)-t(1);

r=k*dt/dx^2;
u=zeros(N,M); %rows through space, columns through time

%boundary conditions
u(1,:)=0;
u(N,:)=0; 


%build the matrix A
A=zeros(N-2,N-2);
A(1,1)= 1+2*r;
A(1,2)= -r;
A(N-2,N-3)= -r;
A(N-2,N-2)= 1+2*r;
tic
for i=2:N-3
    A(i,i-1)= -r;
    A(i,i)= 1+2*r;
    A(i,i+1)= -r;
end
toc
%initial conditions
u(:,1)=sin(pi*x)+.2*sin(10*pi*x); 

b = u(2:N-1,1);
tic
for q=1:M-1
    b(1) = b(1) + r * u(1,q+1);
    b(N-2) = b(N-2) + r * u(N,q+1);
    
    btilda = A\b;
    u(2:N-1,q+1) = btilda;
    b = btilda;
end
toc


%exact solution
ue=zeros(N,M);
for j=1:N
    for i=1:M
    ue(j,i)= exp(-pi^2*k*t(i))*sin(pi*x(j))+0.2*exp(-(10*pi)^2*k*t(i))*sin(10*pi*x(j));
    end
end
error = abs(u-ue);

figure(2)

[X,T]=meshgrid(x,t);
pcolor(X,T,(error')); shading interp; colorbar;
colormap jet
set(gca,'FontSize',18)
ylabel('time','interpreter','latex','FontSize',32)
xlabel('space','interpreter','latex','FontSize',32)
title(sprintf('k = %0.3f',k),'interpreter','latex','fontsize',26)
set(gcf,'position',[403   210   435   440]) 