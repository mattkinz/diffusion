clc
clear
close all

L=1; %rod
N=200; %point on rod
T=1.5; %time interval
M=2000; %points in time
k=.01; %diffusivity constant

x=linspace(0,L,N);
dx = x(2) - x(1);
t=linspace(0,T,M);
dt=t(2)-t(1);
u=zeros(N,M);   %rows through space, columns through time
K1=zeros(N-2,1);
K2=zeros(N-2,1);
Utilda=zeros(N,1);
%initial conditions
u(:,1)=sin(pi*x)+.2*sin(10*pi*x); 
%boundary conditions
u(1,:) = 0;  
u(N,:)= 0;

%loop through time at each point to get the slope at that point
for j=1:M-1
    K1 = (k/dx^2) * (u(1:N-2,j) - 2*u(2:N-1,j) + u(3:N,j)); %slope at point xi at time j
    Utilda(2:N-1) = u(2:N-1,j) + dt * K1; %approximation for value at point xi at time j+1
    K2 = (k/dx^2) * (Utilda(1:N-2) - 2*Utilda(2:N-1) + Utilda(3:N)); %slope at next time step
    u(2:N-1,j+1) = u(2:N-1,j) + dt*((K1+K2)/2); %improved euler's
end


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
pcolor(X,T,error'); shading interp; colorbar;
colormap jet
set(gca,'FontSize',18)
ylabel('time','interpreter','latex','FontSize',32)
xlabel('space','interpreter','latex','FontSize',32)
title(sprintf('k = %0.3f',k),'interpreter','latex','fontsize',26)
set(gcf,'position',[403   210   435   440])
