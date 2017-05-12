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

ue=zeros(M,N);
for j=1:N
    for i=1:M
    ue(i,j)= exp(-pi^2*k*t(i))*sin(pi*x(j))+0.2*exp(-(10*pi)^2*k*t(i))*sin(10*pi*x(j));
    end
end

figure(2)

[X,T]=meshgrid(x,t);
pcolor(X,T,ue); shading interp; colorbar;
colormap jet
set(gca,'FontSize',18)
ylabel('time','interpreter','latex','FontSize',32)
xlabel('space','interpreter','latex','FontSize',32)
title(sprintf('k = %0.3f',k),'interpreter','latex','fontsize',26)
set(gcf,'position',[403   210   435   440])