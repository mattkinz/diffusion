clc 
clear
close all
ti=0; %time interval
tf=1.5;
L=1;   %length of row
N=400;  %number of points on the rod in both x and y directions
M=8000;  %number of time points
k=.01; %diffusibity constant

x=linspace(0,L,N);
dx=x(2)-x(1);
t=linspace(ti,tf,M);
dt=t(2)-t(1);

disp(k*dt/dx^2)

u=zeros(M,N); %rows over time, columns over space
u(1,:)=sin(pi*x)+.2*sin(10*pi*x); %initial conditions


u(:,1) = 0; %boundary conditions 
u(:,N)= 0;
r=k*dt/(dx^2);  %r needs to be less than .5

tic;
for j=1:M-1
    
    for n=2:N-1
        u(j+1,n)= u(j,n) + r * (u(j,n+1) - 2*u(j,n) + u(j,n-1));
    end 
    
end
toc


%exact solution
% ue=zeros(M,N);
% for j=1:N
%     for i=1:M
%     ue(i,j)= exp(-pi^2*k*t(i))*sin(pi*x(j))+0.2*exp(-(10*pi)^2*k*t(i))*sin(10*pi*x(j));
%     end
% end



[X,T]=meshgrid(x,t);

ue = exp(-pi^2*k*T).*sin(pi*X)+0.2*exp(-(10*pi)^2*k*T).*sin(10*pi*X);
error = abs(u-ue);
figure(2)

pcolor(X,T,error); shading interp; colorbar;
colormap jet
set(gca,'FontSize',18)
ylabel('time','interpreter','latex','FontSize',32)
xlabel('space','interpreter','latex','FontSize',32)
title(sprintf('r = %0.3f',r),'interpreter','latex','fontsize',26)
set(gcf,'position',[403   210   435   440])
