clc 
clear
close all
ti=0; %time interval
tf=.1;
L=1;   %length of row
N=100;  %number of points on the rod

k=.01; %diffusibity constant
%M=linspace(2000,100000,99);
M=round(linspace(N^2*k*2*tf,5000,500));

Utilda=zeros(N,1);



x=linspace(0,L,N);
dx=x(2)-x(1);
K1=zeros(1,N-2);

ue=zeros(N,1);
error=zeros(N,numel(M)); %error for each version of M, over all points N
error_norm = zeros(1,numel(M)); %an entry for each choice of time points, M

for i=1:numel(M)
    t=linspace(ti,tf,M(i));
    dt=t(2)-t(1);
    u=zeros(N,M(i)); %rows over time, columns over space
    u(:,1)=sin(pi*x)+.2*sin(10*pi*x); %initial conditions
    
    %exact solutions
    for p=1:N
        ue(p)=exp(-pi^2*k*t(M(i)))*sin(pi*x(p))+0.2*exp(-(10*pi)^2*k*t(M(i)))*sin(10*pi*x(p));
    end
    

    %boundary conditions 

    u(1,:) = 0;
    u(N,:)= 0;
    

    for j=1:(M(i)-1)
    
        K1 = (k/dx^2) * (u(1:N-2,j) - 2*u(2:N-1,j) + u(3:N,j)); %slope at point xi at time j
        Utilda(2:N-1) = u(2:N-1,j) + dt * K1; %approximation for value at point xi at time j+1
        K2 = (k/dx^2) * (Utilda(1:N-2) - 2*Utilda(2:N-1) + Utilda(3:N)); %slope at next time step
        u(2:N-1,j+1) = u(2:N-1,j) + dt*((K1+K2)/2); %improved euler's
    end
    error(:,i)= (ue - u(:,M(i)));
    error_norm(i) = norm(ue - u(:,M(i)),inf);
    %error_norm(i) = norm(error(:,i),inf); 
    
end
plot(tf./M*k/dx^2,error_norm,'*-')
set(gca,'FontSize',18)
xlabel('$r=\frac{k\Delta t}{\Delta x^2}$','interpreter','latex','FontSize',24)
ylabel('Maximum error','interpreter','latex','FontSize',24)
legend(sprintf('N = %d',N))
title(sprintf('$k$ = %0.3f',k),'interpreter','latex','fontsize',26)
set(gcf,'position',[403   210   665   440])

