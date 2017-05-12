clc 
clear
close all
ti=0; %time interval
tf=.1;
L=1;   %length of row
N=400;  %number of points on the rod in both x and y directions
k=.01; %diffusibity constant


M=round(linspace(N^2*k*2*tf,5000,500));


x=linspace(0,L,N);
dx=x(2)-x(1);

ue=zeros(1,N);
error=zeros(numel(M),N); %error for each version of M, over all points N
error_norm = zeros(1,numel(M)); %rows over number of time points, columns over num of time points

for i=1:numel(M)
    t=linspace(ti,tf,M(i));
    dt=t(2)-t(1);
    u=zeros((M(i)),N); %rows over time, columns over space
    u(1,:)=sin(pi*x)+.2*sin(10*pi*x); %initial conditions
    
    %exact solutions
    %for p=1:N
    %    ue(p)=exp(-pi^2*k*t(.02*M(i)))*sin(pi*x(p))+0.2*exp(-(10*pi)^2*k*t(.02*M(i)))*sin(10*pi*x(p));
    %end
    ue = exp(-pi^2*k*t(end)).*sin(pi*x)+0.2*exp(-(10*pi)^2*k*t(end)).*sin(10*pi*x);
    
    %boundary conditions
    
    u(:,1) = 0;
    u(:,N)= 0;
    r=k*dt/(dx^2);  %r needs to be less than .5
    
    
    for j=1:M(i)-1
        
        for n=2:N-1
            u(j+1,n)= u(j,n) + r * (u(j,n+1) - 2*u(j,n) + u(j,n-1));
        end
    end
    error_norm(i) = norm(ue-u(end,:),inf);
    
end



loglog(tf./M*k/dx^2,error_norm,'*-')
set(gca,'FontSize',18)
xlabel('$r=\frac{k\Delta t}{\Delta x^2}$','interpreter','latex','FontSize',24)
ylabel('Maximum error','interpreter','latex','FontSize',24)
title(sprintf('$k$ = %0.3f',k),'interpreter','latex','fontsize',26)
legend(sprintf('N = %d',N))
set(gcf,'position',[403   210   665   440])

