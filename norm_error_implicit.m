clc 
clear
close all
ti=0; %time interval
tf=.1;
L=1;   %length of row
N=100;  %number of points on the rod


k=.01; %diffusibity constant
M=round(linspace(N^2*k*2*tf,5000,500));


x=linspace(0,L,N);
dx=x(2)-x(1);

ue=zeros(N,1);
error=zeros(N,numel(M)); %error for each version of M, over all points N
error_norm = zeros(1,numel(M)); %an entry for each choice of time points, M


%build the matrix A
A=zeros(N-2,N-2);

for i=1:numel(M)
    t=linspace(ti,tf,M(i));
    dt=t(2)-t(1);
    r=k*dt/dx^2;
    u=zeros(N,M(i)); %rows over time, columns over space
    u(:,1)=sin(pi*x)+.2*sin(10*pi*x); %initial conditions
    
    %exact solutions
    for p=1:N
        ue(p)=exp(-pi^2*k*t(M(i)))*sin(pi*x(p))+0.2*exp(-(10*pi)^2*k*t(M(i)))*sin(10*pi*x(p));
    end
    

    %boundary conditions 

    u(1,:) = 0;
    u(N,:)= 0;
    
    A(1,1)= 1+2*r;
    A(1,2)= -r;
    A(N-2,N-3)= -r;
    A(N-2,N-2)= 1+2*r;
    
    for q=2:N-3
        A(q,q-1)= -r;
        A(q,q)= 1+2*r;
        A(q,q+1)= -r;
    end
    
    b = u(2:N-1,1); %all of the interior points at time 1

    %for j=1:(.02*M(i)-1)
    
     %   K1 = (k/dx^2) * (u(1:N-2,j) - 2*u(2:N-1,j) + u(3:N,j)); %slope at point xi at time j
      %  Utilda(2:N-1) = u(2:N-1,j) + dt * K1; %approximation for value at point xi at time j+1
      %  K2 = (k/dx^2) * (Utilda(1:N-2) - 2*Utilda(2:N-1) + Utilda(3:N)); %slope at next time step
      %  u(2:N-1,j+1) = u(2:N-1,j) + dt*((K1+K2)/2); %improved euler's
        
    for j=1:(M(i)-1)
        b(1) = b(1) + r * u(1,j+1);  %the first and last interior points depend on the end points
        b(N-2) = b(N-2) + r * u(N,j+1);
    
        btilda = A\b;
        u(2:N-1,j+1) = btilda;
        b = btilda;
    end
    error(:,i)= (ue - u(:,M(i)));
    error_norm(i) = norm(error(:,i),inf); 
    
end
plot(tf./M*k/dx^2,error_norm,'*-')
set(gca,'FontSize',18)
xlabel('$r=\frac{k\Delta t}{\Delta x^2}','interpreter','latex','FontSize',24)
ylabel('Maximum error','interpreter','latex','FontSize',24)
title(sprintf('k = %0.3f',k),'interpreter','latex','fontsize',26)
set(gcf,'position',[403   210   665   440])