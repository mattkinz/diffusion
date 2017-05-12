clc
clear

L=10; %rod
N=100; %point on rod
T=1; %time interval
M=2000; %points in time
D=15; %diffusivity constant

x=linspace(0,L,N);
dx = x(2) - x(1);
t=linspace(0,T,M);
dt=t(2)-t(1);

r=D*dt/dx^2;
u=zeros(N,M); %rows through space, columns through time

%boundary conditions
u(1,:)=sin((2*pi/T)*t);
u(N,:)=sin((2*pi/T)*t); 
%u(1,:)=0;
%u(N,:)=0;

%build the matrix A
A=zeros(N-2,N-2);
A(1,1)= 1+2*r;
A(1,2)= -r;
A(N-2,N-3)= -r;
A(N-2,N-2)= 1+2*r;

for k=2:N-3
    A(k,k-1)= -r;
    A(k,k)= 1+2*r;
    A(k,k+1)= -r;
end

%initial conditions
u(:,1)=sin((2*pi/L)*x);

b = u(2:N-1,1);

for k=2:M-1
    b(1) = b(1) + r * u(1,k+1);
    b(N-2) = b(N-2) + r * u(N,k+1);
    
    btilda = A\b;
    u(2:N-1,k) = btilda;
    b = btilda;
end

[X,T]=meshgrid(x,t);

pcolor(X,T,u'); shading interp; colorbar;
colormap jet
xlabel('rod')
ylabel('time')
    


