%
% Linear advection in 1D
%
%  dT/dt =  -K d2T/dx2
%
%
clear all
close all
%
nu=20000;		% Diffusion parameter [m2.s-1]
dx=5.e3;		% X resolution [m]
dt=1250;		% time step [s]
daymax=10;		% duration of simulation [days]
xmax=100.e3;    % Length of the basin [m]
%
% Test of the CFL criterion
%
disp(['DT = ',num2str(dt)])
disp(['DX/DT = ',num2str(dx/dt)])
disp(['COURANT NUMBER (2 nu dt2/dx2) = ',num2str(2*nu*dt/dx^2)])
disp(' ')
%
% Grid definition
%
x=-100e3:dx:100e3;
L=length(x);
disp(['grid size : ',num2str(L)])
disp(' ')
%
% Initial condition
%
T=15+5*exp(-(x/20e3).^2);
%
% First plot
%
plot(x/1000,T,'r')
xlabel('X [km]')
ylabel('T [^oC]')
axis([-100 100 13 22]) 
hold on
%
% Main computing loop
%
nstep=0;
for time=0:dt/(24*3600):daymax
  nstep=nstep+1;
%
% Diffusion term
% 
  rhs=nu*(T(1:end-2)-2*T(2:end-1)+T(3:end))./(dx^2);
%
% Euler time stepping
%
  T(2:end-1)=T(2:end-1)+dt*rhs;
%
% Boundary conditions (Dirichlet)
%
  T(1)=15;
  T(end)=15;
%
% Figure
%
  plot(x/1000,T,'m')
  Mov(nstep) = getframe;
end
hold off
movie(Mov,2)
