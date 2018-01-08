%
% Linear advection in 1D
%
% Resolution of dT/dt =  u0 dT/dx
% using an explicit forward time step 
%
%
clear all
close all
% xxxx
% GC between these line
%    paramters 
u0=5;			 % Advection speed [m.s-1]
dx=5.e3;		 % X resolution [m]
dt=800;	     % time step [s]
daymax=.5;		 % duration of simulation [days]
xmax=300.e3;    % Length of the basin [m]
% xxxxx
%
% Test of the CFL criterion
%
disp(['DT = ',num2str(dt)])
disp(['DX/DT = ',num2str(dx/dt)])
disp(['CMAX = ',num2str(u0)])
% xxx
% GC between these line
%    compute and print the advection CFL criteria
disp(['COURANT NUMBER for advection = ', num2str(u0*dt/dx)])
disp(' ')
% xxx
%
% Grid definition
%
% xxx
% GC between these line
%    define the 1D grid
%    then print the size of the grid
x=-100e3:dx:100e3;
L=length(x);
disp(['grid size : ',num2str(L)])
disp(' ')
% xxx
%
% Initial condition
% xxxxx
% GC: between these line
%     define the intial condition
% xxxxxx
T=15+5*exp(-(x/20e3).^2);
% First plot
%
%Plot the initial condition
% xxxxx
% GC : between this line 
%      plot the initial condition
% xxxxx
plot(x/1000,T,'r')
xlabel('X [km]')
ylabel('T [^oC]')
axis([-100 100 13 22]) 
hold on
% Main computing loop
%
nstep=0;
for time=0:dt/(24*3600):daymax
%
  nstep=nstep+1;
%
% Advection term
%
% xxxx
% GC : between this line define the right hand side term
%  1 - % first order "downstream" scheme - unconditionally unstable
%  d1x = zeros(size(T));
%  d1x(1:end-1) = T(2:end)-T(1:end-1);
%  d1x(end) = T(1) - T(end-1);
%  or 
%  2 - first order upstream scheme - unconditionally stable, diffusive
% IB: I think this is the same as downstream, you just have to adjust the
% different values in time stepping
%  rhs=-u0*(T(2:end)-T(1:end-1))/dx;
%  or
%  3 - % second order centered scheme - unconditionally unstable
 d1x = zeros(size(T));
 d1x(2:end-1) = T(3:end)-T(1:end-2);
 d1x(1) = T(2)-T(end);
 d1x(end) = T(1) - T(end-1);  % periodic at boundaries
%
% xxxxxxxxxx
 rhs = d1x * -u0/(2*dx);
 
% Time stepping
% Euler time stepping
%   T=T+dt*rhs;
% or
% Lax-Freidrich scheme
% technically, this scheme is space-and-time, using centered in space
% so for official Lax-Friedrich, we should use case 3 above
T_update_lw = zeros(size(T));
T_update_lw(2:end-1) = T(1:end-2)+T(3:end);
T_update_lw(1) = T(end)+T(2);
T_update_lw(end) = T(end-1)+T(1);
T = .5*T_update_lw + dt*rhs; 
% Boundary conditions
%
% xxx
% GC : between this line 
%      define the boundary conditions
% PERIODIC IS DEFAULT - INCLUDED IN SPATIAL SCHEME SO NONE NEEDED
% xxx
%
% Figure
%
% xxx
% GC : between this line 
%      plot the numerical solution 
%      at eache time step
  h1 = plot(x/1000,T,'m');
  Mov(nstep) = getframe;
% GC
%
% Analytical solution
%
  t=time*24*3600;
  h2=plot(x/1000,15+5*exp(-((x-u0*t)/20e3).^2),'r--');
  legend([h1 h2],{'Model','Analytical'})
%
  Mov(nstep)=getframe;
end
hold off
%
% Play the movie
%
% movie(Mov,2)
