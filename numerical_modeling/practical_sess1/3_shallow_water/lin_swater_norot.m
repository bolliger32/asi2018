%
% Linear shallow water equations : no rotation + walls
%
%  Du/Dt = -gDzeta/Dx
%  Dv/Dt = -gDzeta/Dy
%  Dzeta/Dt = -H(Du/Dx-Dv/Dy)
%
clear all
close all
%
% Parameters
%
dx=100.e3;		% X resolution [m]
dy=dx;		    % Y resolution [m]
dt=432;			% time step [s]
daymax=3;		% duration of simulation [days]
xmax=20000.e3;		% Length of the basin [m]
ymax=20000.e3;		% Width of the basin [m]
Hmax=2500.;		% depth [m]
g=9.81;			% gravity acceleration [m.s-2]
rho0=1025;		% Density [kg.m-3]
%
R0=1000e3;		% Diameter of the initial perturbation [m]
zeta0=0.2;		% Amplitude initial perturbation [m]
%
daydiag=0.1;		% Number of days between diagnostics
dayplot=0.1;		% Number of days between images
skip=0;			% Plot 1 arrow every 'skip' grid points
%
% Test of the CFL criterion for barotropic gravity waves
%
disp(['DT = ',num2str(dt)])
disp(['DX/DT = ',num2str(dx/dt)])
disp(['CMAX = ',num2str(sqrt(g*Hmax))])
disp(['CURRENT NUMBER = ',num2str(sqrt(g*Hmax)*dt/dx)])
disp(' ')
%
% Build a horizontal C-grid
%
x_rho=-xmax/2-dx/2:dx:xmax/2+dx/2;
y_rho=-ymax/2-dy/2:dy:ymax/2+dy/2;
x_u=0.5*(x_rho(2:end)+x_rho(1:end-1));
y_u=y_rho;
x_v=x_rho;
y_v=0.5*(y_rho(2:end)+y_rho(1:end-1));
[x_rho,y_rho]=meshgrid(x_rho,y_rho);
[x_u,y_u]=meshgrid(x_u,y_u);
[x_v,y_v]=meshgrid(x_v,y_v);
[M,L]=size(x_rho);
disp(['grid size : ',num2str(L),' X ',num2str(M)])
%
% Create the bottom topography
%
H=Hmax + 0*x_rho;
%
% Define the initial conditions
%
zeta=zeta0*exp(-(x_rho.^2+y_rho.^2)/R0^2);
u=0*y_u;
v=0*y_v;
%
nstep=0;
nframe=0;
%
% Main computing loop
%
for time=0:dt/(24*3600):daymax
%
  nstep=nstep+1;
%
% Compute u-rhs terms : pressure gradient
%
  ru=-g.*(zeta(:,2:end)-zeta(:,1:end-1))./dx;
%
% Compute v-rhs terms : pressure gradient
%
  rv=-g.*(zeta(2:end,:)-zeta(1:end-1,:))./dy;
%
% Step U : forward step U(t+1)=U(t)+dt*F(t)
%
  u(2:end-1,2:end-1)=u(2:end-1,2:end-1)+dt*ru(2:end-1,2:end-1);
%
% Boundary conditions for u: wall
%
  u(:,1)=0;
  u(:,end)=0;
  u(1,:)=u(2,:);
  u(end,:)=u(end-1,:);
%
% Step V : forward step V(t+1)=V(t)+dt*F(t)
%
  v(2:end-1,2:end-1)=v(2:end-1,2:end-1)+dt*rv(2:end-1,2:end-1);
%
% Boundary conditions for v: wall
%
  v(:,1)=v(:,2);
  v(:,end)=v(:,end-1);
  v(1,:)=0;
  v(end,:)=0;
%
% Compute zeta-rhs terms 
%
  rz=0.*y_rho;
%
% Divergence of u 
%
  rz(:,2:end-1)=rz(:,2:end-1)+(u(:,2:end)-u(:,1:end-1))/dx;
%
% Divergence of v 
%
  rz(2:end-1,:)=rz(2:end-1,:)+(v(2:end,:)-v(1:end-1,:))/dy;
%
% Backward step for zeta (don't forget H)
%
  zeta(2:end-1,2:end-1)=zeta(2:end-1,2:end-1)-...
                        dt.*H(2:end-1,2:end-1).*rz(2:end-1,2:end-1);
%
% Boundary conditions for zeta (no gradient across a wall)
%
  zeta(:,1)=zeta(:,2);
  zeta(:,end)=zeta(:,end-1);
  zeta(1,:)=zeta(2,:);
  zeta(end,:)=zeta(end-1,:);
%
% Diagnostics if necessary
%
  if rem(time,daydiag)==0
    disp(['Day: ',num2str(time),' - Step: ',num2str(nstep)])
    mz=100*mean(mean(zeta));
    mu=mean(mean(u));  
    mv=mean(mean(v));  
    disp(['Mean zeta: ',num2str(mz),...
          ' - u: ',num2str(mu),...
          ' - v: ',num2str(mv)])
    mz=100*max(max(zeta));
    mu=max(max(u));  
    mv=max(max(v));  
    disp(['Max zeta: ',num2str(mz),...
          ' - u: ',num2str(mu),...
          ' - v: ',num2str(mv)])
    mz=100*min(min(zeta));
    mu=min(min(u));  
    mv=min(min(v));  
    disp(['Min zeta: ',num2str(mz),...
          ' - u: ',num2str(mu),...
          ' - v: ',num2str(mv)])
    disp(' ')
  end
%
% Plot if necessary
%
  if rem(time,dayplot)==0 
    nframe=nframe+1;
    pcolor(x_rho/1000,y_rho/1000,100*zeta)
    shading interp
    axis image
    caxis([-10 10])
    colorbar('horiz')
    if skip>0
      hold on
      ur=0.5*(u(2:end-1,2:end)+u(2:end-1,1:end-1));
      vr=0.5*(v(2:end,2:end-1)+v(1:end-1,2:end-1));
      quiver(x_rho(2:skip:end-1,2:skip:end-1)/1000,...
             y_rho(2:skip:end-1,2:skip:end-1)/1000,...
             ur(1:skip:end,1:skip:end),vr(1:skip:end,1:skip:end),'k')
      hold off
    end
    title(['zeta [cm] - day: ',num2str(time)])
    xlabel('X [km]')
    ylabel('Y [km]')
    warning off
    Mov(nframe) = getframe;
  end
end
movie(Mov,5)
