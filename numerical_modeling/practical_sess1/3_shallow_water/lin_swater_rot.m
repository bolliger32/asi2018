%
% Linear shallow water equations : constant rotation
%
%  Du/Dt =fv-gDzeta/Dx
%  Dv/Dt =-fu-gDzeta/Dy
%  Dzeta/Dt = -H(Du/Dx-Dv/Dy)
%
clear all
close all
%
dx=100.e3;		% X resolution [m]
dy=dx;		        % Y resolution [m]
dt=216;			% time step [s]
daymax=3;		% duration of simulation [days]
xmax=20000.e3;		% Length of the basin [m]
ymax=20000.e3;		% Width of the basin [m]
Hmax=5000.;		% depth [m]
g=9.81;			% gravity acceleration [m.s-2]
rho0=1025;		% Density [kg.m-3]
f0=-1e-4;		% coriolis parameter [s-1] 
                        % (negative in the southern hemisphere)
%
R0=1000e3;		% Diameter of the initial perturbation [m]
zeta0=0.2;		% Amplitude initial perturbation [m]
%
daydiag=0.1;		% Number of days between diagnostics
dayplot=0.1;		% Number of days between images
skip=4;			% Plot 1 arrow every 'skip' grid points
%
rsponge=1e-4;		% Maximum value of rayleigh damping in the sponge [m.s-1]
xsponge=1000e3;		% Size of the sponge [m]
%
% Test of the CFL criterion for barotropic gravity waves
%
c=sqrt(g*Hmax)*dt/dx;
disp(['DT = ',num2str(dt)])
disp(['DX/DT = ',num2str(dx/dt)])
disp(['CMAX = ',num2str(sqrt(g*Hmax))])
disp(['COURANT NUMBER = ',num2str(c)])
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
% Get the Coriolis parameter
%
f= f0 + 0*y_rho;
f_u=0.5*(f(:,1:end-1)+f(:,2:end));
f_v=0.5*(f(1:end-1,:)+f(2:end,:));
%
% Compute the sponge
%
rspg=0*x_rho;
if rsponge>0
  distE=-x_rho+xmax/2;
  distW=x_rho+xmax/2;
  distN=-y_rho+ymax/2;
  distS=y_rho+ymax/2;
  dist=cat(3,distE,distW,distN,distS);
  dist=min(dist,[],3);
  dist(dist<0)=0;
  rspg=0.5*(cos(pi*dist/xsponge)+1);
  rspg(dist>=xsponge)=0;
  rspg=rsponge*rspg;
end
rspg=rspg./H;
rspg_u=0.5*(rspg(:,2:end)+rspg(:,1:end-1));
rspg_v=0.5*(rspg(2:end,:)+rspg(1:end-1,:));
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
% Compute u-rhs terms 
%
  ru=0* y_u;
%
% Coriolis term
%
  ru(2:end-1,:)=f_u(2:end-1,:).*...
                  0.25.*(v(1:end-1,1:end-1)+v(1:end-1,2:end)+...
                         v(2:end  ,1:end-1)+v(2:end  ,2:end));
%
% Pressure gradient term
%
  ru=ru-g.*(zeta(:,2:end)-zeta(:,1:end-1))./dx;
%
% Compute v-rhs terms 
%
  rv=0* y_v;
%
% Coriolis term
%
  rv(:,2:end-1)=-f_v(:,2:end-1).*...
                   0.25.*(u(1:end-1,1:end-1)+u(1:end-1,2:end)+...
                          u(2:end  ,1:end-1)+u(2:end,2:end));
%
% Pressure gradient term
%
  rv=rv-g.*(zeta(2:end,:)-zeta(1:end-1,:))./dy;
%
% Step U : forward step U(t+1)=U(t)+dt*F(t)
%
  u(2:end-1,2:end-1)=u(2:end-1,2:end-1)+dt*ru(2:end-1,2:end-1);
%
% Open boundary conditions for u
%
  u(:,1)=(u(:,1)+c.*u(:,2))./(1+c);
  u(:,end)=(u(:,end)+c.*u(:,end-1))./(1+c);
  u(1,:)=(u(1,:)+c.*u(2,:))./(1+c);
  u(end,:)=(u(end,:)+c.*u(end-1,:))./(1+c);
%
% Step V : forward step V(t+1)=V(t)+dt*F(t)
%
  v(2:end-1,2:end-1)=v(2:end-1,2:end-1)+dt*rv(2:end-1,2:end-1);
%
% Open boundary conditions for v
%
  v(:,1)=(v(:,1)+c.*v(:,2))./(1+c);
  v(:,end)=(v(:,end)+c.*v(:,end-1))./(1+c);
  v(1,:)=(v(1,:)+c.*v(2,:))./(1+c);
  v(end,:)=(v(end,:)+c.*v(end-1,:))./(1+c);
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
  zeta(2:end-1,2:end-1)=zeta(2:end-1,2:end-1)-dt.*H(2:end-1,2:end-1).*rz(2:end-1,2:end-1);
%
% Open boundary conditions for zeta
%
  zeta(:,1)=(zeta(:,1)+c.*zeta(:,2))./(1+c);
  zeta(:,end)=(zeta(:,end)+c.*zeta(:,end-1))./(1+c);
  zeta(1,:)=(zeta(1,:)+c.*zeta(2,:))./(1+c);
  zeta(end,:)=(zeta(end,:)+c.*zeta(end-1,:))./(1+c);
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
