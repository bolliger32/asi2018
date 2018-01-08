function phi_num=chaleur(schema_temp,schema_spat)
%%
%%
%% TESTS de differents schemas de discretisation temporelle
%% sur l equation de la chaleur: d(phi)/dt=nu*d^2(phi)dx^2
%%  
%% schema_temp : [1:7] choix du schema de discretisation temporelle
%%  1: Euler Forward
%%  2: Leap-Frog
%%  3: Adams-Bashforth-2
%%  4: Adams-Bashforth-3
%%  5: Runge-Kutta-4
%%
%% schema_spat : [1:2] choix du schema de discretisation spatiale  
%%  1: schema centre d ordre 2
%%  2: schema centre d ordre 4
  
nt=300;      % nombre d'iterations
dt=0.2;      % pas de temps
nu=2;        % Coefficient de diffusion
neuler=20    % Frequence a laquelle on intercale un Euler avant dans le LF
nx=201;
dx=1;
x=[0:dx:(nx-1)*dx];
time=[0:dt:(nt-1)*dt];

%
% Conditions initiales
%
x0=50*dx;
sig=10*dx;
phi0(1:nx)=exp(-(x-x0).^2./(2*sig^2));

%
% INITIALISATION DES TABLEAUX
%
phi_num=zeros(nt,nx);
phib=[1:nx].*0;
phin=[1:nx].*0;
phia=[1:nx].*0;

dxxphi=[1:nx].*0;

%
% Variables
%
%    phi0 : profil initial
%
%    phib : phi BEFORE (pas de temps n-1)
%    phin : phi NOW    (pas de temps n  )
%    phia : phi AFTER  (pas de temps n+1)
%


tic % enclenchement du chronometre
  switch schema_temp
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
   case 1 %Euler Forward
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
    %-- Conditions Initiales
    schema='Euler Forward'
    phin=phi0;
    phi_num(1,:)=phin;
    %-- Boucle Temporelle
    tic % debut du chronometre
    for t=2:nt
      dxxphi=dxx(phin,schema_spat);      
      % Resolution des equations algebriques
      phia=phin+dt*nu/(dx*dx).*dxxphi;
      phin=phia;
      phi_num(t,:)=phia;
    end
    elapsed_time = toc; % arret du chronometre    
%-------------------------------------------------------------------
%-------------------------------------------------------------------  
   case 2, % Leap-Frog
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
    schema='Leap-Frog';
    % -- Conditions Initiales
    phib=phi0;
    dxxphi=dxx(phib,schema_spat);      
    phin=phib+2*dt*nu/(dx*dx).*dxxphi;
    phi_num(1,:)=phib;
    phi_num(2,:)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=3:nt    
      dxxphi=dxx(phin,schema_spat);      
      phia=phib+2*dt*nu/(dx*dx).*dxxphi;
      phib=phin;
      phin=phia;
      phi_num(t,:)=phia;
    end
    elapsed_time = toc; % arret du chronometre    

    %-------------------------------------------------------------------
%-------------------------------------------------------------------
   case 3, % Adams-Bashforth-2
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
    schema='Adams-Basforths-2';
    % -- Conditions Initiales
    phib=phi0;
    dxxphi=dxx(phib,schema_spat);      
    phin=phib+dt*nu/(dx*dx).*dxxphi;
    phi_num(1,:)=phib;
    phi_num(2,:)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=3:nt    
      dxxphi=dxx(phib,schema_spat);            
      phia=phin+(dt/2)*nu/(dx*dx)*(3*dxx(phin,schema_spat)-dxx(phib,schema_spat));
      phib=phin;
      phin=phia;
      phi_num(t,:)=phia;
    end
    elapsed_time = toc; % arret du chronometre    

    
%-------------------------------------------------------------------
%-------------------------------------------------------------------
   case 4, % Adams-Bashforth-3
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
    schema='Adams-Basforths-3';
    % -- Conditions Initiales
    phibb=phi0;
    phib=phibb+dt*nu/(dx*dx)*dxx(phibb,schema_spat);
    phin=phib+dt*nu/(dx*dx)*dxx(phib,schema_spat);
    phi_num(1,:)=phibb;
    phi_num(2,:)=phib;
    phi_num(3,:)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=4:nt    
      phia=phin+(dt/12)*nu/(dx*dx)*(23*dxx(phin,schema_spat)-16*dxx(phib,schema_spat)+5*dxx(phibb,schema_spat));
      phibb=phib;
      phib=phin;
      phin=phia;
      phi_num(t,:)=phia;
    end
    elapsed_time = toc; % arret du chronometre    

%-------------------------------------------------------------------
%-------------------------------------------------------------------    
   case 5, % RungeKutta 4eme ordre en temps
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
    schema='RungeKutta';
    % -- Conditions Initiales
    phin=phi0;
    phi_num(1,:)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=2:nt
    f1 = dt*nu/(dx*dx)*dxx(phin,schema_spat);
    f2 = dt*nu/(dx*dx)*dxx(phin+f1/2,schema_spat);
    f3 = dt*nu/(dx*dx)*dxx(phin+f2/2,schema_spat);
    f4 = dt*nu/(dx*dx)*dxx(phin+f3,schema_spat);    
    phia=phin+(f1+2*f2+2*f3+f4)/6;
    phin=phia;
    phi_num(t,:)=phia;
    end
    elapsed_time = toc; % arret du chronometre    

%-------------------------------------------------------------------
%-------------------------------------------------------------------    
   case 6, % RungeKutta 4eme ordre en temps pour laplacien en polaire
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
    schema='RungeKutta';
    % -- Conditions Initiales
    nt=3;    % nombre d'iterations
    dt=0.1;   % pas de temps
    nu=0.01;      % Coefficient de diffusion
    time=[0:dt:(nt-1)*dt];
    dr=0.01;
    r=[0.0001:dr:1];
    nr=length(r);
    rvortex=10*dr;
    phi_num=zeros(nt,nr);
    phib=[1:nr].*0;
    phin=[1:nr].*0;
    phia=[1:nr].*0;
    dxxphi=[1:nr].*0;


    phi0=exp(-(r/rvortex).^2);
    phin=phi0;
    phi_num(1,:)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=2:nt
    f1 = dt*nu/(dr*dr)*dxx(phin,2) + dt*nu*(1./r).*d1x(phin,dr,4);
    f2 = dt*nu/(dr*dr)*dxx(phin+f1/2,2) + dt*nu*(1./r).*d1x(phin+f1/2,dr,4);
    f3 = dt*nu/(dr*dr)*dxx(phin+f2/2,2) + dt*nu*(1./r).*d1x(phin+f2/2,dr,4);
    f4 = dt*nu/(dr*dr)*dxx(phin+f3,2) + dt*nu*(1./r).*d1x(phin+f3,dr,4);
    phia=phin+(f1+2*f2+2*f3+f4)/6;
    phin=phia;
    phi_num(t,:)=phia;
    end
    elapsed_time = toc; % arret du chronometre    
        
 %-------------------------------------------------------------------
%-------------------------------------------------------------------    
   case 7, % Euler Backward (Schema Implicite)
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
    schema='Euler Backward';
    % -- Conditions Initiales
    tic    
    phin=phi0;

        phi_num(1,:)=phin;
        e = ones(nx,1);
        A = spdiags([-e*nu*dt/(dx*dx) e+2*nu*dt/(dx*dx)*e -e*nu*dt/(dx*dx)], -1:1, nx, nx);
        for t=2:nt
        t
        %phia=transpose(inv(A)*phin');
        phia=transpose(A\phin');
        phin=phia;
        phi_num(t,:)=phia;
        end
        elapsed_time = toc; % arret du chronometre   
   

%-------------------------------------------------------------------
%-------------------------------------------------------------------    
% calcul de la solution exacte
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
%for it=1:length(time)
%  t=time(it);
%  phi_exact(it,:)=(rvortex^2)/(rvortex^2+4*nu*t)*exp(-(r.^2)/(rvortex^2+4*nu*t));
%end  
%   figure;
%   contourf(x,time,phi_num);
%   hold on
%   title(strcat(schema));  
end
  
figure(1)
hold on
plot(phi0,'r','Linewidth',6);
%plot(x,phi_num(t,:),'b','Linewidth',1);
title(strcat(schema,' time stepping'));  

disp('*************************')
disp('    TEMPS de CALCUL ')
disp('')
elapsed_time
disp('*************************')

