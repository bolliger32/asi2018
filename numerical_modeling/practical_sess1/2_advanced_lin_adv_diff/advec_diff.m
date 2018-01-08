function phi_num=advec_diff(schema_temp,schema_spat)
%%
%%
%% TESTS de differents schemas de discretisation temporelle
%% sur l equation d advection: d(phi)/dt+c*d(phi)/dx = 0
%%  
%% schema_temp : [1:7] choix du schema de discretisation temporelle
%%  1: Leap-Frog + Dissipation Euler
%%  2: Leap-Frog + Dissipation Leap-frog
%%  3: Adams-Bashforth-2  + Dissipation Euler
%%  4: Adams-Bashforth-3
%%  5: Runge-Kutta-4
%%  6: Leap-Frog + Euler + Dissipation avec Euler
%%  7: Leap-Frog+Asselin + Dissipation avec Euler
%%  8: Adams-Bashforth-Moulton-Predictor-Corrector
%%   
%% schema_spat : [1:2] choix du schema de discretisation spatiale  
%%  1: schema centre d ordre 2
%%  2: schema centre d ordre 4
  
nt=600;     % nombre d'iterations
dt=10.;     % pas de temps
c=0.05;     %vitesse d advection
nu=0.001;   %coefficient de diffusion
neuler=20; % Frequence a laquelle on intercale un Euler avant dans le LF

nx=101;
dx=10;
x=[0:dx:(nx-1)*dx];

time=[0:dt:(nt-1)*dt];

%
% Conditions initiales
%
x0=50*dx;
sig=10*dx;
phi0(1:nx)=exp(-(x-x0).^2./(2*sig^2));
phi0(1:nx)=0.;
phi0(fix(nx/2)-5:fix(nx/2)+5)=1.;
%phi0(1:nx)=cos(2*pi/((nx-1)*dx/8)*x);

%ATTENTION
%si condition aux limites periodiques, alors phi0(1)=phi0(nx)
phi0(1)=phi0(nx);

%
% INITIALISATION DES TABLEAUX
%
phi_num=zeros(nt,nx);
phib=[1:nx].*0;
phin=[1:nx].*0;
phia=[1:nx].*0;

dxphi=[1:nx].*0;
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
   case 1 % Leap-Frog + Dissipation avec Euler
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
    schema='Leap-Frog + Dissipation Euler';
    % -- Conditions Initiales
    phib=phi0;
    phin=phib-c*dt*d1x(phib,dx,schema_spat)+...
               nu*dt/(dx*dx)*dxx(phib,1);
    phi_num(1,:)=phib;
    phi_num(2,:)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=3:nt    
      phia=phib-2*dt*c*d1x(phin,dx,schema_spat)+...
               nu*2*dt/(dx*dx)*dxx(phib,1);
      phib=phin;
      phin=phia;
      phi_num(t,:)=phia;
    end
    elapsed_time = toc; % arret du chronometre    


%-------------------------------------------------------------------
%-------------------------------------------------------------------  
   case 2, % Leap-Frog + Dissipation avec Leap Frog
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
    schema='Leap-Frog + Dissipation Leap-Frog';
    % -- Conditions Initiales
    phib=phi0;
    phin=phib-c*dt*d1x(phib,dx,schema_spat)+...
               nu*dt/(dx*dx)*dxx(phib,1);
    phi_num(1,:)=phib;
    phi_num(2,:)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=3:nt    
      phia=phib-2*dt*c*d1x(phin,dx,schema_spat)+...
               nu*2*dt/(dx*dx)*dxx(phin,1);
      phib=phin;
      phin=phia;
      phi_num(t,:)=phia;
    end
    elapsed_time = toc; % arret du chronometre    

%-------------------------------------------------------------------
%-------------------------------------------------------------------
   case 3, % Adams-Bashforth-2 + Dissipation Euler
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
    schema='Adams-Basforths-2 + Dissipation Euler';
    % -- Conditions Initiales
    phib=phi0;
    phin=phib-c*dt*d1x(phib,dx,schema_spat)+...
               nu*dt/(dx*dx)*dxx(phib,1);
    phi_num(1,:)=phib;
    phi_num(2,:)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=3:nt    
      phia=phin-c*(dt/2)*(3*d1x(phin,dx,schema_spat)-d1x(phib,dx,schema_spat))+...
               nu*dt/(dx*dx)*dxx(phin,1);
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
    phib=phibb-c*dt*d1x(phibb,dx,schema_spat)+...
               nu*dt/(dx*dx)*dxx(phibb,1);
    phin=phibb-c*2*dt*d1x(phib, dx,schema_spat)+...
               nu*2*dt/(dx*dx)*dxx(phibb,1);
    phi_num(1,:)=phibb;
    phi_num(2,:)=phib;
    phi_num(3,:)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=4:nt    
      phia=phin-c*(dt/12)*(23*d1x(phin,dx,schema_spat)-16*d1x(phib,dx,schema_spat)+5*d1x(phibb,dx,schema_spat))+...
               nu*dt/(dx*dx)*dxx(phin,1);
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
    f1 = -c*dt*d1x(phin,dx,schema_spat);
    f2 = -c*dt*d1x(phin+f1/2,dx,schema_spat);
    f3 = -c*dt*d1x(phin+f2/2,dx,schema_spat);
    f4 = -c*dt*d1x(phin+f3,dx,schema_spat);    
    phia=phin+(f1+2*f2+2*f3+f4)/6+...
               nu*dt/(dx*dx)*dxx(phin,1);
    phin=phia;
    phi_num(t,:)=phia;
    end
    elapsed_time = toc; % arret du chronometre    

%-------------------------------------------------------------------
%-------------------------------------------------------------------  
   case 6, % Leap-Frog + Euler + Dissipation avec euler
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
    schema='Leap-Frog+Euler + Dissipation avec euler';
    % -- Conditions Initiales
    phib=phi0;
    phin=phib-c*dt*d1x(phib,dx,schema_spat)+...
         +dt*nu/(dx*dx)*dxx(phib,1);
    phi_num(1,:)=phib;
    phi_num(2,:)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=3:nt    
      if(mod(t,neuler)==0.)
          phia=phin-c*dt*d1x(phin,dx,schema_spat)+...
               nu*dt/(dx*dx)*dxx(phin,1);
      else
          phia=phib-2*dt*c*d1x(phin,dx,schema_spat)+...
               nu*2*dt/(dx*dx)*dxx(phib,1);
      end
      phib=phin;
      phin=phia;
      phi_num(t,:)=phia;
    end
    elapsed_time = toc; % arret du chronometre    

%-------------------------------------------------------------------
%-------------------------------------------------------------------  
   case 7, % Leap-Frog + Filtre Asselin + Dissipation avec Euler
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
    schema='Leap-Frog+Asselin + Dissipation avec Euler';
    gamma=0.1;
    % -- Conditions Initiales
    phib=phi0;
    phin=phib-c*dt*d1x(phib,dx,schema_spat)+dt*nu/(dx*dx)*dxx(phib,1);
    phi_num(1,:)=phib;
    phi_num(2,:)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=3:nt    
      phia=phib-2*dt*c*d1x(phin,dx,schema_spat)+...
               nu*2*dt/(dx*dx)*dxx(phib,1);
      phin=phin+gamma*(phib-2*phin+phia);
      phib=phin;
      phin=phia;
      phi_num(t,:)=phia;
    end
    elapsed_time = toc; % arret du chronometre    

%-------------------------------------------------------------------
%-------------------------------------------------------------------  
   case 8, % ABM Predictor-Corrector
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
    schema='Adams-Bashforth-Moulton-Predictor-Corrector';
    % -- Conditions Initiales
    phib=phi0;
    phin=phib-c*dt*d1x(phib,dx,schema_spat);
    phi_num(1,:)=phib;
    phi_num(2,:)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=3:nt   
      f1=phin-c*dt/2*(3*d1x(phin,dx,schema_spat)-d1x(phib,dx,schema_spat))
      phia=phin-c*dt/12*(5*d1x(f1,dx,schema_spat)...
                          +8*d1x(phin,dx,schema_spat)...
                          -d1x(phib,dx,schema_spat))...
                          +dt*nu/(dx*dx)*dxx(phin,1);
      phib=phin;
      phin=phia;
      phi_num(t,:)=phia;
    end
    elapsed_time = toc; % arret du chronometre    

    
    end 


%---------------------------------------------------------------------
%
%           Quelques plots
%
%---------------------------------------------------------------------
  
%figure;
%contourf(x,time,phi_num);
%hold on
%title(strcat(schema));  

  
disp '*************************'
disp '    TEMPS de CALCUL '
disp ''
disp 'elapsed_time'
disp '*************************'

