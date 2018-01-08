function phi_num=oscillation(numero);
%%
%%
%% TESTS de differents schemas de discretisation temporelle
%% sur l equation d ocillation dt(phi)=i*omega*phi
%%  
%% numero : [1:7] choix du schema de discretisation temporelle
%%  1: Euler Forward
%%  2: Leap-Frog
%%  3: Adams-Bashforth-2
%%  4: Adams-Bashforth-3
%%  5: Runge-Kutta-4
%%  6: Leap-Frog+Euler ts les neuler pas de temps
%%  7: Leap-Frog+Filtre d Asselin
  
nt=200;       % nombre d'iterations
dt=5.;        % pas de temps
T=100;        % Periode theorique de l oscillation
omega=2*pi/T; % pulsation
neuler=5; % frequence a laquelle on intercale un Euler avant dans le LF

time=[0:dt:(nt-1)*dt];
%
% Conditions initiales
%
phi0=1;

%
% Tableau de sortie
phi_num=zeros(nt);
phi_num(1)=phi0;

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
  switch numero
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
   case 1 %Euler Forward
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
    %-- Conditions Initiales
    schema='Euler Forward'
    phin=phi0;
    phi_num(1)=phin;
    %-- Boucle Temporelle
    tic % debut du chronometre
    for t=2:nt
      phia=phin+dt*i*omega*phin;
      phin=phia;
      phi_num(t)=phia;
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
    %phin=phib+dt*i*omega*phib+0.1;
    phin=phib+dt*i*omega*phib;
    phi_num(1)=phib;
    phi_num(2)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=3:nt    
      phia=phib+2*dt*i*omega*phin;
      phib=phin;
      phin=phia;
      phi_num(t)=phia;
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
    phin=phib+dt*i*omega*phib;
    phi_num(1)=phib;
    phi_num(2)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=3:nt    
      phia=phin+dt/2*(3*i*omega*phin-i*omega*phib);
      phib=phin;
      phin=phia;
      phi_num(t)=phia;
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
    phib=phibb+dt*i*omega*phibb;
    phin=phibb+2*dt*i*omega*phib;
    phi_num(1)=phibb;
    phi_num(2)=phib;
    phi_num(3)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=4:nt    
      phia=phin+dt/12*(23*i*omega*phin-16*i*omega*phib+5*i*omega*phibb);
      phibb=phib;
      phib=phin;
      phin=phia;
      phi_num(t)=phia;
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
    phi_num(1)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=2:nt
    f1 = i*omega*dt*(phin);
    f2 = i*omega*dt*(phin+f1/2);
    f3 = i*omega*dt*(phin+f2/2);
    f4 = i*omega*dt*(phin+f3);    
    phia=phin+(f1+2*f2+2*f3+f4)/6;
    phin=phia;
    phi_num(t)=phia;
    end
    elapsed_time = toc; % arret du chronometre    

%-------------------------------------------------------------------
%-------------------------------------------------------------------    
   case 6, % Leap-Frog+Euler avant tous les neuler pas de temps
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
    schema='Leap-Frog+Euler';
    % -- Conditions Initiales
    phib=phi0;
    phin=phib+dt*i*omega*phib;
    phi_num(1)=phib;
    phi_num(2)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=3:nt
      if(mod(t,neuler)==0.)
        phia=phin+dt*i*omega*phin;
      else
        phia=phib+2*dt*i*omega*phin;
      end
      phib=phin;
      phin=phia;
      phi_num(t)=phia;
    end
    elapsed_time = toc; % arret du chronometre    
    %nt
%-------------------------------------------------------------------
%-------------------------------------------------------------------        
   case 7, % Leap-Frog + Asselin filter
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
     schema='Leap-Frog + Asselin';
     gamma=0.1;
     %omega=0; % pour faire apparaitre le mode numerique
    % -- Conditions Initiales
    phib=phi0;
    %phin=phib+0.2;
    phin=phib+dt*i*omega*phib;
    phi_num(1)=phib;
    phi_num(2)=phin;
    % -- Boucle Temporelle
    tic % debut du chronometre
    for t=3:nt
      phia=phib+2*dt*i*omega*phin;
      phin=phin+gamma*(phib-2*phin+phia);
      phib=phin;
      phin=phia;
      phi_num(t)=phia;
    end
    elapsed_time = toc; % arret du chronometre    
     
  end 

%-------------------------------------------------------------------
%-------------------------------------------------------------------    
% calcul de la solution exacte
%-------------------------------------------------------------------
%-------------------------------------------------------------------    
  phi_exact=phi0*exp(i*omega*time);
  
  real_phi_num=real(phi_num);
  real_phi_exact=phi_exact;
  
  figure;
  aa1=plot(time,real_phi_num,'b');
  hold on
  aa2=plot(time,real_phi_exact,'r');  
  %legend([aa2],{'Numerical approximation'})
  legend([aa2],{'Exact solution'})
  size(real_phi_num)
  size(real_phi_exact)
  title(['Temporal scheme is : ',schema]);  

disp('*************************')
disp('    TEMPS de CALCUL ')
disp('')
elapsed_time
disp('*************************')
