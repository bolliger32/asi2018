function phi_num=advection(temp_scheme,spat_scheme)
%%
%% Author: Steven Herbette (UBO-Brest- France)
%%
%% We TEST several numerical combinations of finite difference numerical
%% schemes on the linear advection equation:
%% d(phi)/dt + c*d(phi)/dx = 0
%%
%% phi is the scalar field we are interested in.
%%  ex: phi can be potential temperature
%%
%% There are two partial derivatives that we must approximate thanks to
%% finite difference schemes:
%%   - a time derivative
%%   - a spatial derivative
%%
%% You can select various schemes for the time derivative by giving values
%% from 1 to 7 to the variable temp_scheme:
%% tempschema_temp :
%%  1: Euler Forward
%%  2: Leap-Frog
%%  3: Adams-Bashforth-2
%%  4: Adams-Bashforth-3
%%  5: Runge-Kutta-4
%%  6: Leap-Frog + Euler
%%  7: Leap-frog + Asselin
%%  8: Adams-Bashforth-Moulton-Predictor-Corrector
%%
%% The partial derivative is solved in the separate function:
%% d1x.m
%% The function admits three arguments:
%% dxphi=d1x(phi,dx,spat_scheme)
%% phi is an array of values of the scalar field.
%% It is discrete in space
%% dx is the grid space
%% spat_scheme selects which schemes you want to use to approximate the
%% spatial derivative
%%  1: upwind (non centered - order 1) (if c>0)
%%  2: downwind (non centered - order 1) (if c>0)
%%  3: centered - order 2)
%%  4: downwind order 3 (non centered - order 3)
%%  5: quick (centrered - order 4)


%% A) First things to do:
%% 1) define a grid in space
nx=200;              % number of points in x
dx=10;               % grid space
x=[0:dx:(nx-1)*dx]; % discrete values of x

%% 2) Disctretize time
nt=600;    % number of total iterations
dt=1.6;   % time step
time=[0:dt:(nt-1)*dt]; %discrete values of time

c=5;      % advection speed

%% 3) Initial conditions
% phi0 : initial profile
% Gaussian
x0=nx/2*dx;
sig0=10*dx;
phi0(1:nx)=exp(-(x-x0).^2./(2*sig0^2));

% Square
%phi0(1:nx)=0;
%phi0(fix(nx/2)-5:fix(nx/2)+5)=1.;

% Gaussian + Square
%phi0(50-5:50)=2/sig0*[0:5].*dx;
%phi0(50:55)=2/sig0*[5:-1:0].*dx;

% Cosine
%phi0(1:nx)=cos(2*pi/((nx-1)*dx/8)*x);
%phi0(1:nx)=cos(0.8*pi/dx*x);

% WARNING
% Periodic boundary conditions : phi0(1)=phi0(nx)
phi0(1)=phi0(nx);

%
% 4) INITIALISATION of arrays
%
% Variables
%

%
%    phib : phi BEFORE (time n-1)
%    phin : phi NOW    (time n  )
%    phia : phi AFTER  (time n+1)
%
%
phi_num=zeros(nt,nx);
phib=[1:nx].*0;
phin=[1:nx].*0;
phia=[1:nx].*0;

dxphi=[1:nx].*0;
dxxphi=[1:nx].*0;


%%
%% B ) Algorithm
%% Integrate in time the finite difference equations
%%

switch temp_scheme
    %-------------------------------------------------------------------
    %-------------------------------------------------------------------
    case 1 %Euler Forward
        %-------------------------------------------------------------------
        %-------------------------------------------------------------------
        %-- First steps - Start
        schema='Euler Forward'
        phin=phi0;
        phi_num(1,:)=phin;
        %-- Time Loop
        tic % start counting time
        for t=2:nt
            phia=phin-c*dt*d1x(phin,dx,spat_scheme);
            phin=phia;
            phi_num(t,:)=phia;
        end
        elapsed_time = toc; % stop counting time
        %-------------------------------------------------------------------
        %-------------------------------------------------------------------
    case 2, % Leap-Frog
        %-------------------------------------------------------------------
        %-------------------------------------------------------------------
        schema='Leap-Frog';
        % -- First steps - Start
        phib=phi0;
        phin=phib-c*dt*d1x(phib,dx,spat_scheme);
        phi_num(1,:)=phib;
        phi_num(2,:)=phin;
        % -- Time Loop
        tic % Start Time counting
        for t=3:nt
            phia=phib-2*dt*c*d1x(phin,dx,spat_scheme);
            phib=phin;
            phin=phia;
            phi_num(t,:)=phia;
        end
        elapsed_time = toc; % Stop Time counting
        
        %-------------------------------------------------------------------
        %-------------------------------------------------------------------
    case 3, % Adams-Bashforth-2
        %-------------------------------------------------------------------
        %-------------------------------------------------------------------
        schema='Adams-Basforths-2';
        % -- First steps - Start
        phib=phi0;
        phin=phib-c*dt*d1x(phib,dx,spat_scheme);
        phi_num(1,:)=phib;
        phi_num(2,:)=phin;
        % -- Time Loop
        tic % Start Time counting
        for t=3:nt
            phia=phin-c*(dt/2)*(3*d1x(phin,dx,spat_scheme)-d1x(phib,dx,spat_scheme));
            phib=phin;
            phin=phia;
            phi_num(t,:)=phia;
        end
        elapsed_time = toc; % Stop Time counting
        
        
        %-------------------------------------------------------------------
        %-------------------------------------------------------------------
    case 4, % Adams-Bashforth-3
        %-------------------------------------------------------------------
        %-------------------------------------------------------------------
        schema='Adams-Basforths-3';
        % -- First steps - Start
        phibb=phi0;
        phib=phibb-c*dt*d1x(phibb,dx,spat_scheme);
        phin=phib -c*dt*d1x(phib, dx,spat_scheme);
        phi_num(1,:)=phibb;
        phi_num(2,:)=phib;
        phi_num(3,:)=phin;
        % -- Time Loop
        tic % Start Time counting
        for t=4:nt
            phia=phin-c*(dt/12)*(23*d1x(phin,dx,spat_scheme)-16*d1x(phib,dx,spat_scheme)+5*d1x(phibb,dx,spat_scheme));
            phibb=phib;
            phib=phin;
            phin=phia;
            phi_num(t,:)=phia;
        end
        elapsed_time = toc; % Stop Time counting
        
        %-------------------------------------------------------------------
        %-------------------------------------------------------------------
    case 5, % RungeKutta 4eme ordre en temps
        %-------------------------------------------------------------------
        %-------------------------------------------------------------------
        schema='RungeKutta';
        % -- First steps - Start
        phin=phi0;
        phi_num(1,:)=phin;
        % -- Time Loop
        tic % Start Time counting
        for t=2:nt
            f1 = -c*dt*d1x(phin,dx,spat_scheme);
            f2 = -c*dt*d1x(phin+f1/2,dx,spat_scheme);
            f3 = -c*dt*d1x(phin+f2/2,dx,spat_scheme);
            f4 = -c*dt*d1x(phin+f3,dx,spat_scheme);
            phia=phin+(f1+2*f2+2*f3+f4)/6;
            phin=phia;
            phi_num(t,:)=phia;
        end
        elapsed_time = toc; % Stop Time counting
        
        %-------------------------------------------------------------------x
        %-------------------------------------------------------------------
    case 6, % Leap-Frog + Euler
        %-------------------------------------------------------------------
        %-------------------------------------------------------------------
        schema='Leap-Frog+Euler';
        % -- First steps - Start
        phib=phi0;
        phin=phib-c*dt*d1x(phib,dx,spat_scheme);
        phi_num(1,:)=phib;
        phi_num(2,:)=phin;
        
        neuler=50; % Frequence a laquelle on intercale un Euler avant dans le LF
        % -- Time Loop
        tic % Start Time counting
        for t=3:nt
            if(mod(t,neuler)==0.)
                phia=phin-c*dt*d1x(phin,dx,spat_scheme);
            else
                phia=phib-2*dt*c*d1x(phin,dx,spat_scheme);
            end
            phib=phin;
            phin=phia;
            phi_num(t,:)=phia;
        end
        elapsed_time = toc; % Stop Time counting
        
        %-------------------------------------------------------------------
        %-------------------------------------------------------------------
    case 7, % Leap-Frog + Filtre Asselin
        %-------------------------------------------------------------------
        %-------------------------------------------------------------------
        schema='Leap-Frog+Asselin';
        gamma=0.2;
        % -- First steps - Start
        phib=phi0;
        phin=phib-c*dt*d1x(phib,dx,spat_scheme);
        phi_num(1,:)=phib;
        phi_num(2,:)=phin;
        % -- Time Loop
        tic % Start Time counting
        for t=3:nt
            phia=phib-2*dt*c*d1x(phin,dx,spat_scheme);
            phin=phin+gamma*(phib-2*phin+phia);
            phib=phin;
            phin=phia;
            phi_num(t,:)=phia;
        end
        elapsed_time = toc; % Stop Time counting
        
        
        %-------------------------------------------------------------------
        %-------------------------------------------------------------------
    case 8, % Adams-Moulton-Bashforth-Predictor-Corrector
        %-------------------------------------------------------------------
        %-------------------------------------------------------------------
        schema='Adams-Bashforth-Moulton-Predictor-Corrector';
        % -- First steps - Start
        phib=phi0;
        phin=phib-c*dt*d1x(phib,dx,spat_scheme);
        phi_num(1,:)=phib;
        phi_num(2,:)=phin;
        % -- Time Loop
        tic % Start Time counting
        for t=3:nt
            f1=phin-c*dt/2*(3*d1x(phin,dx,spat_scheme)-d1x(phib,dx,spat_scheme));
            phia=phin-c*dt/12*(5*d1x(f1,dx,spat_scheme)...
                +8*d1x(phin,dx,spat_scheme)...
                -d1x(phib,dx,spat_scheme));
            phib=phin;
            phin=phia;
            phi_num(t,:)=phia;
        end
        elapsed_time = toc; % Stop Time counting

end
%---------------------------------------------------------------------
%
%           Quelques plots
%
%--------------------------------------------------------------------- 
figure(1)
hold on
plot(phi0,'r','Linewidth',4);
%plot(x,phi_num(t,:),'b','Linewidth',1);
title(strcat(schema)); 

disp('*************************')
disp('    TEMPS de CALCUL ')
disp('')
elapsed_time
disp('*************************')

