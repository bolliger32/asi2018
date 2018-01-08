%------------------------------------------------------------------
%-----------------------------------------------------------------
function dxphi=d1x(phi,dx,schema_spat)
%
% Calcule une derivee spatiale avec differents schemas de
% discretisation temporelle
%
%  phi: Champ 1D
%  dx: Pas d espace ou "Pas de grille" (spatiale)
%  schema_spat: numero du schema de discretisation spatiale
%      1: schema decentre a gauche (upstream si c>0)
%      2: schema decentre a droite (upstream si c>0)
%      3: schema centre d ordre 2
%      4: schema centre d ordre 4
%
%
% COURS MA2 2006-2007
% Steven Herbette
%------------------------------------------------------------------
%-----------------------------------------------------------------

dxphi=zeros(size(phi));
switch schema_spat
    %---------------------------------------------------------------
    %   Schema decentre a gauche
    %---------------------------------------------------------------
    case 1
        % Operateur derivee spatiale
        dxphi(2:end)=phi(2:end)-phi(1:end-1);
        % Conditions aux limites (Plusieurs choix sont possibles)
        % Extrapolation ("Slip")
        dxphi(1)=2*dxphi(2)-dxphi(3);
        % Conditions aux limites periodiques
        % remarque: Si phi0(1)=phi0(Nx), alors les CL ecrites
        %           ci-dessous doivent garrantir la periodicite
        dxphi(1)=phi(1)-phi(end);
        dxphi=dxphi/dx;
        
        %---------------------------------------------------------------
        %   Schema decentre a droite
        %---------------------------------------------------------------
    case 2
        % Operateur derivee spatiale
        dxphi(1:end-1)=phi(2:end)-phi(1:end-1);
        % Conditions aux limites (Plusieurs choix sont possibles)
        % Extrapolation ("Slip")
        dxphi(end)=2*dxphi(end-1)-dxphi(end-2);
        %Condition aux limites periodiques
        % remarque: Si phi0(1)=phi0(Nx), alors les CL ecrites
        %           ci-dessous doivent garrantir la periodicite
        dxphi(end)=phi(1)-phi(end);
        dxphi=dxphi/dx;
        %---------------------------------------------------------------
        %   Schema centre d ordre 2
        %---------------------------------------------------------------
    case 3
        % Operateur derivee spatiale
        dxphi(2:end-1)=phi(3:end)-phi(1:end-2);
        % Conditions aux limites (Plusieurs choix sont possibles)
        % Extrapolation ("Slip")
        %dxphi(1)=2*dxphi(2)-dxphi(3);
        %dxphi(end)=2*dxphi(end-1)-dxphi(end-2);
        %Condition aux limites periodiques
        % remarque: Si phi0(1)=phi0(Nx), alors les CL ecrites
        %           ci-dessous doivent garrantir la periodicite
        dxphi(1)=phi(2)-phi(end);
        dxphi(end)=phi(1)-phi(end-1);
        dxphi=dxphi/(2*dx);
        
    case 4
        dxphi(3:end-1)=2*phi(4:end)+3*phi(3:end-1)-6*phi(2:end-2)+phi(1:end-3);
        dxphi(2)=2*phi(3)+3*phi(2)-6*phi(1)+phi(end);
        dxphi(1)=2*phi(2)+3*phi(1)-6*phi(end)+phi(end-1);
        dxphi(end)=2*phi(1)+3*phi(end)-6*phi(end-1)+phi(end-2);
        dxphi=dxphi./(6*dx);
        
        %---------------------------------------------------------------
        %   Schema centre d ordre 4
        %---------------------------------------------------------------
    case 5
        % Operateur derivee spatiale
        dxphi(3:end-2)=-1/12*phi(5:end)+2/3*phi(4:end-1)-2/3*phi(2:end-3)+1/12*phi(1:end-4);
        % Conditions aux limites (Plusieurs choix sont possibles)
        dxphi(2)=2*dxphi(3)-dxphi(4);
        dxphi(1)=2*dxphi(2)-dxphi(3);
        dxphi(end-1)=2*dxphi(end-2)-dxphi(end-3);
        dxphi(end)=2*dxphi(end-1)-dxphi(end-2);
        %Condition aux limites periodiques
        % remarque: Si phi0(1)=phi0(Nx), alors les CL ecrites
        %           ci-dessous doivent garrantir la periodicite
        
        dxphi(2)=-1/12*phi(4)+2/3*phi(3)-2/3*phi(1)+1/12*phi(end);
        dxphi(end-1)=-1/12*phi(1)+2/3*phi(end)-2/3*phi(end-2)+1/12*phi(end-3);
        dxphi(1)=-1/12*phi(3)+2/3*phi(2)-2/3*phi(end)+1/12*phi(end-1);
        dxphi(end)=-1/12*phi(2)+2/3*phi(1)-2/3*phi(end-1)+1/12*phi(end-2);
        
        
        %dxphi(2)=-1/12*phi(4)+2/3*phi(3)-2/3*phi(1)+1/12*phi(2);
        %dxphi(1)=0;
        %dxphi(end-1)=2*dxphi(end-2)-dxphi(end-3);
        %dxphi(end)=2*dxphi(end-1)-dxphi(end-2);
        
        
        
        dxphi=dxphi/dx;
end
