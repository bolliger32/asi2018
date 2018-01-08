function dxxphi=dxx(phin,schema_spat)
  dxxphi=zeros(size(phin));
  switch schema_spat 
   case 1
   %---------------------------------------------------------------
   %   Schema centre d ordre 2
   %---------------------------------------------------------------
      % Operateur derivee spatiale
      dxxphi(2:end-1)=phin(3:end)-2*phin(2:end-1)+phin(1:end-2);
      % Conditions aux limites (Plusieurs choix sont possibles)
      dxxphi(1)=2*dxxphi(2)-dxxphi(3);
      dxxphi(end)=2*dxxphi(end-1)-dxxphi(end-2);
      %Ici periodique
      dxxphi(1)=phin(2)-2*phin(1)+phin(end);
      dxxphi(end)=phin(1)-2*phin(end)+phin(end-1);
   
   %---------------------------------------------------------------
   %   Schema centre d ordre 4
   %---------------------------------------------------------------
   case 2
      % Operateur derivee spatiale
      dxxphi(3:end-2)=-1/12*phin(5:end)+4/3*phin(4:end-1)-15/6*phin(3:end-2)+4/3*phin(2:end-3)-1/12*phin(1:end-4);
      % Conditions aux limites (Plusieurs choix sont possibles)
      dxxphi(2)=2*dxxphi(3)-dxxphi(4);
      dxxphi(1)=2*dxxphi(2)-dxxphi(3);
      dxxphi(end-1)=2*dxxphi(end-2)-dxxphi(end-3);
      dxxphi(end)=2*dxxphi(end-1)-dxxphi(end-2);
      % Conditions symetrie polaire
      dxxphi(2)=-1/12*phin(4)+4/3*phin(3)-15/6*phin(2)+4/3*phin(1)-1/12*phin(2);
      dxxphi(1)=-1/12*phin(3)+4/3*phin(2)-15/6*phin(1)+4/3*phin(2)-1/12*phin(3);
      dxxphi(end-1)=2*dxxphi(end-2)-dxxphi(end-3);
      dxxphi(end)=2*dxxphi(end-1)-dxxphi(end-2);
      
      
      end
