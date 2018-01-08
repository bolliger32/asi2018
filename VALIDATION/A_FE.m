%
% ANALISIS ESTAD�STICO EN CLIMATOLOG�A: 513428
% A. Montecinos
%
load F.dat
plot(F(:,1),F(:,2),'o')
%
% USAMOS ANOMALIAS
Af(:,1)=F(:,1)-mean(F(:,1));
Af(:,2)=F(:,2)-mean(F(:,2));

%
% anomalias estandarizadas
Asf(:,1)=(F(:,1)-mean(F(:,1)))/std(F(:,1));
Asf(:,2)=(F(:,2)-mean(F(:,2)))/std(F(:,2));

% 
% AZ = F E
% 
% Se quiere encontrar A tal que la varianza de A(:,1) sea la
% mayor posible...
%
VAR_F=std(F).*std(F)
VAR_total=sum(VAR_F)
%
% Encontramos la matriz T que transforma linealmente a F
% conviertiendola en Z
%
% alfa=45; 
alfa=55; 

a=alfa*pi/180;
E=[cos(a) -sin(a)
   sin(a)  cos(a)];

A=Af*E;
VAR_F=std(A).*std(A)
VAR_totalF=sum(VAR_F)

subplot(121)
plot(Af(:,1),Af(:,2),'o')
line([-6 8],[0 0])
line([0 0],[-6 8])
h=line([-6 6],[-6*tan(a) 6*tan(a)])
set(h,'LineStyle','--')
h=line([-6 6],[-6*tan(pi-a) 6*tan(pi-a)])
set(h,'LineStyle','--')
title('Diagrama de dispersion de anomalia de F')
subplot(122)
h=plot(A(:,1),A(:,2),'.')
set(h,'MarkerSize',20)
line([-8 8],[0 0])
line([0 0],[-8 8])
title('Diagrama de dispersion de anomalia de A')

%
% busquemos alfa tal que var(Z) sea maxima...
%
for alfa=1:89; 
 a=alfa*pi/180;
 E=[cos(a) -sin(a)
   sin(a)  cos(a)];
 A=Af*E;
 VARia(alfa,:)=std(A).*std(A);
end

plot(1:89,VARia(:,1))
grid





%
N=length(Af(:,1));
C=Af'*Af/(N-1);
[E,D]=eig(C);
L=diag(D);
L/sum(L)
A=Af*E;



