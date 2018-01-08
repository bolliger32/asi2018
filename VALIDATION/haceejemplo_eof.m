%% ersst 1950 a 2016

clear all

lat=ncread('sst.mnmean.v4.nc','lat');
lon=ncread('sst.mnmean.v4.nc','lon');
aux=ncread('sst.mnmean.v4.nc','sst',[1 1 1153],[180 89 804],[1 1 1]);

anhos=1950:2016;

for i=1:length(anhos)
    sst(:,:,i)=mean(aux(:,:,1:12),3);
    aux(:,:,1:12)=[];
end

xx=find(sst<-1000);
sst(xx)=NaN;
figure(1),contourf(lon,lat,squeeze(sst(:,:,1))'),colorbar

puntero=find(isnan(reshape(squeeze(sst(:,:,1)),length(lon)*length(lat),1))==0);
for i=1:length(anhos)
    aux=reshape(squeeze(sst(:,:,i)),length(lon)*length(lat),1)';
    tsm(i,:)=aux(puntero);
end

campo=reshape(NaN(length(lon),length(lat)),length(lon)*length(lat),1);
campo(puntero)=tsm(1,:);
figure(2),contourf(lon,lat,reshape(campo,length(lon),length(lat))'),colorbar

% Usamos anomalias promedio anual, definimos matriz F(N,M)

mtsm=mean(tsm);
for i=1:length(anhos)
    F(i,:)=tsm(i,:)-mtsm;
end

N=length(F(:,1));            % tiempo
M=length(F(1,:));             % espacio
[L,E,A,error]=EOF(F);

L(1:10)/sum(L)

% en caso de ir lento...
%[L,A,E,error]=EOF(F);
%L(1:10)/sum(L)


%%
modo=3;
% A estandarizado
Ac=A(:,modo)/std(A(:,modo));
% E como correlacion
for i=1:M
    Ec(i)=corr(A(:,modo),F(:,i));
end

figure
subplot(211)
plot(anhos,Ac)
title(['Componente principal modo ' num2str(modo) ', ' num2str(L(modo)/sum(L)*100) '%'])
subplot(212)
campo=reshape(NaN(length(lon),length(lat)),length(lon)*length(lat),1);
campo(puntero)=Ec;
contourf(lon,lat,reshape(campo,length(lon),length(lat))'),colorbar
title('Patron espacial como correlacion')

