%%

clear all

load SSTobs
load SSTsim

plot(tiempo,squeeze(dato(100,100,:))),datetick('x','yymm')
hold on
plot(tiempo,squeeze(dats(100,100,:)),'r'),datetick('x','yymm')

% analisis ciclo anual (promedio en el tiempo)
% sesgo
Oprom=mean(dato,3);
Yprom=mean(dats,3);
figure(1)
subplot(311)
contourf(lono,lato,Oprom'),colorbar
title('SST promedio observada')
caxis([13 21])
subplot(312)
contourf(lono,lato,Yprom'),colorbar
title('SST promedio simulada')
caxis([13 21])
subplot(313)
clabel(contourf(lono,lato,(Yprom-Oprom)',[-2:1:3])),colorbar
title('Sesgo')

% amplitud
Ostd=std(dato,0,3);
Ystd=std(dats,0,3);
figure(2)
subplot(311)
contourf(lono,lato,Ostd'),colorbar
title('SST desviacion estandar observada')
caxis([1 3])
subplot(312)
contourf(lono,lato,Ystd'),colorbar
title('SST desviacion estandar simulada')
caxis([1 3])
subplot(313)
clabel(contourf(lono,lato,(Ystd./Ostd)'),[.5 1 1.5]),colorbar
title('Amplitud')
%%
for j=1:length(lono)
    for i=1:length(lato);
        r(i,j)=corr(squeeze(dato(j,i,:)),squeeze(dats(j,i,:)));
    end
end
figure(3)
contourf(lono,lato,r),colorbar
title('correlacion (fase)')

%% rmse (raiz del error cuadratico medio)
for j=1:length(lono)
    for i=1:length(lato);
        y=squeeze(dats(j,i,:));
        o=squeeze(dato(j,i,:));
        rmse(i,j)=sqrt(sum((y-o).*(y-o))/(length(tiempo)-1));
    end
end
figure(4)
subplot(221)
contourf(lono,lato,rmse),colorbar
title('RMSE')
subplot(222)
contourf(lono,lato,r),colorbar
title('correlacion')
subplot(223)
contourf(lono,lato,(Yprom-Oprom)'),colorbar
title('sesgo')
subplot(224)
contourf(lono,lato,(Ystd./Ostd)'),colorbar
title('amplitud')

%% LECTURA SST observada y simulada

clear all

dias=[31 28 31 30 31 30 31 31 30 31 30 31];
c=0;
for i=2009:2010
    for j=1:12
        for k=1:dias(j)
            c=c+1;
            tiempo(c)=datenum(i,j,k);
            fecha(c,:)=[i j k];
        end
   end
end
tiempo=tiempo(1:3:end-10);
fecha=fecha(1:3:end-10,:);

% OBSERVA (en K)
lato=ncread('ext_2009_2010_All_SST-CCI.nc','lat',[151],[322],[1]);
lono=ncread('ext_2009_2010_All_SST-CCI.nc','lon',[51],[274],[1]);

dato=ncread('ext_2009_2010_All_SST-CCI.nc','sst',[51 151 1],[274 322 240],[1 1 3])-273.15;

% SIMULADO

lats=ncread('surf_roms_avg_All_Y2009_Y2010.nc','lat_rho');
lats=lats(1,:);
lons=ncread('surf_roms_avg_All_Y2009_Y2010.nc','lon_rho');
lons=lons(:,1);
for i=1:240
    aux=ncread('surf_roms_avg_All_Y2009_Y2010.nc','temp',[1 1 1 i],[410 561 1 1],[1 1 1 1]); % Celsius
    datt=squeeze(aux(:,:,1,:));
    xx=find(datt==0);
    datt(xx)=NaN;
    dats(:,:,i)=interp2(lats,lons,datt,lato',lono);
end

save SSTobs dato lato lono tiempo fecha
save SSTsim dats lato lono tiempo fecha

%%
