%%

clear all

load SSTobs
load SSTsim

% analisis ciclo anual (promedio en el tiempo)
% sesgo
Oprom=squeeze(nanmean(nanmean(dato,1),2));
Yprom=squeeze(nanmean(nanmean(dats,1),2));
figure(1)
subplot(311)
plot(tiempo,Oprom),datetick('x','yymm')
title('SST promedio observada')
subplot(312)
plot(tiempo,Yprom),datetick('x','yymm')
title('SST promedio simulada')
subplot(313)
plot(tiempo,(Yprom-Oprom)),datetick('x','yymm')
title('Sesgo')

% amplitud
nlat = size(dato,1);
nlon = size(dato,2);
ntime = size(dato,3);
dato_vec = reshape(dato,nlat*nlon,ntime);
dats_vec = reshape(dats,nlat*nlon,ntime);
Ostd=nanstd(dato_vec,0,1);
Ystd=nanstd(dats_vec,0,1);
figure(2)
subplot(311)
plot(tiempo,Ostd),datetick('x','yymm')
title('SST desviacion estandar observada')
subplot(312)
plot(tiempo,Ystd),datetick('x','yymm')
title('SST desviacion estandar simulada')
subplot(313)
plot(tiempo,(Ystd./Ostd)),datetick('x','yymm')
title('Amplitud')
%%
r = zeros(ntime,1);
for j=1:ntime
    no_nans = ~any([isnan(dato_vec(:,j)) isnan(dats_vec(:,j))],2);
    dato_tmp = squeeze(dato_vec(no_nans,j));
    dats_tmp = squeeze(dats_vec(no_nans,j));
    r(j)=corr(dato_tmp,dats_tmp);
end
figure(3)
plot(tiempo,r),datetick('x','yymm')
title('correlacion (fase)')

%% rmse (raiz del error cuadratico medio)
rmse = zeros(ntime,1);
for j=1:ntime
    no_nans = ~any([isnan(dato_vec(:,j)) isnan(dats_vec(:,j))],2);
    y=squeeze(dats_vec(no_nans,j));
    o=squeeze(dato_vec(no_nans,j));
    rmse(j)=sqrt(sum((y-o).*(y-o))/(length(y)-1));
end
figure(4)
subplot(221)
plot(tiempo,rmse),datetick('x','yymm')
title('RMSE')
subplot(222)
plot(tiempo,r),datetick('x','yymm')
title('correlacion')
subplot(223)
plot(tiempo,(Yprom-Oprom)),datetick('x','yymm')
title('sesgo')
subplot(224)
plot(tiempo,(Ystd./Ostd)),datetick('x','yymm')
title('amplitud')

