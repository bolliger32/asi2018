%%

% clear all

load SSTobs
load SSTsim

% make double for EOF function to succeed
dats = double(dats);

% dimension lengths
nlat = length(lato);
nlon = length(lono);
nt = length(tiempo);

% create time X space matrix
dato_vec = reshape(dato,nlat*nlon,nt);
dats_vec = reshape(dats,nlat*nlon,nt);
no_nans = ~any([isnan(dato_vec(:,1)) isnan(dats_vec(:,1))],2);

Fo = dato_vec(no_nans,:)';
Fs = dats_vec(no_nans,:)';

% get anomalies
Fo = Fo - mean(Fo,1);
Fs = Fs - mean(Fs,1);

% run EOF function
No=length(Fo(:,1));            % tiempo
Mo=length(Fo(1,:));             % espacio
[Lo,Eo,Ao,erroro]=EOF(Fo,100);

Ns=length(Fs(:,1));            % tiempo
Ms=length(Fs(1,:));             % espacio
[Ls,Es,As,errors]=EOF(Fs,100);

%% Plotting
modo=3;

% A estandarizado
Ac_o=Ao(:,modo)/std(Ao(:,modo));

% E como correlacion
for i=1:Mo
    Ec_o(i)=corr(Ao(:,modo),Fo(:,i));
end

figure
subplot(231)
plot(tiempo,Ac_o),datetick('x','yymm')
title(['Componente principal modo Obs' num2str(modo) ', ' num2str(Lo(modo)/sum(Lo)*100) '%'])
subplot(234)
campo=reshape(NaN(length(lono),length(lato)),length(lono)*length(lato),1);
campo(no_nans)=Ec_o;
contourf(lono,lato,reshape(campo,length(lono),length(lato))'),colorbar
caxis([-1 0]);
title('Patron espacial como correlacion Obs')

% A estandarizado
Ac_s=As(:,modo)/std(As(:,modo));
% E como correlacion
for i=1:Ms
    Ec_s(i)=corr(As(:,modo),Fs(:,i));
end

subplot(232)
plot(tiempo,Ac_s),datetick('x','yymm')
title(['Componente principal modo Sim' num2str(modo) ', ' num2str(Ls(modo)/sum(Ls)*100) '%'])
subplot(235)
campo=reshape(NaN(length(lono),length(lato)),length(lono)*length(lato),1);
campo(no_nans)=Ec_s;
contourf(lono,lato,reshape(campo,length(lono),length(lato))'),colorbar
caxis([-1 0]);
title('Patron espacial como correlacion Sim')

% calc difference
Ac_diff = Ac_s - Ac_o;
Ec_diff = Ec_s - Ec_o;

subplot(233)
plot(tiempo,Ac_diff),datetick('x','yymm')
title(['Componente principal modo Sim-Obs'])
subplot(236)
campo=reshape(NaN(length(lono),length(lato)),length(lono)*length(lato),1);
campo(no_nans)=Ec_diff;
contourf(lono,lato,reshape(campo,length(lono),length(lato))'),colorbar
title('Patron espacial como correlacion Sim-Obs')

%% Remove seasonal trend
noseasonal_o = Fo' - Eo(:,1)*Ao(:,1)';
noseasonal_s = Fs' - Es(:,1)*As(:,1)';

Oprom = mean(noseasonal_o,1);
Yprom = mean(noseasonal_s,1);
orig_space_avg_o = nanmean(dato_vec,1);
figure
subplot(121)
plot(Oprom);
subplot(122)
plot(orig_space_avg_o);

% analisis ciclo anual (promedio en el tiempo)
% sesgo

figure(1)
subplot(223)
plot(tiempo,(Yprom-Oprom)),datetick('x','yymm')
title('Bias')

% amplitud
Ostd=nanstd(noseasonal_o,0,1);
Ystd=nanstd(noseasonal_s,0,1);
subplot(224)
plot(tiempo,(Ystd./Ostd)),datetick('x','yymm')
title('Amplitud')
%
r = zeros(nt,1);
for j=1:nt
    r(j)=corr(noseasonal_o(:,j),noseasonal_s(:,j));
end
subplot(222)
plot(tiempo,r),datetick('x','yymm')
title('correlacion (fase)')

% rmse (raiz del error cuadratico medio)
rmse = zeros(nt,1);
for j=1:nt
    y=squeeze(noseasonal_s(:,j));
    o=squeeze(noseasonal_o(:,j));
    rmse(j)=sqrt(sum((y-o).*(y-o))/(length(y)-1));
end
subplot(221)
plot(tiempo,rmse),datetick('x','yymm')
title('RMSE')
