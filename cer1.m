%% Certamen 1 
clear variables; close all; clc
addpath(genpath('C:\Users\br1\OneDrive - Universidad de Concepción\Desktop\AEC\Certamen 1\CERTAMEN 1\UTILITIES'))
%% Cargando los datos
% Dimensiones: HE(Longitud, Latitud, Tiempo)
%vienen en sigle, por lo que tenemos que utilizar pasarlo a doble
lat = double(ncread('shum.nc','lat'));
lon = double(ncread('shum.nc','lon'));
t = double(ncread('shum.nc','time'));
indx_lon1=find(lon==150);
indx_lon2=find(lon==300);
[~,indx_lat1]=min(abs(lat - 20));
[~,indx_lat2]=min(abs(lat + 50));
%% extrayendo los datos
shum=double(ncread('shum.nc','shum',[indx_lon1 indx_lat1 1],[length(indx_lon1:indx_lon2) length(indx_lat1:indx_lat2) inf],[1 1 1]));
lat1=double(ncread('shum.nc','lat',37,74-36,1));
lon1=double(ncread('shum.nc','lon',81,161-80,1));
%% Promedio invierno verano
invd=[18 19 20];
verd=[24 25 26];
[X,Y,T] = size(shum); % Dimensiones
k=1; % contador
for i=12:12:T-6%para que calzar los años de inv y ver
    shum_inv(:,:,k)=mean(shum(:,:,invd+i-12),3);
    shum_ver(:,:,k)=mean(shum(:,:,verd+i-12),3);
    k=k+1;
end
clear k 
%% Promedio temporal y trimean
shum_ver_p=mean(shum_ver,3);%verano
shum_inv_p=mean(shum_inv,3);%invierno
% TRIMEAN
for y=1:1:38
    for x=1:1:81
        trimean_inv(x,y)=(prctile(squeeze(shum_inv(x,y,:)),25)+...
            2*prctile(squeeze(shum_inv(x,y,:)),50)+...
            prctile(squeeze(shum_inv(x,y,:)),75))/4;
        trimean_ver(x,y)=(prctile(squeeze(shum_ver(x,y,:)),25)+...
            2*prctile(squeeze(shum_ver(x,y,:)),50)+...
            prctile(squeeze(shum_ver(x,y,:)),75))/4;
    end
end

%% Grafico de los promedios invierno y verano
%verano trimediana
figure()
subplot(211)
m_proj('mercator','lon',[min(lon1) max(lon1)],'lat',[min(lat1) max(lat1)])
m_contourf(lon1,lat1,permute(trimean_ver,[2 1]))
colorbar
colormap('redbluecmap')
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Trimediana de la humedad específica para verano, región del Pacífico Tropical-Sur')
%verano mediana
subplot(212)
m_proj('mercator','lon',[min(lon1) max(lon1)],'lat',[min(lat1) max(lat1)])
m_contourf(lon1,lat1,permute(shum_ver_p,[2 1]))
colorbar
colormap('redbluecmap')
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Promedio de la humedad específica para vernano, región del Pacífico Tropical-Sur')
%Invierno trimediana y promedio
subplot(2,1,1)
m_proj('mercator','lon',[min(lon1) max(lon1)],'lat',[min(lat1) max(lat1)])
m_contourf(lon1,lat1,permute(trimean_inv,[2 1]))
colorbar
colormap('redbluecmap')
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Trimediana de la humedad específica para invierno, región del Pacífico Tropical-Sur')
%invierno
subplot(2,1,2)
m_proj('mercator','lon',[double(min(lon1)) double(max(lon1))],'lat',[double(min(lat1)) double(max(lat1))])
m_contourf(lon1,lat1,permute(shum_inv_p,[2 1]))
colorbar   
colormap(redbluecmap);
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Promedio de la humedad específica para invierno, región del Pacífico Tropical-Sur')
%diferencia invierno
figure()    
subplot(211)
m_proj('mercator','lon',[double(min(lon1)) double(max(lon1))],'lat',[double(min(lat1)) double(max(lat1))])
m_contourf(lon1,lat1,permute(shum_inv_p,[2 1])-permute(trimean_inv,[2 1]) )
colorbar 
caxis([-0.2 0.2]);
colormap(redbluecmap);
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',3);
hold off
title('Diferencia entre promedio y trimediana de la humedad específica para invierno, región del Pacífico Tropical-Sur')
subplot(212)
m_proj('mercator','lon',[double(min(lon1)) double(max(lon1))],'lat',[double(min(lat1)) double(max(lat1))])
m_contourf(lon1,lat1,permute(shum_ver_p,[2 1])-permute(trimean_ver,[2 1]) )
colorbar   
caxis([-0.2 0.2]);
colormap(redbluecmap);
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',3);
hold off
title('Diferencia entre promedio y trimediana de la humedad específica para verano, región del Pacífico Tropical-Sur')

%% Anomalias
ashum_inv = (shum_inv - mean(shum_inv,3));
ashum_ver = (shum_ver - mean(shum_ver,3));
%sin tendencias
sashum_inv=detrend(ashum_inv);
sashum_ver=detrend(ashum_ver);
%figuras
subplot(2,1,1)
m_proj('mercator','lon',[double(min(lon1)) double(max(lon1))],'lat',[double(min(lat1)) double(max(lat1))])
m_contourf(lon1,lat1,permute(ashum_inv(:,:,3),[2 1 3]))
colorbar
colormap('redbluecmap')
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Anomalia de humedad especifica invierno, región del Pacífico Tropical-Sur')
subplot(2,1,2)
m_proj('mercator','lon',[double(min(lon1)) double(max(lon1))],'lat',[double(min(lat1)) double(max(lat1))])
m_contourf(lon1,lat1,permute(ashum_ver(:,:,3),[2 1 3]))
colorbar
colormap('redbluecmap')
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Anomalia de humedad especifica verano, región del Pacífico Tropical-Sur')
%% Distribución de weibul
% De matiz a vector
shum_ver_vec=shum_ver(:);
shum_inv_vec=shum_inv(:);
%boxplot
boxplot([shum_inv_vec , shum_ver_vec],'Notch','on','Labels',{'Inverno','Verano'},'Whisker',1) 
title('Distribución de datos para la humedad especifica para invierno y verano')

a=0;
n=length(shum_inv_vec);
for i=1:n
wb(i)=(i-a)/(n+1-(2*a));
end
figure()
plot(sort(shum_ver_vec),wb,'r','LineWidth',2)
hold on
xlabel('Humedad específica [gr/kg]', 'FontName', 'times','FontSize', 14)
ylabel('Frecuencia','FontName', 'times','FontSize', 14)
title('Distribución de Weibull para la serie invierno y verano, región Pacífico Tropical-Sur','FontName', 'times', 'FontSize', 14, 'FontWeight' , 'bold')
  set(gca, 'Box', 'off','TickDir', 'out','TickLength', [.02 .02],'XMinorTick',...
    'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', [.3 .3 .3], 'YColor',...
    [.3 .3 .3],'LineWidth', 1)
plot(sort(shum_inv_vec),wb,'b','LineWidth',2)
legend('Serie de verano','Serie de invierno','Location','northwest','FontName', 'times');
hold off
grid minor
axis tight
% histograma
subplot(2,1,1)
 histogram(shum_inv)
 xlabel('Humedad específica [gr/kg]', 'FontName', 'times','FontSize', 14)
 ylabel('Frecuencia','FontName', 'times','FontSize', 14)
 title('Histograma para la humedad específica en invierno, región del Pacífico Tropical-Sur','FontName', 'times', 'FontSize', 14, 'FontWeight' , 'bold')
  set(gca, 'Box', 'off','TickDir', 'out','TickLength', [.02 .02],'XMinorTick',...
    'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', [.3 .3 .3], 'YColor',...
    [.3 .3 .3],'LineWidth', 1)
 grid minor
 subplot(2,1,2)
 histogram(shum_ver)
 xlabel('Humedad específica [gr/kg]', 'FontName', 'times','FontSize', 14)
 ylabel('Frecuencia','FontName', 'times','FontSize', 14)
 title('Histograma para la humedad específica en verano, región del Pacífico Tropical-Sur','FontName', 'times', 'FontSize', 14, 'FontWeight' , 'bold')
  set(gca, 'Box', 'off','TickDir', 'out','TickLength', [.02 .02],'XMinorTick',...
    'on', 'YMinorTick', 'on', 'YGrid', 'on','XColor', [.3 .3 .3], 'YColor',...
    [.3 .3 .3],'LineWidth', 1)
 grid minor
%% desviación standar vs rango intercuantil
shum_std_inv=std(shum_inv,0,3);
shum_std_ver=std(shum_ver,0,3);
% Rango intercuartil
for x = 1:81
    for y = 1:38
        rintercuartil   = prctile(squeeze(shum_inv(x,y,:)),75) - prctile(squeeze(shum_inv(x,y,:)),25);
        rintercuartil1   = prctile(squeeze(shum_ver(x,y,:)),75) - prctile(squeeze(shum_ver(x,y,:)),25);
        shum_inv_itq(x,y)     = rintercuartil;
        shum_ver_itq(x,y)     = rintercuartil1;
    end
end
%verano
subplot(211)
m_proj('mercator','lon',[double(min(lon1)) double(max(lon1))],'lat',[double(min(lat1)) double(max(lat1))])
m_contourf(lon1,lat1,permute(shum_ver_itq,[2 1]))
colorbar
caxis([-0.5 0.5])
colormap('redbluecmap') 
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Rango intercuartil de la humedad especifica para verano, región Pacífico Tropical-Sur',14,'FontWeight','bold','FontName', 'times')
subplot(212)
m_proj('mercator','lon',[double(min(lon1)) double(max(lon1))],'lat',[double(min(lat1)) double(max(lat1))])
m_contourf(lon1,lat1,permute(shum_inv_itq,[2 1]))
colorbar
caxis([-0.5 0.5])
colormap('redbluecmap')
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Rango intercuartil de la humedad especifica para invierno, región Pacífico Tropical-Sur',14,'FontWeight','bold','FontName', 'times')
%diferencias

%% Yulekendal vs skewness
figure()
subplot(211)
m_proj('mercator','lon',[double(min(lon1)) double(max(lon1))],'lat',[double(min(lat1)) double(max(lat1))])
m_contourf(lon1,lat1,permute(skewness(shum_inv,0,3),[2 1]))
colorbar
caxis([-1 1])
cmap = cmocean('balance');
colormap(cmap);
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Índice de simetría skewness de la humedad especifica para invierno, región Pacífico Tropical-Sur',14,'FontWeight','bold')
subplot(212)
m_proj('mercator','lon',[double(min(lon1)) double(max(lon1))],'lat',[double(min(lat1)) double(max(lat1))])
m_contourf(lon1,lat1,permute(skewness(shum_ver,0,3),[2 1]))
colorbar
caxis([-1 1])
cmap = cmocean('balance');
colormap(cmap);
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Índice de simetría skewness de la humedad especifica para verano, región Pacífico Tropical-Sur',14,'FontWeight','bold')

% Calculo del indice yule-kendall
[X,Y,~] = size(shum_inv)
for x = 1:X
    for y = 1:Y
        yule_kn_inv         = (prctile(squeeze(shum_inv(x,y,:)),25) - 2*prctile(squeeze(shum_inv(x,y,:)),50) + ...
            prctile(squeeze(shum_inv(x,y,:)),75))/shum_inv_itq(x,y);
        shum_inv_yuken(x,y)  = yule_kn_inv;
        yule_kn_ver= (prctile(squeeze(shum_ver(x,y,:)),25) - 2*prctile(squeeze(shum_ver(x,y,:)),50) + ...
            prctile(squeeze(shum_ver(x,y,:)),75))/shum_ver_itq(x,y);
        shum_ver_yuken(x,y) = yule_kn_ver;
    end
end
figure()
subplot(211)
m_proj('mercator','lon',[double(min(lon1)) double(max(lon1))],'lat',[double(min(lat1)) double(max(lat1))])
m_pcolor(lon1,lat1,permute(shum_inv_yuken,[2 1]))
shading Interp
colorbar
colormap(gca, redbluecmap)
caxis([-0.5 0.5])
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Índice de simetria Yule-Kendall para la humedad específica en invierno','FontSize',10,'FontWeight','bold')
subplot(212)
m_proj('mercator','lon',[double(min(lon1)) double(max(lon1))],'lat',[double(min(lat1)) double(max(lat1))])
m_pcolor(lon1,lat1,permute(shum_ver_yuken,[2 1]))
shading Interp
colorbar
colormap(gca, redbluecmap)
caxis([-0.5 0.5])
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Índice de simetria Yule-Kendall para la humedad específica en verano','FontSize',10,'FontWeight','bold')

%% Definiendo eventos para persistencia
% Calculo de la anomalia sin tendencia de shum
ashum_inv = detrend(shum_inv - mean(shum_inv,3));
ashum_ver = detrend(shum_ver - mean(shum_ver,3));
% Persistencia de Invierno
desf=5; % 1 año
crit=0; % sobre o bajo el promedio
%Definicion de eventos
%monte carlo
E01=ashum_inv(:,:,1+desf:end); % hoy
E02=ashum_inv(:,:,1:end-desf); %ayer   
for x=1:81
    for y=1:38
    E2=find(E02(x,y,:)>crit);
    E1E2=find(E01(x,y,E2)>crit);
    pers_inv(x,y)=length(E1E2)/length(E2);
    end
end
E01_r=ashum_ver(:,:,1+desf:end); % hoy
E02_r=ashum_ver(:,:,1:end-desf); %ayer   
for x=1:81
    for y=1:38
    E2_r=find(E02_r(x,y,:)>crit);
    E1E2_r=find(E01_r(x,y,E2_r)>crit);
    pers_ver(x,y)=length(E1E2_r)/length(E2_r);
    end
end

%visualizacion
subplot(2,1,1)
m_proj('mercator','lon',[double(min(lon1)) double(max(lon1))],'lat',[double(min(lat1)) double(max(lat1))])
m_pcolor(lon1,lat1,permute(pers_inv,[2 1]))
cmap = cmocean('amp',10);
colorbar
caxis([0 1]);
colormap(jet(10));
shading Interp
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Persistencia a 5 años de la humedad específica para invierno , región del Pacífico Tropical-Sur')

subplot(2,1,2)
m_proj('mercator','lon',[double(min(lon1)) double(max(lon1))],'lat',[double(min(lat1)) double(max(lat1))])
m_pcolor(lon1,lat1,pers_ver')
colorbar
caxis([0 1]);
colormap(jet(10));
shading Interp
hold on
m_gshhs_i('color','k','linewidth',2);
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
hold off
title('Persistencia a 5 años de la humedad específica para verano , región del Pacífico Tropical-Sur')
%monte carlo persistencia
%% SIGNIFICANCIA DE LA PERSISTENCIA
%% 
desf=1; % 1 año
crit=0; % sobre o bajo el promedio
%Definicion de eventos
E01=ashum_inv(:,:,1+desf:end); % hoy
E02=ashum_inv(:,:,1:end-desf); %ayer
%verano
E001=ashum_ver(:,:,1+desf:end); % hoy
E002=ashum_ver(:,:,1:end-desf); %ayer
%calculo de la persistencia
 for x = 1:size(ashum_inv,1)
     for y = 1:size(ashum_inv,2)
         for i = 1:1000
             E2=find(remuestreo(E02(x,y,:))>crit);
             E1E2=find(E01(x,y,E2)>crit);
             re_pers_inv(x,y,i)=length(E1E2)/length(E2);
             %verano
             E22=find(remuestreo(E002(x,y,:))>crit);
             E11E22=find(E001(x,y,E22)>crit);
             re_pers_ver(x,y,i)=length(E11E22)/length(E22);
         end
     end
 end
 
 %% intervalo de confianza
for x = 1:size(ashum_inv,1)
    for y = 1:size(ashum_inv,2)
        ashum_pers_inv_ic_sup(x,y) = prctile(squeeze(re_pers_inv(x,y,:)),97.5);
        ashum_pers_inv_ic_inf(x,y) = prctile(squeeze(re_pers_inv(x,y,:)),2.5);
        %verano
        ashum_pers_ver_ic_sup(x,y) = prctile(squeeze(re_pers_ver(x,y,:)),97.5);
        ashum_pers_ver_ic_inf(x,y) = prctile(squeeze(re_pers_ver(x,y,:)),2.5);
    end
end
% Aplicando la signficancia:
ashum_pers_sig_inv = pers_inv;
ashum_pers_sig_ver = pers_ver;
for x = 1:size(ashum_inv,1)
    for y = 1:size(ashum_inv,2)
        if pers_inv(x,y) < ashum_pers_inv_ic_sup(x,y) && pers_inv(x,y) > ashum_pers_inv_ic_inf(x,y)
            ashum_pers_sig_inv(x,y) = NaN;
        else
            ashum_per_sig_inv(x,y) =1;
        end
        if pers_ver(x,y) < ashum_pers_ver_ic_sup(x,y) && pers_ver(x,y) > ashum_pers_ver_ic_inf(x,y)
            ashum_pers_sig_ver(x,y) = NaN;
        else
            ashum_per_sig_ver(x,y) =1;
        end
    end
end
subplot(211)
m_proj('mercator','lon',[min(lon1(:)) max(lon1(:))],'lat',[min(lat1(:)) ...
    max(lat1(:))]);
m_pcolor(lon1,lat1,ashum_pers_sig_inv')
shading Interp
colorbar(gca,'Location','SouthOutside')
colormap(jet(10))
caxis([0 1]);
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Significancia de persistencia a 5 años de la humedad específica para verano , región del Pacífico Tropical-Sur')
subplot(212)
m_proj('mercator','lon',[min(lon1(:)) max(lon1(:))],'lat',[min(lat1(:)) ...
    max(lat1(:))]);
m_pcolor(lon1,lat1,ashum_pers_sig_ver')
shading Interp
colorbar(gca,'Location','SouthOutside')
colormap(jet(10))
caxis([0 1]);
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Significancia de persistencia a 5 años de la humedad específica para verano , región del Pacífico Tropical-Sur')



