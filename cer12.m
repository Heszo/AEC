%%parte 2
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
[X,Y,T] = size(shum_inv); % Dimensiones

%% Cargando datos para ENSSO
%parte desde 1871 para que parta del 1949 son 78x12, se borran los datos.
a=load('N34_oni_1871_2020.txt');
b=a(:,6);
b=b(a(:,1)>=1949 & a(:,1)<=2020);
aux=reshape(b,12,2020-1949+1);
n34_inv=mean(aux(6:8,:))'; % junio-agosto : 1948-2020
k=1;
JJA = [6 7 8];
DEF=[12 13 14];

for i=12:12:864
    if i<864
        n34_ver(k)=mean(b(DEF-12+i));
    else
        n34_ver(k)=b(end);
    end
    k=k+1;
end
n34_ver=n34_ver';
%% Calculo de la anomalia estandarizada de shum
ashum_inv = (shum_inv - mean(shum_inv,3))./std(shum_inv,0,3);
ashum_ver = (shum_ver - mean(shum_ver,3))./std(shum_ver,0,3);
%% Limpieza
clear T k i b 

for x = 1:size(shum_inv,1)
    for y = 1:size(shum_inv,2)
        ashum_inv_cor(x,y) = corr(squeeze(ashum_inv(x,y,:)),n34_inv);
    end
end

for x = 1:size(shum_ver,1)
    for y = 1:size(shum_ver,2)
        ashum_ver_cor(x,y) = corr(squeeze(ashum_ver(x,y,:)),n34_ver);
    end
end
figure()
subplot(2,1,1)
m_proj('mercator','lon',[min(lon1) max(lon1)],'lat',[min(lat1) max(lat1)])
m_contourf(lon1,lat1,permute(ashum_inv_cor(:,:),[2 1]))
colorbar
colormap(cmocean('balance',20))
caxis([-1 1])
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Mapa de correlación de humedad especifica y modo climatico ENSO para invierno, región del Pacífico Tropical-Sur')

subplot(2,1,2)
m_proj('mercator','lon',[min(lon1) max(lon1)],'lat',[min(lat1) max(lat1)])
m_contourf(lon1,lat1,permute(ashum_ver_cor(:,:),[2 1]))
colorbar
colormap(cmocean('balance',20))
caxis([-1 1])
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Mapa de correlación de humedad especifica y modo climatico ENSO para verano, región del Pacífico Tropical-Sur')


%% test de significancia monte carlo
 for x = 1:size(ashum_inv,1)
    for y = 1:size(ashum_inv,2)
        for i = 1:500
            ashum_corsig_inv(x,y,i) = corr(remuestreo(squeeze(ashum_inv(x,y,:)))',n34_inv);
            ashum_corsig_ver(x,y,i) = corr(remuestreo(squeeze(ashum_ver(x,y,:)))',n34_ver);
            
        end
    end
 end
% intervalo de confianza
for x = 1:size(ashum_inv,1)
    for y = 1:size(ashum_inv,2)
        ashum_cor2_ic_sup(x,y) = prctile(squeeze(ashum_corsig_inv(x,y,:)),97.5);
        ashum_cor2_ic_inf(x,y) = prctile(squeeze(ashum_corsig_inv(x,y,:)),2.5);
        ashum_cor2_ic_supv(x,y)= prctile(squeeze(ashum_corsig_inv(x,y,:)),97.5);
        ashum_cor2_ic_infv(x,y)= prctile(squeeze(ashum_corsig_inv(x,y,:)),2.5);
    end
end
% Aplicando la signficancia:
% Para temperatura superficial del mar:
for x = 1:size(ashum_inv,1)
    for y = 1:size(ashum_inv,2)
        if ashum_inv_cor(x,y) < ashum_cor2_ic_sup(x,y) && ashum_inv_cor(x,y) > ashum_cor2_ic_inf(x,y)
            ashum_inv_cor(x,y) = NaN;
        else
        end
    end
end
for x = 1:size(ashum_ver,1)
    for y = 1:size(ashum_ver,2)
        if ashum_ver_cor(x,y) < ashum_cor2_ic_supv(x,y) && ashum_ver_cor(x,y) > ashum_cor2_ic_infv(x,y)
            ashum_ver_cor(x,y) = NaN;
        else
        end
    end
end
%plots correlación significancia

subplot(2,2,1)
m_proj('mercator','lon',[min(lon1) max(lon1)],'lat',[min(lat1) max(lat1)])
m_contourf(lon1,lat1,permute(ashum_inv_cor(:,:),[2 1]))
colorbar
colormap(cmocean('balance',20))
caxis([-1 1])
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Mapa de correlación de humedad especifica y modo climatico ENSO para invierno, región del Pacífico Tropical-Sur')

subplot(2,2,2)
m_proj('mercator','lon',[min(lon1) max(lon1)],'lat',[min(lat1) max(lat1)])
m_contourf(lon1,lat1,permute(ashum_ver_cor(:,:),[2 1]))
colorbar
colormap(cmocean('balance',20))
caxis([-1 1])
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Mapa de correlación de humedad especifica y modo climatico ENSO para verano, región del Pacífico Tropical-Sur')

subplot(2,1,1)
m_proj('mercator','lon',[min(lon1) max(lon1)],'lat',[min(lat1) max(lat1)])
m_contourf(lon1,lat1,permute(ashum_inv_cor,[2 1]))
shading Interp
colorbar
colormap(cmocean('balance',20))
caxis([-1 1])
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Significancia para correlación en invierno')

subplot(2,1,2)
m_proj('mercator','lon',[min(lon1) max(lon1)],'lat',[min(lat1) max(lat1)])
m_contourf(lon1,lat1,permute(ashum_ver_cor,[2 1]))
shading Interp
colorbar
colormap(cmocean('balance',20))
caxis([-1 1])
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Significancia para correlación en verano')
% Probabilidad de que la humedad específica sea sobre el promedio
% debido a que ocurra el niñp en invierno
%% EL NIÑO
EN_inv=find(n34_inv>=0.5);
EN_ver=find(n34_ver>=0.5);

for x=1:X
for y=1:Y
e1=find(squeeze(shum_inv(x,y,EN_inv))>prctile(squeeze(shum_inv(x,y,:))',75));
proba(x,y)=length(e1)/length(EN_inv);
%verano
e2=find(squeeze(shum_ver(x,y,EN_ver))>prctile(squeeze(shum_ver(x,y,:))',75));
proba_ver(x,y)=length(e2)/length(EN_ver);
end
end
subplot(212)
m_proj('mercator','lon',[min(lon1(:)) max(lon1(:))],'lat',[min(lat1(:)) ...
    max(lat1(:))]);
m_pcolor(lon1,lat1,proba_ver')
shading Interp
colorbar
colormap(jet(10))
caxis([0 1]);
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Probabilidad de que la humedad específica esté sobre el percentil 75 debido debido a que ocurra el niño en verano')
subplot(211)
m_proj('mercator','lon',[min(lon1(:)) max(lon1(:))],'lat',[min(lat1(:)) ...
    max(lat1(:))]);
m_pcolor(lon1,lat1,proba')
shading Interp
colorbar
colormap(jet(10))
caxis([0 1]);
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Probabilidad de que la humedad específica esté sobre el percentil 75 debido debido a que ocurra el niño en invierno')
clear e1 e2
%Monte carlo
EN_inv=find(n34_inv>=0.5);
EN_ver=find(n34_ver>=0.5);
for x=1:X
for y=1:Y
for i=1:500
aux1=remuestreo(shum_ver(x,y,:));
aux2=remuestreo(shum_inv(x,y,:));    
ee1=find(aux2(EN_inv)>prctile(squeeze(shum_inv(x,y,:))',75));
probaa_rem(x,y,i)=length(ee1)/length(EN_inv);
%verano
ee2=find(aux1(EN_ver)>prctile(squeeze(shum_ver(x,y,:))',75));
probaa_rem_ver(x,y,i)=length(ee2)/length(EN_ver);
end
end
end

% intervalo de confianza
for x = 1:size(shum_inv,1)
for y = 1:size(shum_inv,2)
proba_rem_sup(x,y) = prctile(squeeze(probaa_rem(x,y,:)),97.5);
proba_rem_inf(x,y) = prctile(squeeze(probaa_rem(x,y,:)),2.5);
%verano
proba_rem_sup_ver(x,y) = prctile(squeeze(probaa_rem_ver(x,y,:)),97.5);
proba_rem_inf_ver(x,y) = prctile(squeeze(probaa_rem_ver(x,y,:)),2.5);
end
end
% Aplicando la signficancia:
for x = 1:size(shum_inv,1)
    for y = 1:size(shum_inv,2)
    if proba(x,y) < proba_rem_sup(x,y) && proba(x,y) > proba_rem_inf(x,y)
    proba(x,y) = NaN;
    else
    end
    if proba_ver(x,y) < proba_rem_sup_ver(x,y) && proba_ver(x,y) > proba_rem_inf_ver(x,y)
    proba_ver(x,y) = NaN;
    else
    end
end
end
subplot(211)
m_proj('mercator','lon',[min(lon1(:)) max(lon1(:))],'lat',[min(lat1(:)) ...
    max(lat1(:))]);
m_pcolor(lon1,lat1,proba')
shading Interp
colorbar
colormap(jet(5))
caxis([0 1])
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Significancia campo compuesto invierno para el niño')

subplot(212)
m_proj('mercator','lon',[min(lon1(:)) max(lon1(:))],'lat',[min(lat1(:)) ...
    max(lat1(:))]);
m_pcolor(lon1,lat1,proba_ver')
shading Interp
colorbar
caxis([0 1])
colormap(jet(10))
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Significancia campo compuesto verano para el niño') 
 
clear EN_inv EN_ver e1 e2 proba proba_ver proba_rem proba_rem_ver 
clear proba_rem_sup proba_rem_inf proba_rem_sup_ver proba_rem_inf_ver proba_sig proba_sig_ver

%% La niña
E12_inv_ni=find(n34_inv<=-0.5);
E12_ver_ni=find(n34_ver<=-0.5);
for x=1:X
    for y=1:Y
    e1_ni= find(shum_inv(x,y,E12_inv_ni)<prctile(squeeze(shum_inv(x,y,:)),25));
    proba_inv_ni(x,y)=length(e1_ni)/length(E12_inv_ni);
  
    e2_ni= find(shum_ver(x,y,E12_ver_ni)<prctile(squeeze(shum_ver(x,y,:)),25));
    proba_ver_ni(x,y)=length(e2_ni)/length(E12_ver_ni);
    end
end
figure()
subplot(211)
m_proj('mercator','lon',[min(lon1(:)) max(lon1(:))],'lat',[min(lat1(:)) ...
    max(lat1(:))]);
m_pcolor(lon1,lat1,proba_inv_ni')
shading Interp
colorbar
colormap(jet(10))
caxis([0 1]);
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Probabilidad de que la humedad específica esté bajo el percentil 25 debido debido a que ocurra la niña en invierno')
subplot(212)
m_proj('mercator','lon',[min(lon1(:)) max(lon1(:))],'lat',[min(lat1(:)) ...
    max(lat1(:))]);
m_pcolor(lon1,lat1,proba_ver_ni')
shading Interp
colorbar
colormap(jet(10))
caxis([0 1]);
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Probabilidad de que la humedad específica esté bajo el percentil 25 debido debido a que ocurra la niña en verano')
% Probabilidad de que la humedad específica sea sobre el promedio
% debido a que ocurra la niña en invierno
% Carlos montes
clear aux1 aux2 probaa_rem probaa_rem_ver ee1 ee2
EN_inv_ni=find(n34_inv<=-0.5);
EN_ver_ni=find(n34_ver<=-0.5);
for x=1:X
for y=1:Y
for i=1:100
aux1=remuestreo(shum_ver(x,y,:));
aux2=remuestreo(shum_inv(x,y,:));    
ee1=find(aux2(EN_inv_ni)<prctile(squeeze(shum_inv(x,y,:))',25));
probaa_rem(x,y,i)=length(ee1)/length(EN_inv_ni);
%verano
ee2=find(aux1(EN_ver_ni)<prctile(squeeze(shum_ver(x,y,:))',25));
probaa_rem_ver(x,y,i)=length(ee2)/length(EN_ver_ni);
end
end
end
clear proba_rem_sup proba_rem_inf proba_rem_sup_ver proba_rem_inf_ver
% intervalo de confianza
for x = 1:size(shum_inv,1)
for y = 1:size(shum_inv,2)
proba_rem_sup(x,y) = prctile(squeeze(probaa_rem(x,y,:)),97.5);
proba_rem_inf(x,y) = prctile(squeeze(probaa_rem(x,y,:)),2.5);
%verano
proba_rem_sup_ver(x,y) = prctile(squeeze(probaa_rem_ver(x,y,:)),97.5);
proba_rem_inf_ver(x,y) = prctile(squeeze(probaa_rem_ver(x,y,:)),2.5);
end
end
% Aplicando la signficancia:
clear proba_sig proba_sig_ver
for x = 1:size(shum_inv,1)
    for y = 1:size(shum_inv,2)
    if proba_inv_ni(x,y) < proba_rem_sup(x,y) && proba_inv_ni(x,y) > proba_rem_inf(x,y)
    proba_inv_ni(x,y) = NaN;
    else
    end
    if proba_ver_ni(x,y) < proba_rem_sup_ver(x,y) && proba_ver_ni(x,y) > proba_rem_inf_ver(x,y)
    proba_ver_ni(x,y) = NaN;
    else
    end
end
end
subplot(211)
m_proj('mercator','lon',[min(lon1(:)) max(lon1(:))],'lat',[min(lat1(:)) ...
    max(lat1(:))]);
m_pcolor(lon1,lat1,proba_inv_ni')
shading Interp
colormap(jet(10))
colorbar
caxis([0 1]);
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Significancia campo compuesto invierno para la niña')

subplot(212)
m_proj('mercator','lon',[min(lon1(:)) max(lon1(:))],'lat',[min(lat1(:)) ...
    max(lat1(:))]);
m_pcolor(lon1,lat1,proba_ver_ni')
shading Interp
colormap(jet(10))
colorbar
caxis([0 1]);
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Significancia campo compuesto verano para la niña')



 
%NEUTRO
E12_inv_ne=find(n34_inv>-0.5 & n34_inv<0.5);
E12_ver_ne=find(n34_ver>-0.5 & n34_ver<0.5);
for x=1:X
    for y=1:Y
    %invierno
    E12E2_inv_ne= find(shum_inv(x,y,E12_inv_ne)<prctile(squeeze(shum_inv(x,y,:)),75) & shum_inv(x,y,E12_inv_ne)>prctile(squeeze(shum_inv(x,y,:)),25));
    proba_inv_ne(x,y)=length(E12E2_inv_ne)/length(E12_inv_ne);
    %verano
    E12E2_ver_ne= find(shum_ver(x,y,E12_ver_ne)<prctile(squeeze(shum_ver(x,y,:)),75) & shum_ver(x,y,E12_ver_ne)>prctile(squeeze(shum_ver(x,y,:)),25));
    proba_ver_ne(x,y)=length(E12E2_ver_ne)/length(E12_ver_ne);
    end
end

subplot(211)
m_proj('mercator','lon',[min(lon1(:)) max(lon1(:))],'lat',[min(lat1(:)) ...
    max(lat1(:))]);
m_pcolor(lon1,lat1,proba_inv_ne')
shading Interp
colorbar
colormap(jet(10))
caxis([0 1]);
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Probabilidad de que la humedad específica esté sobre el percentil 25 y bajo el percentil 75 debido debido a que ocurra un evento neutro en invierno')

subplot(212)
m_proj('mercator','lon',[min(lon1(:)) max(lon1(:))],'lat',[min(lat1(:)) ...
    max(lat1(:))]);
m_pcolor(lon1,lat1,proba_ver_ne')
shading Interp
colorbar
colormap(jet(10))
caxis([0 1]);
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Probabilidad de que la humedad específica esté sobre el percentil 25 y bajo el percentil 75 debido debido a que ocurra un evento neutro en verano')

% Carlos monte
clear aux1 aux2 probaa_rem probaa_rem_ver ee1 ee2
EN_inv_ne=find(n34_inv>-0.5 & n34_inv<0.5);
EN_ver_ne=find(n34_ver>-0.5 & n34_ver<0.5);
for x=1:X
for y=1:Y
for i=1:100
aux1=remuestreo(shum_ver(x,y,:));
aux2=remuestreo(shum_inv(x,y,:));    
ee1=find(aux2(EN_inv_ne)<prctile(squeeze(shum_inv(x,y,:))',75) & aux2(EN_inv_ne)>prctile(squeeze(shum_inv(x,y,:))',25));
probaa_rem(x,y,i)=length(ee1)/length(EN_inv_ne);
%verano
ee2=find(aux1(EN_ver_ne)<prctile(squeeze(shum_ver(x,y,:))',75) & aux1(EN_ver_ne)>prctile(squeeze(shum_ver(x,y,:))',25));
probaa_rem_ver(x,y,i)=length(ee2)/length(EN_ver_ne);
end
end
end
clear proba_rem_sup proba_rem_inf proba_rem_sup_ver proba_rem_inf_ver
% intervalo de confianza
for x = 1:size(shum_inv,1)
for y = 1:size(shum_inv,2)
proba_rem_sup(x,y) = prctile(squeeze(probaa_rem(x,y,:)),97.5);
proba_rem_inf(x,y) = prctile(squeeze(probaa_rem(x,y,:)),2.5);
%verano
proba_rem_sup_ver(x,y) = prctile(squeeze(probaa_rem_ver(x,y,:)),97.5);
proba_rem_inf_ver(x,y) = prctile(squeeze(probaa_rem_ver(x,y,:)),2.5);
end
end
% Aplicando la signficancia:
clear proba_sig proba_sig_ver
for x = 1:size(shum_inv,1)
    for y = 1:size(shum_inv,2)
    if proba_inv_ne(x,y) < proba_rem_sup(x,y) && proba_inv_ne(x,y) > proba_rem_inf(x,y)
    proba_inv_ne(x,y) = NaN;
    else
    end
    if proba_ver_ne(x,y) < proba_rem_sup_ver(x,y) && proba_ver_ne(x,y) > proba_rem_inf_ver(x,y)
    proba_ver_ne(x,y) = NaN;
    else
    end
end
end
subplot(211)
m_proj('mercator','lon',[min(lon1(:)) max(lon1(:))],'lat',[min(lat1(:)) ...
    max(lat1(:))]);
m_pcolor(lon1,lat1,proba_inv_ne')
shading Interp
colorbar
colormap(jet(10))
caxis([0 1]);
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Significancia campo compuesto invierno para evento neutro')

subplot(212)
m_proj('mercator','lon',[min(lon1(:)) max(lon1(:))],'lat',[min(lat1(:)) ...
    max(lat1(:))]);
m_pcolor(lon1,lat1,proba_ver_ne')
shading Interp
colorbar
colormap(jet(10))
caxis([0 1]);
hold on
m_grid('Box','Fancy','LineStyle','none','FontSize',14);
m_gshhs_i('color','k','linewidth',2);
hold off
title('Significancia campo compuesto verano para evento neutro')
%fin
