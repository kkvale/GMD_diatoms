base_path=fileparts(pwd);

gridFile=fullfile(base_path,'grid');
load(gridFile,'bathy','nx','ny','nz')

bathys=bathy(:,:,1);
bathys=repmat(bathys,[1 1 12]);

aice=double(ncread('tavg.nc','O_icefra')); % fraction
aice=bathys.*aice;

hice=double(ncread('tavg.nc','O_icethk')); % m
hice=hice*100; % convert m to cm
hice=bathys.*hice;

hsno=double(ncread('tavg.nc','O_snothk')); % m
hsno=hsno*100; % convert m to cm
hsno=bathys.*hsno;

swrad=double(ncread('tavg.nc','F_dnswr')); % W/m^2
swrad=swrad*1000; % convert W/m^2 (kg/s^3) to g/s^3
swrad=bathys.*swrad;

wind=double(ncread('tavg.nc','A_windspd')); % m/s
wind=wind*100; % convert m/s to cm/s
wind=bathys.*wind;

% dissolved Fe concentration when not using prognostic iron model
tmp=double(ncread(fullfile('UVOK_in','O_fe_dissolved.nc'),'O_dissolved_fe'));
Fe=zeros([nx ny nz 12]);
n3=size(tmp,3);
for it=1:12
  Fe(:,:,1:n3,it)=tmp(:,:,:,it);
  Fe(:,:,:,it)=Fe(:,:,:,it).*bathy;
end

% data for prognostic iron cycle
% Atmospheric Fe deposition
feadepFile=fullfile('UVOK_in','O_feflux.nc');
Fe_adep=squeeze(double(ncread(feadepFile,'O_feflux')));

% data for prognostic Silica cycle
% Atmospheric Si deposition
siadepFile=fullfile('UVOK_in','O_sil_dep.nc');
Si_adep=squeeze(double(ncread(siadepFile,'O_SIL_DEP')));

% Detrital Fe flux
Fe_detr_flux=zeros([nx ny 12]);

% Hydrothermal Fe input
fehydrFile=fullfile('UVOK_in','O_fe_hydr.nc');
Fe_hydr=squeeze(double(ncread(fehydrFile,'O_fe_hydr')));

% Hydrothermal Si input
sihydrFile=fullfile('UVOK_in','O_si_hydr.nc');
Si_hydr=squeeze(double(ncread(sihydrFile,'O_SIL_HYDR')));

sgbathy=double(ncread(fullfile('UVOK_in','G_subgrid_bathy.nc'),'SUBGRID_BATHY'));

save UVOK_input_data aice hice hsno wind Si_adep Si_hydr Fe Fe_adep Fe_hydr Fe_detr_flux sgbathy swrad
%save UVOK_input_data aice hice hsno wind Fe swrad
