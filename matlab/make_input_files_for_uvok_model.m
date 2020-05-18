% Set toplevel path to GCMs configuration
base_path='/sfs/fs1/work-geomar6/smomw258/UVOK_fct_ext';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_2.8deg';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO';
% base_path='/data2/spk/TransportMatrixConfigs/MITgcm_ECCO_v4';

periodicForcing=1
periodicMatrix=1

dt=28800; % time step to use

rearrangeProfiles=1
bigMat=0
writeFiles=1
writeTMs=1
useCoarseGrainedMatrix=0
writePCFiles=0
writeCalibrInput=0

O_npzd_iron=1
O_npzd_fe_limitation=0
O_kk_si=1
READ_SWRAD=1
O_sbg=1

% Set path names, etc.
load(fullfile(base_path,'config_data'))

matrixPath=fullfile(base_path,matrixPath);

explicitMatrixFileBase=fullfile(base_path,explicitMatrixFileBase);
implicitMatrixFileBase=fullfile(base_path,implicitMatrixFileBase);

explicitAnnualMeanMatrixFile=fullfile(base_path,explicitAnnualMeanMatrixFile);
implicitAnnualMeanMatrixFile=fullfile(base_path,implicitAnnualMeanMatrixFile);

preconditionerMatrixFile=fullfile(base_path,preconditionerMatrixFile);

gcmDataPath=fullfile(base_path,'GCM');
bgcDataPath=fullfile(base_path,'BiogeochemData');
% changed path by KK to reflect addition of river discharge
freshWaterForcingFile=fullfile(gcmDataPath,'FreshWaterForcing_gcm');
empFixFile=fullfile(gcmDataPath,empFixFile);

% model specific data
uvokInputDataFile=fullfile(bgcDataPath,'UVOK_input_data');

%
gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','dz','dznom','da','x','y','z','deltaT','gridType')

dtMultiple=dt/deltaT;
if rem(dt,deltaT)
  error('ERROR: Incorrect time step specified! dt must be divisible by deltaT.')
end
disp(['dtMultiple is set to ' num2str(dtMultiple)])

load(boxFile,'Xboxnom','Yboxnom','Zboxnom','izBox','nb','volb')

Ib=find(izBox==1);
nbb=length(Ib);

if rearrangeProfiles || bigMat
  load(profilesFile,'Ip_pre','Ir_pre','Ip_post','Ir_post','Irr')
  Ip=Ip_pre;
  Ir=Ir_pre;
end

if useCoarseGrainedMatrix
  error('NOT FULLY IMPLEMENTED YET!')
end

if periodicForcing
  nm=12;
else
  nm=1;
end

% Use steady state T/S from GCM. Note we always load seasonal data here.
load(fullfile(gcmDataPath,'Theta_gcm'),'Tgcm')
load(fullfile(gcmDataPath,'Salt_gcm'),'Sgcm')
Theta=gridToMatrix(Tgcm,[],boxFile,gridFile);
Salt=gridToMatrix(Sgcm,[],boxFile,gridFile);

% surface area
dab_surf=gridToMatrix(da,Ib,boxFile,gridFile);

% calibration things
if writeCalibrInput
   volFrac=zeros(nb,nbb);
   volFrac=volb/sum(volb); 
   weights=ones(size(volFrac));
% Changed to surface only !!
%   weights=zeros(size(volFrac));
%   weights(Ib)=1;
%   weights(Ibp1)=1;
%   weights(Ibp2)=1;
end

% Surface forcing
useEmP=1; % we always use this
if useEmP
% Compute total E-P for virtual flux:
% Fv = TRg*(E-P), E-P in m/s
% d[TR]/dt = ... + Fv/dz
  load(freshWaterForcingFile,'EmPgcm','Srelaxgcm','saltRelaxTimegcm')
  zeroNetEmP=1;
  EmP=get_surface_emp_for_virtual_flux(gridFile,boxFile,EmPgcm,Srelaxgcm,saltRelaxTimegcm,Salt,zeroNetEmP);
%   gV=EmP./dzb;
  volFracSurf=zeros(nb,1);
  volFracSurf(Ib)=volb(Ib)/sum(volb(Ib));  
%KK change for testing top 3 layers
%volFracNearSurf=zeros(nb,1);
%volFracNearSurf(Ib)=volb(Ib)/(sum(volb(Ib))+sum(volb(Ibp1))+sum(volb(Ibp2)));
%volFracNearSurf(Ibp1)=volb(Ibp1)/(sum(volb(Ib))+sum(volb(Ibp1))+sum(volb(Ibp2)));
%volFracNearSurf(Ibp2)=volb(Ibp2)/(sum(volb(Ib))+sum(volb(Ibp1))+sum(volb(Ibp2)));

  if fixEmP
	load(empFixFile,'empFixX','empFixY')
	nEmPFix=length(empFixX);
%   Points to fix: each SET is individually fixed.
	for k=1:nEmPFix
	  if length(empFixX{k})>1 % polygon  
		Ifix{k}=find(inpolygon(Xboxnom(Ib),Yboxnom(Ib),empFixX{k},empFixY{k})); % indexed to Ib
	  else % single point
		Ifix{k}=find(Xboxnom(Ib)==empFixX{k} & Yboxnom(Ib)==empFixY{k}); % referenced to Ib
	  end
	end
%   Rest of ocean (not fixed)
    Infix=find(~ismember([1:nbb]',cat(1,Ifix{:}))); % Everything else indexed to Ib	
%   Now fix E-P so that annual and spatial integral is 0.
	for k=1:length(Ifix) % loop over each SET of problematic points
	  if useAreaWeighting
		areabfix=dab_surf(Ifix{k});
		EmP(Ifix{k},:) = EmP(Ifix{k},:) - mean(areabfix'*EmP(Ifix{k},:))/sum(areabfix);
	  else
		EmP(Ifix{k},:) = EmP(Ifix{k},:) - mean(mean(EmP(Ifix{k},:)));
	  end
	end
	if useAreaWeighting
	  areabnfix=dab_surf(Infix);
	  EmP(Infix,:) = EmP(Infix,:) - mean(areabnfix'*EmP(Infix,:))/sum(areabnfix);
	else
	  EmP(Infix,:) = EmP(Infix,:) - mean(mean(EmP(Infix,:)));
	end
  end  
end

load(freshWaterForcingFile,'dischgcm')
dischb=gridToMatrix(dischgcm,Ib,boxFile,gridFile,1);

% Grid variables
dzb=gridToMatrix(dz,[],boxFile,gridFile);

if rescaleForcing
  load(fullfile(base_path,rescaleForcingFile))
end

% now take annual mean if necessary
if ~periodicForcing
  Theta=mean(Theta,2);
  Salt=mean(Salt,2);
  if useEmP
    EmP=mean(EmP,2);
  end  
end

if ~periodicMatrix
  if rescaleForcing
	Rfs=mean(Rfs,2);    
  end  
end

% Surface forcing data
load(uvokInputDataFile,'aice')
% acie should be a fraction
aiceb=gridToMatrix(aice,Ib,boxFile,gridFile,1);
if ~periodicForcing
  aiceb=mean(aiceb,2);
end

load(uvokInputDataFile,'hice')
hiceb=gridToMatrix(hice,Ib,boxFile,gridFile,1);
if ~periodicForcing
  hiceb=mean(hiceb,2);
end

load(uvokInputDataFile,'hsno')
hsnob=gridToMatrix(hsno,Ib,boxFile,gridFile,1);
if ~periodicForcing
  hsnob=mean(hsnob,2);
end

load(uvokInputDataFile,'wind')
windb=gridToMatrix(wind,Ib,boxFile,gridFile,1);
if ~periodicForcing
  windb=mean(windb,2);
end

atmospb=load_ocmip_variable([],'P',Xboxnom(Ib),Yboxnom(Ib));  
if ~periodicForcing
  atmospb=mean(atmospb,2);
end

if O_npzd_fe_limitation==1
  load(uvokInputDataFile,'Fe')
  Feb=gridToMatrix(Fe,[],boxFile,gridFile);
  if ~periodicForcing
	Feb=mean(Feb,2);
  end
end
if O_npzd_iron==1
  load(uvokInputDataFile,'Fe_adep')
  Fe_adepb=gridToMatrix(Fe_adep,Ib,boxFile,gridFile,1);
  load(uvokInputDataFile,'Fe_detr_flux')
  Fe_detr_fluxb=gridToMatrix(Fe_detr_flux,Ib,boxFile,gridFile,1);
  load(uvokInputDataFile,'Fe_hydr')
  Fe_hydrb=gridToMatrix(Fe_hydr,[],boxFile,gridFile);
  if ~periodicForcing
	Fe_adepb=mean(Fe_adepb,2);
	Fe_detr_fluxb= mean(Fe_detr_fluxb,2);
  end
end
if O_kk_si==1
  load(uvokInputDataFile,'Si_adep')
  Si_adepb=gridToMatrix(Si_adep,Ib,boxFile,gridFile,1);
  load(uvokInputDataFile,'Si_hydr')
  Si_hydrb=gridToMatrix(Si_hydr,[],boxFile,gridFile);
  if ~periodicForcing
	Si_adepb=mean(Si_adepb,2);
  end
end
if O_sbg==1
  load(uvokInputDataFile,'sgbathy')
  sgbathyb=gridToMatrix(sgbathy,[],boxFile,gridFile);
  Ibot=cellfun(@max,Ip); % indices of bottom cell
  sgbathyb(Ibot)=1; % required for PO4 conservation
end
if READ_SWRAD
% Read SW radiation from file
  load(uvokInputDataFile,'swrad')
  swradb=gridToMatrix(swrad,Ib,boxFile,gridFile,1);
  if ~periodicForcing
	swradb=mean(swradb,2);
  end  
end

latb=Yboxnom(Ib);

% Initial condition
trNames={'idealage','dic','c14','alk','o2','po4','phyt','zoop','detr','no3',...
'diaz','dfe','detrfe','cocc','caco3','detr_B','diat','sil'};  

numTracers=length(trNames);
TRini=zeros(nb,numTracers);
for itr=1:numTracers
  fn=fullfile('InitialConditionProfiles',[trNames{itr} '.dat']);
  [hdr,dat]=hdrload(fn);
  TRini(:,itr)=interp1(dat(:,1),dat(:,2),Zboxnom);
  kk=find(Zboxnom<dat(1,1));
  TRini(kk,itr)=dat(1,2);
  kk=find(Zboxnom>dat(end,1));
  TRini(kk,itr)=dat(end,2);  
%%%%%%%%%%%
  if writeCalibrInput
   if itr==1
    fn=fullfile('obs_for_calibration','O_totcarb.nc');
    TRini(:,1) = gridToMatrix(double(ncread(fn,'O_totcarb')),[],boxFile,gridFile); 
    TRDIC = TRini(:,1);
   elseif itr==2
    fn=fullfile('obs_for_calibration','O_dc14.nc');
    TRini(:,2) = gridToMatrix(double(ncread(fn,'O_dc14')),[],boxFile,gridFile);
% convert units
    TRini(:,2) = (TRini(:,2)*1e-3 + 1)*1.176e-12.*TRDIC;
    elseif itr==17
    fn=fullfile('obs_for_calibration','O_silica.nc');
    TRini(:,17) = gridToMatrix(double(ncread(fn,'O_silica')),[],boxFile,gridFile); 
   else
    fn=fullfile('obs_for_calibration',['O_' trNames{itr} '.nc']);
    if exist(fn,'file') 
     TRini(:,itr) = gridToMatrix(double(ncread(fn,['O_' trNames{itr}])),[],boxFile,gridFile);
    end
    fn = 0.;
   end
  end
%%%%%%%%%%
  if any(isnan(TRini(:,itr)))
    error('ERROR interpolating initial conditions!')
  end
end

if rearrangeProfiles
  Xboxnom=Xboxnom(Ir);
  Yboxnom=Yboxnom(Ir);
  Zboxnom=Zboxnom(Ir);
  izBox=izBox(Ir);
  Theta=Theta(Ir,:);
  Salt=Salt(Ir,:);
  if O_npzd_fe_limitation==1
	Feb=Feb(Ir,:);
  end
  if O_npzd_iron==1
    Fe_hydrb=Fe_hydrb(Ir,:);
  end  	
  if O_kk_si==1
    Si_hydrb=Si_hydrb(Ir,:);
  end
  if O_sbg==1  	
   sgbathyb=sgbathyb(Ir);
  end
  TRini=TRini(Ir,:);
  if useEmP
    volFracSurf=volFracSurf(Ir);
%    volFracNearSurf=volFracNearSurf(Ir);
  end  
  if writeCalibrInput
   volFrac=volFrac(Ir);
   weights=weights(Ir);
  end
  dzb=dzb(Ir);
  if rescaleForcing
    Rfs=Rfs(Ir,:);
  end    
  Ib=find(izBox==1);
%
  Ip=Ip_post;
  Ir=Ir_post;
end  

if useCoarseGrainedMatrix
% Coarse grain initial conditions and forcing data
end

if writeFiles
  calc_periodic_times_for_tmm('monthly-365-day year','periodic_times_365d.bin');
  calc_periodic_times_for_tmm('monthly-360-day year','periodic_times_360d.bin');  
% Transport matrices
  if writeTMs
%   Explicit transport matrix
	I=speye(nb,nb);
	if ~periodicMatrix
      disp('loading annual mean explicit TM')	
      load(explicitAnnualMeanMatrixFile,'Aexpms')	
	  if rearrangeProfiles
		Aexpms=Aexpms(Ir_pre,Ir_pre); % rearrange
	  end      
	  % make discrete
	  Aexpms=dt*Aexpms;
	  Aexpms=I+Aexpms;
	  if useCoarseGrainedMatrix
		Aexpms=Beta*Aexpms*M; % coarse-grained explicit transport matrix
	  end        
	  writePetscBin('Ae.petsc',Aexpms,[],1)
	else
      % load each month from separate file
      disp('loading monthly mean explicit TMs')	      
	  for im=1:12 
		fn=[explicitMatrixFileBase '_' sprintf('%02d',im)];
		load(fn,'Aexp')
		if rearrangeProfiles
		  Aexp=Aexp(Ir_pre,Ir_pre); % rearrange
		end
		% make discrete
		Aexp=dt*Aexp;
		Aexp=I+Aexp;
		if useCoarseGrainedMatrix
%         Not sure if this is really kosher!		  
		  Aexp=Beta*Aexp*M; % coarse-grained explicit transport matrix
		end
		writePetscBin(['Ae_' sprintf('%02d',im-1)],Aexp,[],1)
		clear Aexp
	  end
	end
%   Implicit transport matrix
	if ~periodicMatrix
      disp('loading annual mean implicit TM')		
      load(implicitAnnualMeanMatrixFile,'Aimpms')
      if dtMultiple~=1
		if bigMat % big matrix. do it a block at a time.
		  for is=1:nbb % change time step multiple
			Aimpms(Ip_pre{is},Ip_pre{is})=Aimpms(Ip_pre{is},Ip_pre{is})^dtMultiple;
		  end
		else
		  Aimpms=Aimpms^dtMultiple;
		end  
	  end	
	  if rearrangeProfiles
		Aimpms=Aimpms(Ir_pre,Ir_pre); % rearrange
	  end
	  if useCoarseGrainedMatrix
		Aimpms=Beta*Aimpms*M; % coarse-grained implicit transport matrix
	  end
	  writePetscBin('Ai.petsc',Aimpms,[],1)
	else
	  % load each month from separate file
      disp('loading monthly mean implicit TMs')	      	  
	  for im=1:12
		fn=[implicitMatrixFileBase '_' sprintf('%02d',im)];		
		load(fn,'Aimp')
		if dtMultiple~=1
		  if bigMat % big matrix. do it a block at a time.		
			for is=1:nbb % change time step multiple
			  Aimp(Ip_pre{is},Ip_pre{is})=Aimp(Ip_pre{is},Ip_pre{is})^dtMultiple;
			end
		  else
			Aimp=Aimp^dtMultiple;		
		  end
		end  
		if rearrangeProfiles
		  Aimp=Aimp(Ir_pre,Ir_pre); % rearrange
		end
		if useCoarseGrainedMatrix
		  Aimp=Beta*Aimp*M; % coarse-grained implicit transport matrix		
		end
		writePetscBin(['Ai_' sprintf('%02d',im-1)],Aimp,[],1)
		clear Aimp
	  end
	end
  end	  	  
% Initial conditions  
  for itr=1:numTracers
    fn=[trNames{itr} 'ini.petsc'];
	writePetscBin(fn,TRini(:,itr))  
  end

  dischb=dischb/10; % convert to UVOK units kg/m^2/s -> g/cm^2/s  

  if ~periodicForcing
	write_binary('aice.bin',aiceb,'real*8')
	write_binary('hice.bin',hiceb,'real*8')
	write_binary('hsno.bin',hsnob,'real*8')	
	write_binary('wind.bin',windb,'real*8')	
	write_binary('atmosp.bin',atmospb,'real*8')
	write_binary('disch.bin',dischb,'real*8')	
  else
    for im=1:nm
	  write_binary(['aice_' sprintf('%02d',im-1)],aiceb(:,im),'real*8')
	  write_binary(['hice_' sprintf('%02d',im-1)],hiceb(:,im),'real*8')
	  write_binary(['hsno_' sprintf('%02d',im-1)],hsnob(:,im),'real*8')
	  write_binary(['wind_' sprintf('%02d',im-1)],windb(:,im),'real*8')	  
	  write_binary(['atmosp_' sprintf('%02d',im-1)],atmospb(:,im),'real*8')	  
	  write_binary(['disch_' sprintf('%02d',im-1)],dischb(:,im),'real*8')
	end
  end
  if READ_SWRAD
	if ~periodicForcing  
	  write_binary('swrad.bin',swradb,'real*8')
	else
	  for im=1:nm	
		write_binary(['swrad_' sprintf('%02d',im-1)],swradb(:,im),'real*8')
	  end  
    end	
  end
  Salt=(Salt-35)/1000; % change to UVOK units
  if ~periodicForcing
	writePetscBin('Ts.petsc',Theta)
	writePetscBin('Ss.petsc',Salt)
  else
    for im=1:nm
	  writePetscBin(['Ts_' sprintf('%02d',im-1)],Theta(:,im))
	  writePetscBin(['Ss_' sprintf('%02d',im-1)],Salt(:,im))
    end    
  end    
  if O_npzd_fe_limitation==1
	if ~periodicForcing
	  writePetscBin('Fe_dissolved.petsc',Feb)
	else
	  for im=1:nm
		writePetscBin(['Fe_dissolved_' sprintf('%02d',im-1)],Feb(:,im))
	  end    
	end      
  end  
  if O_npzd_iron==1
	writePetscBin('Fe_hydr.petsc',Fe_hydrb)
	if ~periodicForcing
	   write_binary('Fe_adep.bin',Fe_adepb,'real*8')
	   write_binary('Fe_detr_flux.bin',Fe_detr_fluxb,'real*8')
	else
	  for im=1:nm
	    write_binary(['Fe_adep_' sprintf('%02d',im-1)],Fe_adepb(:,im),'real*8')
	    write_binary(['Fe_detr_flux_' sprintf('%02d',im-1)],Fe_detr_fluxb(:,im),'real*8')
	  end    
	end      
  end  
  if O_kk_si==1
	writePetscBin('Si_hydr.petsc',Si_hydrb)
	if ~periodicForcing
	   write_binary('Si_adep.bin',Si_adepb,'real*8')
	else
	  for im=1:nm
	    write_binary(['Si_adep_' sprintf('%02d',im-1)],Si_adepb(:,im),'real*8')
	  end    
	end      
  end  
  if O_sbg==1
    writePetscBin('sgbathy.petsc',sgbathyb)  
  end
  if useEmP
    EmP=EmP*100; % change to UVOK units (cm/s)
    writePetscBin('surface_volume_fraction.petsc',volFracSurf)
    if ~periodicForcing
  	  write_binary('EmP.bin',EmP,'real*8')
    else
	  for im=1:nm
  	    write_binary(['EmP_' sprintf('%02d',im-1)],EmP(:,im),'real*8')
	  end
	end    
  end
    if writeCalibrInput
   writePetscBin('volume_fraction.petsc',volFrac)
   writePetscBin('weights.petsc',weights)
  end
  if rescaleForcing
	if ~periodicMatrix
	  writePetscBin('Rfs.petsc',Rfs)
	else
	  for im=1:nm
		writePetscBin(['Rfs_' sprintf('%02d',im-1)],Rfs(:,im))
	  end    
	end    
  end  
  
% Grid data
  z=z*100; % convert to cm
  dznom=dznom*100; % convert to cm
  dzb=dzb*100; % convert to cm  
  dab_surf=dab_surf*1e4; % convert to cm^2
  write_binary('zt.bin',nz,'int')  
  write_binary('zt.bin',z,'real*8',1)
  write_binary('drF.bin',nz,'int')  
  write_binary('drF.bin',dznom,'real*8',1)
  writePetscBin('dz.petsc',dzb)
  write_binary('latitude.bin',latb,'real*8')
  write_binary('dA.bin',dab_surf,'real*8')    
% Profile data  
  if rearrangeProfiles
    if ~useCoarseGrainedMatrix
	  gStartIndices=cellfun(@(x)x(1),Ip);
	  gEndIndices=cellfun(@(x)x(end),Ip);
    else % useCoarseGrainedMatrix
	  gStartIndices=cellfun(@(x)x(1),Ipcg);
	  gEndIndices=cellfun(@(x)x(end),Ipcg);
    end  
    write_binary('gStartIndices.bin',[length(gStartIndices);gStartIndices],'int')
    write_binary('gEndIndices.bin',[length(gEndIndices);gEndIndices],'int')
  end
end

if useCoarseGrainedMatrix
  numProfiles=nbbcg;
else  
  numProfiles=nbb;
end
disp(['Number of Profiles in this Configuration: ' int2str(numProfiles)])

if writePCFiles
  pc=load(preconditionerMatrixFile,'Aexpms');
  if rearrangeProfiles
    A=pc.Aexpms(Ir_pre,Ir_pre);
  else
    A=pc.Aexpms;
  end
  clear pc  
  if useCoarseGrainedMatrix
    A=Beta*A*M;
    save pc_cg_data A nbbcg CGgrid CG Ipcg Ibcg dt
  else
    save pc_data A nbb nz nb Ip Ib dt
  end
end
