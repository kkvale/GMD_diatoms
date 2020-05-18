base_path='/sfs/fs1/work-geomar6/smomw258/UVic_matrix_iron_test';

load(fullfile(base_path,'config_data'))

matrixPath=fullfile(base_path,matrixPath);

gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','x','y','z','gridType')

load(profilesFile,'Irr')

trNames={'dic','c14','alk','o2','po4','phyt','zoop','detr','no3',...
	 'diaz','dfe','detrfe','cocc','caco3','detr_B','diat','sil'};  
dischfile={'disch_00','disch_01','disch_02','disch_03','disch_04',...
	   'disch_05','disch_06','disch_07','disch_08','disch_09',...
	   'disch_10','disch_11'};
numTr=length(trNames);

if strcmp(gridType,'llc_v4')
  load(fullfile(base_path,'llc_v4_grid'))
  gcmfaces_global
end

for itr=1:numTr
  varName=upper(trNames{itr})
  fn=[trNames{itr} 'ini.petsc'];
  tr=readPetscBinVec(fn,1,-1);
  TR=matrixToGrid(tr(Irr,end),[],boxFile,gridFile);

  if strcmp(gridType,'llc_v4')
    varName=[varName '_plot'];
	tmp=gcmfaces(TR);
	[x,y,TRplot]=convert2pcol(mygrid.XC,mygrid.YC,tmp);
	[n1,n2]=size(TRplot);
	eval([trNames{itr} '=zeros([n1 n2 nz]);']);
	for iz=1:nz
	  eval(['[x,y,' varName '(:,:,iz)]=convert2pcol(mygrid.XC,mygrid.YC,tmp(:,:,iz));']);
	end
  else
	eval([varName '=TR;']);
  end
  save(varName,varName,'x','y','z')
  write2netcdf([varName 'ini.nc'],TR,x,y,z,[],upper(varName))
end
    
    load(boxFile,'izBox','nb','volb')
    Ib=find(izBox==1);
    nbb=length(Ib);
for itr=1:12
  varName='dischb';
  fn=[dischfile{itr}];
  tr=read_binary(fn,[],'float64');
  nt=length(tr)/nbb;
  tr=reshape(tr,[nbb,nt]);
  TR=matrixToGrid(tr,Ib,boxFile,gridFile);
	eval([varName '=TR;']);
  save(varName,varName,'x','y')
  write2netcdf([fn 'ini.nc'],TR,x,y,[],[],upper(varName),[],[])
end
