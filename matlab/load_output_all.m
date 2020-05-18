base_path='/sfs/fs1/work-geomar6/smomw258/UVic_matrix_iron_test';

load(fullfile(base_path,'config_data'))

matrixPath=fullfile(base_path,matrixPath);

gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','x','y','z','gridType')

load(profilesFile,'Irr')

timeFile='output_time.txt';

trNames={'dic','c14','alk','o2','po4','phyt','zoop','detr','no3',...
'diaz','dfe','detrfe','cocc','caco3','detr_B','diat','sil'};   

numTr=length(trNames);

if strcmp(gridType,'llc_v4')
  load(fullfile(base_path,'llc_v4_grid'))
  gcmfaces_global
end

[hdr,tdat]=hdrload(timeFile);
T=tdat(:,2);

for itr=1:numTr
  varName=upper(trNames{itr})
  fn=[trNames{itr} '.petsc'];
  tr=readPetscBinVec(fn,-1);
  TR=matrixToGrid(tr(Irr,:),[],boxFile,gridFile);
  nt=size(TR,4); 
  if strcmp(gridType,'llc_v4')
    varName=[varName '_plot'];
	tmp=gcmfaces(TR);
	[x,y,TRplot]=convert2pcol(mygrid.XC,mygrid.YC,tmp);
	[n1,n2]=size(TRplot);
	eval([trNames{itr} '=zeros([n1 n2 nz nt]);']);
	for it=1:nt
	  for iz=1:nz
		eval(['[x,y,' varName '(:,:,iz,it)]=convert2pcol(mygrid.XC,mygrid.YC,tmp(:,:,iz,it));']);
	  end
	end  
  else
	eval([varName '=TR;']);
  end
  save(varName,varName,'x','y','z','T')
  write2netcdf([varName '.nc'],TR,x,y,z,T,upper(varName))
end
