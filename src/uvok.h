extern void uvok_copy_data_(PetscInt *nzloc, PetscInt *itr, PetscScalar localTR[], PetscInt *direction);

extern void uvok_sed_copy_data_(PetscScalar localSedMixedTR[], PetscScalar localSedBuriedTR[], PetscScalar localSedStore[], PetscInt *direction);

extern void uvok_ini_(PetscInt *nzmax, PetscScalar zt[], PetscScalar drF[], PetscScalar *DeltaT, 
               PetscScalar *Sglobavg,PetscScalar TRglobavg[], 
# if defined O_TMM_opt
                PetscInt *nIter,  PetscInt *nSession, 
# endif
#ifdef O_sed  
               PetscInt *nzmaxSed, PetscInt *ibmaxSed, PetscInt *numSedMixedTracers, 
               PetscInt *numSedBuriedTracers, PetscScalar *globalweathflx,
#endif               
               PetscInt *debugFlag);

extern void uvok_calc_(PetscInt *nzloc, PetscScalar *locallatitude, PetscScalar *day, PetscScalar *relyr,
               PetscScalar localTs[], PetscScalar localSs[], PetscScalar TRglobavg[], PetscScalar localdz[],
               PetscScalar zt[], PetscScalar *localarea,
# if defined O_carbon
#if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
               PetscScalar *pCO2atm,
#endif               
               PetscScalar *localwind,
#endif
# if defined O_c14ccn_data
               PetscScalar *dc14ccnnatm, PetscScalar *dc14ccnsatm, PetscScalar *dc14ccneatm,
#endif            
#  if defined O_npzd_subgridbathy
               PetscScalar localsgbathy[],
#  endif
#  if defined O_npzd_fe_limitation
               PetscScalar localFe_dissolved[],
#  endif
#ifdef O_npzd_iron
	       PetscScalar *localFe_adep, PetscScalar *localFe_detr_flux, PetscScalar localFe_hydr[],
#endif
#ifdef O_kk_si
	       PetscScalar *localSi_adep, PetscScalar localSi_hydr[],
#endif
#  if defined O_embm
               PetscScalar *localswrad,
#  endif
#  if defined O_ice
#   if !defined O_ice_cpts
               PetscScalar *localaice, PetscScalar *localhice, PetscScalar *localhsno,
#   endif
#  endif
               PetscScalar *localEmP, PetscScalar *empglobavg, 
# if defined O_carbon
               PetscScalar *gasexfluxloc, PetscScalar *totfluxloc, 
# endif           
# if defined O_embm
               PetscScalar *globdisch, PetscScalar *localdisch,
#  if defined O_kk_si
               PetscScalar *globsilwflx, PetscScalar *localsilwflx,
#  endif     
# endif
# if defined O_sed
               PetscScalar *globwflx, PetscScalar *localwflx,
# endif
               PetscInt *debugFlag);

extern void uvok_diags_ini_(PetscInt *lNumProfiles, PetscInt *lTotNumPoints, PetscInt *lNum2dDiags,
               PetscInt *lNum3dDiags, PetscInt *debugFlag);

extern void uvok_diags_start_(PetscInt *debugFlag);

extern void uvok_diags_stop_(PetscInt *debugFlag);

extern void uvok_diags_accumulate_(PetscInt *ipro, PetscInt *kl, PetscInt *nzloc, PetscInt *numAvg, PetscInt *avgFlag, PetscInt *debugFlag);

extern void uvok_diags2d_copy_(PetscInt *id, PetscScalar diagArr[], char *fname, PetscInt *debugFlag);

extern void uvok_diags3d_copy_(PetscInt *id, PetscScalar diagArr[], char *fname, PetscInt *debugFlag);

extern void uvok_diags_finalize_(PetscInt *debugFlag);

#if !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) 
#define uvok_copy_data_ uvok_copy_data
#define uvok_sed_copy_data_ uvok_sed_copy_data
#define uvok_ini_ uvok_ini
#define uvok_start_ uvok_start
#define uvok_stop_ uvok_stop
#define uvok_calc_ uvok_calc
#define uvok_diags_ini_ uvok_diags_ini
#define uvok_diags_accumulate_ uvok_diags_accumulate
#define uvok_diags2d_copy_ uvok_diags2d_copy
#define uvok_diags3d_copy_ uvok_diags3d_copy
#define uvok_diags_finalize_ uvok_diags_finalize
#endif 

