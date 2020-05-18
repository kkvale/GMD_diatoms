#ifdef O_TMM
      SUBROUTINE UVOK_SED_COPY_DATA(sedmixed, sedburied, sedstore, dir)
      
      implicit none

#include "size.h"
#include "sed.h"

      real sedmixed(nzmax,20,ipmax)
      real sedburied(ibmax,2,ipmax)
      real sedstore(ipmax,100)
      integer dir

#include "pconst.h"
#include "mw.h"

      if (dir.eq.1) then
!      write(*,*)'SPK1: copying from TMM to SED'
C       copy TMM tracer field to SED tracer array        
        carb(:,1,1:ipmax) = sedmixed(:,1,1:ipmax)
        carb(:,2,1:ipmax) = sedmixed(:,2,1:ipmax)
        carb(:,3,1:ipmax) = sedmixed(:,3,1:ipmax)
        dcpls(:,1,1:ipmax) = sedmixed(:,4,1:ipmax)
        dcpls(:,2,1:ipmax) = sedmixed(:,5,1:ipmax)
        dcpls(:,3,1:ipmax) = sedmixed(:,6,1:ipmax)
        dcmin(:,1,1:ipmax) = sedmixed(:,7,1:ipmax)
        dcmin(:,2,1:ipmax) = sedmixed(:,8,1:ipmax)
        dcmin(:,3,1:ipmax) = sedmixed(:,9,1:ipmax)
C
        pore(:,1:ipmax) = sedmixed(:,10,1:ipmax)
        form(:,1:ipmax) = sedmixed(:,11,1:ipmax)
        o2(:,1:ipmax) = sedmixed(:,12,1:ipmax)
        orggg(:,1:ipmax) = sedmixed(:,13,1:ipmax)
        orgml(:,1:ipmax) = sedmixed(:,14,1:ipmax)
        calgg(:,1:ipmax) = sedmixed(:,15,1:ipmax)
        calml(:,1:ipmax) = sedmixed(:,16,1:ipmax)
        dopls(:,1:ipmax) = sedmixed(:,17,1:ipmax)
        domin(:,1:ipmax) = sedmixed(:,18,1:ipmax)
        dbpls(:,1:ipmax) = sedmixed(:,19,1:ipmax)
        dbmin(:,1:ipmax) = sedmixed(:,20,1:ipmax)
C
C       write(*,*)'SPK-copy:',nzmax,ibmax,ipmax
        buried_mass(:,1:ipmax) = sedburied(:,1,1:ipmax)
        buried_calfrac(:,1:ipmax) = sedburied(:,2,1:ipmax)
C
        zrct(1:ipmax)=sedstore(1:ipmax,1)
        k1(1:ipmax)=sedstore(1:ipmax,2)
        k2(1:ipmax)=sedstore(1:ipmax,3)
        k3(1:ipmax)=sedstore(1:ipmax,4)
        csat(1:ipmax)=sedstore(1:ipmax,5)
        rc(1:ipmax)=sedstore(1:ipmax,6)
        ttrorg(1:ipmax)=sedstore(1:ipmax,7)
        ttrcal(1:ipmax)=sedstore(1:ipmax,8)
        sed_ml_mass(1:ipmax)=sedstore(1:ipmax,9)
        c_advect(1:ipmax)=sedstore(1:ipmax,10)
        rain_org_p(1:ipmax)=sedstore(1:ipmax,11)
        rain_cal_p(1:ipmax)=sedstore(1:ipmax,12)
        co3_p(1:ipmax)=sedstore(1:ipmax,13)        
!      write(*,*)'SPK2: finished copying from TMM to SED'		
      elseif (dir.eq.2) then
C       copy SED tracer field to TMM tracer array
!      write(*,*)'SPK3: copying from SED to TMM'
        sedmixed(:,1,1:ipmax) = carb(:,1,1:ipmax)
        sedmixed(:,2,1:ipmax) = carb(:,2,1:ipmax)
        sedmixed(:,3,1:ipmax) = carb(:,3,1:ipmax)
        sedmixed(:,4,1:ipmax) = dcpls(:,1,1:ipmax)
        sedmixed(:,5,1:ipmax) = dcpls(:,2,1:ipmax)
        sedmixed(:,6,1:ipmax) = dcpls(:,3,1:ipmax)
        sedmixed(:,7,1:ipmax) = dcmin(:,1,1:ipmax)
        sedmixed(:,8,1:ipmax) = dcmin(:,2,1:ipmax)
        sedmixed(:,9,1:ipmax) = dcmin(:,3,1:ipmax)
C
        sedmixed(:,10,1:ipmax) = pore(:,1:ipmax)
        sedmixed(:,11,1:ipmax) = form(:,1:ipmax)
        sedmixed(:,12,1:ipmax) = o2(:,1:ipmax)
        sedmixed(:,13,1:ipmax) = orggg(:,1:ipmax)
        sedmixed(:,14,1:ipmax) = orgml(:,1:ipmax) 
        sedmixed(:,15,1:ipmax) = calgg(:,1:ipmax) 
        sedmixed(:,16,1:ipmax) = calml(:,1:ipmax) 
        sedmixed(:,17,1:ipmax) = dopls(:,1:ipmax) 
        sedmixed(:,18,1:ipmax) = domin(:,1:ipmax) 
        sedmixed(:,19,1:ipmax) = dbpls(:,1:ipmax) 
        sedmixed(:,20,1:ipmax) = dbmin(:,1:ipmax) 
C
        sedburied(:,1,1:ipmax) = buried_mass(:,1:ipmax) 
        sedburied(:,2,1:ipmax) = buried_calfrac(:,1:ipmax)
C
        sedstore(1:ipmax,1)=zrct(1:ipmax)
        sedstore(1:ipmax,2)=k1(1:ipmax)
        sedstore(1:ipmax,3)=k2(1:ipmax)
        sedstore(1:ipmax,4)=k3(1:ipmax)
        sedstore(1:ipmax,5)=csat(1:ipmax)
        sedstore(1:ipmax,6)=rc(1:ipmax)
        sedstore(1:ipmax,7)=ttrorg(1:ipmax)
        sedstore(1:ipmax,8)=ttrcal(1:ipmax)
        sedstore(1:ipmax,9)=sed_ml_mass(1:ipmax)
        sedstore(1:ipmax,10)=c_advect(1:ipmax)
        sedstore(1:ipmax,11)=rain_org_p(1:ipmax)
        sedstore(1:ipmax,12)=rain_cal_p(1:ipmax)
        sedstore(1:ipmax,13)=co3_p(1:ipmax)
         
!      write(*,*)'SPK4: finished copying from SED to TMM'
      else
        stop 'ERROR in UVOK_SED_COPY_DATA: unknown copy option!'
      endif

      END
#endif      