C$Header: /Users/ikriest/CVS/mops/mops_biogeochem_misfit.F,v 1.1 2015/11/17 14:18:51 ikriest Exp $
C$Name: mops-2_0 $

C !INTERFACE: ==========================================================
      SUBROUTINE UVOK_BIOGEOCHEM_MISFIT(kmt_loc,
     &           mbgc1loc,mbgc2loc,mbgc3loc)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE
#include "size.h"
#include "mw.h"
!#include "BGC_MISFIT.h"

C !INPUT PARAMETERS: ===================================================
      INTEGER kmt_loc
C     Output variables
      real*8 mbgc1loc(kmt_loc),mbgc2loc(kmt_loc),mbgc3loc(kmt_loc)

C !LOCAL VARIABLES: ====================================================
!      INTEGER k
!      print*,'KK kmt_loc',kmt_loc
C     Copy diagnostic variables
!      DO k=1,kmt_loc
        mbgc1loc(1:kmt_loc)=t(1,1:kmt_loc,1,ipo4,taum1)
!        print*,'KK made it in',mbgc1loc(1:kmt_loc)
        mbgc2loc(1:kmt_loc)=t(1,1:kmt_loc,1,io2,taum1)
        mbgc3loc(1:kmt_loc)=t(1,1:kmt_loc,1,ino3,taum1)
!      ENDDO

      RETURN
      END
