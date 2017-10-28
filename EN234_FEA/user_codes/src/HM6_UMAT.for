
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      INCLUDE 'ABA_PARAM.INC'
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

!     Local variables

      double precision :: edev(6)
      double precision :: matrix_a(6,6)
      double precision :: matrix_b(6,6)
      double precision :: evol
      double precision :: se,ee
      double precision :: g, poi,e0,pt,bulk_modulus,kb
      double precision :: E,xnu
      integer :: i,j
      
       evol = sum(STRAN(1:3)+DSTRAN(1:3))
       edev(1:3) = STRAN(1:3)+DSTRAN(1:3) - evol/3.d0
       edev(4:6) = 0.5d0*(STRAN(4:6)+DSTRAN(4:6))

       g = PROPS(1)
       poi = PROPS(2)
       e0 = PROPS(3)
       pt = PROPS(4)
       bulk_modulus = pt*(1.d0-2.d0*poi)/(1.d0+poi)*(1.d0+e0)
       kb = pt*(1.d0+e0)/3.d0/bulk_modulus*
     1 EXP(-(1.d0+e0)/bulk_modulus*evol)
     2 -2.d0/3.d0*g
       
       matrix_a(1:6,1:6)=0.d0
       matrix_b(1:6,1:6)=0.d0
       
       matrix_a(1,1)=2
       matrix_a(2,2)=2
       matrix_a(3,3)=2
       matrix_a(4,4)=1
       matrix_a(5,5)=1
       matrix_a(6,6)=1
       
       matrix_b(1:3,1:3) = matrix_b(1:3,1:3)+1.d0
       
       DDSDDE(1:6,1:6) = 0.d0
       
       DDSDDE(1:6,1:6) = g*matrix_a(1:6,1:6)+kb*matrix_b(1:6,1:6)

       stress(1:6) = 2.d0*g*edev(1:6)
       stress(1:3) = stress(1:3)+pt/3.d0*
     1 (1-EXP(-(1.d0+e0)/bulk_modulus*evol))
! for test



!
       return

      RETURN
      END SUBROUTINE UMAT