!    EN2340, Siyuan_Song Nov-16-2017
!    God bless me.
!    ABAQUS format user material subroutine for explicit dynamics
!
!

      subroutine vumat(nblock, ndir, nshr, nstatev, nfieldv, nprops,
     1  lanneal, stepTime, totalTime, dt, cmname, coordMp, charLength,
     2  props, density, strainInc, relSpinInc,
     3  tempOld, stretchOld, defgradOld, fieldOld,
     4  stressOld, stateOld, enerInternOld, enerInelasOld,
     5  tempNew, stretchNew, defgradNew, fieldNew,
     6  stressNew, stateNew, enerInternNew, enerInelasNew )
!
!      include 'vaba_param.inc'
!
      implicit double precision (a-h,o-z)

      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     9  defgradNew(nblock,ndir+nshr+nshr),
     1  fieldNew(nblock,nfieldv),
     2  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3  enerInternNew(nblock), enerInelasNew(nblock)

       character*80 cmname
!
!      Local variables
!
       integer k,ntens,nit,maxit

       double precision dedev(ndir+nshr),sdevstar(ndir+nshr)
       double precision sestar,skkstar
       double precision E,xnu,Y,e0,m,edot0,S,H,td,omega,alfa
       double precision eplas,deplas
       double precision deplas_2
       double precision f,dfde
       double precision f_2
       double precision ta, dta, c1, c2, c3
       double precision ta_2, dta_2, c1_2, c2_2, c3_2
!
!      Conventions for storing tensors:
!
!      Deformation gradients are provided as components in the global basis
!      (ABAQUS also allows the user to define a fixed local coord system defined with *ORIENTATION)
!
!      Stresses and stretches are stored as components in a basis that 'rotates with the material'
!      These are defined as follows:
!          Let e_i be the initial basis vectors (the global ijk directions, or user-specified basis vectors)
!          Let R be the rotation tensor, defined through the polar decomposition of deformation gradient F=RU=VR
!          Then define rotated basis vectors m_i = R e_i.
!          The co-rotational components of a tensor are defined as A_ij = m_i . A . m_j
!          The components A_ij can also be interpreted as the global components of the tensor  R^T A R
!
!
!      The components of symmetric tensors (stress,strainInc and stretch) are stored as vectors in the following order
!                      For 3D problems (ndir=3,nshr=3) [11,22,33,12,23,13]
!                      For 2D problems (ndir=3,nshr=1) [11,22,33,12]
!      The components of unsymmetric tensors (defgrad and relSpinInc) are stored as vectors in the following order
!                      For 3D problems (ndir=3,nshr=3) [11,22,33,12,23,31,21,32,13]
!                      For 2D problems (ndir=3,nshr=1) [11,22,33,12,21]
!
!      The stresses are Cauchy (true) stress
!
!      nblock                   No. integration points in data block (data are provided in 'blocks' of integration points)
!                               EN234FEA always has only 1 int pt in each block
!      ndir                     No. direct tensor components
!      nshr                     No. shear tensor components
!      nstatev                  No. user defined state variables (declared in input file)
!      nprops                   No. material properties
!      lanneal                  Annealing flag - if set to 1, then stresses are zeroed, and state vars should be reset to initial state
!      stepTime                 Elapsed time in this step at end of increment
!      totalTime                Total elapsed time at end of time increment
!      dt                       time step
!      cmname                   Material name
!      coordMp(n,i)             ith coord of nth integration point
!      charLength(i)            Characteristic length of element (in EN234FEA this is (element volume)^1/3)
!      props(k)                 kth material property
!      density                  density
!      strainInc(n,i)           ith strain increment component for nth int pt.  Strain components are in a
!                               basis that rotates with material, stored in order [de11, de22, de33, de12, de23, de13]
!                               NOTE THAT SHEAR STRAIN COMPONENTS ARE NOT DOUBLED IN THE VECTOR
!      relSpinInc(n,i)          ith component of relative spin increment tensor.   The global components of relSpinInc are
!                               dW_ij - dR_ij   where dW is the skew part of displacement gradient increment, and
!                               dR_ij is the rotation increment.
!      tempOld                  Temperature at start of increment
!      stretchOld(n,i)          ith component of right stretch V in co-rotational basis (equivalent to U in global basis) at start of increment
!      defgradOld(n,i)          ith component of deformation gradient in fixed global basis at start of increment
!      fieldOld(n,i)            ith user-defined field variable at nth integration point at start of increment
!      stressOld(n,i)           ith component of Cauchy stress in co-rotational basis (equivalent to R^T sigma R in fixed basis)
!      stateOld(n,i)            ith user defined material state variable at start of increment
!      enerInternOld(n)         specific internal energy at nth integration point, start of increment
!      enerInelasOld(n)         specific irreversible dissipation at nth integration point, start of increment
!      tempNew                  Temperature at end of increment
!      stretchNew(n,i)          ith component of right stretch V in co-rotational basis (equivalent to U in global basis) at start of increment
!      defgradNew(n,i)          ith component of deformation gradient in fixed global basis at end of increment
!      fieldNew(n,i)            ith user-defined field variable at nth integration point at end of increment
!      stressNew(n,i)           ith component of Cauchy stress in co-rotational basis (equivale nt to R^T sigma R in fixed basis)
!      stateNew(n,i)            ith user defined material state variable at end of increment
!      enerInternNew(n)         specific internal energy at nth integration point, end of increment
!      enerInelasNew(n)         specific irreversible dissipation at nth integration point, end of increment
!
!       Coded for 3D problems; small stretch assumption (works approximately for finite rotations)
!=====================================================================================================================
!==============================================Begin Coding The UMAT==================================================       
!========================Material Properties
!       
      E = props(1) ! Young's Modulus (Mpa)
      xnu = props(2) ! Possion's ratio
      Y = props(3) ! Y(Mpa)
      e0 = props(4) ! e0
      m  = props(5) ! m
      edot0 = props(6) ! edot0 (s^-1)
      S = props(7) ! S(Mpa)
      H = props(8) ! H(Mpa)
      td = props(9) ! td(s)
      Omega = props(10) ! Omega
      alfa = props(11) ! alfa

!========================Material Properties
!    ndir = number of e11, e22, e33
!    nsjr = number of e12, e23, e13
      
      ntens = ndir+nshr

!========================Do the loop of integration over different point
! 
      do k = 1,nblock

        eplas = stateOld(k,1)
        ta = stateOld(k,2)
        deplas = stateOld(k,3)
! First Step to determine dedev
        
        ! dedev = devia toric strain tensor, please notice that this is only increment.
        dedev(1:ndir) = strainInc(k,1:ndir)
     1                      - sum(strainInc(k,1:ndir))/3.d0
        ! Specifically record the shear term.
        dedev(ndir+1:ntens) = strainInc(k,ndir+1:ntens)     !   No factor of 2 in dedev
! Particular step
        ! skkstar = stress tensor, step n+1
        skkstar = sum(stressOld(k,1:ndir))
     1          + E*sum(strainInc(k,1:ndir))/(1.d0-2.d0*xnu)
        ! sdevstar = the old devia toric stress tensor, however ,it will update in this procedure
! Second and third step to determine skk_star
        sdevstar(1:ndir) = stressOld(k,1:ndir)
     1                   - sum(stressOld(k,1:ndir))/3.d0
        ! Specifically record the sheare term
        sdevstar(ndir+1:ntens) = stressOld(k,ndir+1:ntens)
        ! update the devia toric strain tensor
        sdevstar(1:ntens) = sdevstar(1:ntens)
     1                    + E*dedev(1:ntens)/(1.d0+xnu)
! Fourth step
        ! dsqrt is the sqrt for the double precision
        ! sestar = stress_e, step n+1
        sestar = dsqrt(1.5d0)*
     1    dsqrt(dot_product(sdevstar(1:ndir),sdevstar(1:ndir)) +
     2    2.d0*dot_product(sdevstar(ndir+1:ntens),
     3                     sdevstar(ndir+1:ntens)) )

        ! Elastic increment (either no stress, or no plastic strain rate)
        ! This is to test whether is plastic deformation or not
        if (sestar/Y<1.d-09.or.edot0==0.d0) then
           stressNew(k,1:ntens) = sdevstar(1:ntens)
           ! Here, it will add the (s11+s22+s33)/3
           stressNew(k,1:ndir) = stressNew(k,1:ndir) + skkstar/3.d0
           ! state(New) remain unchanged 
           stateNew(k,1) = eplas
           stateNew(k,2) = 0.d0
           cycle
        endif

!========================Newton-Rapson Interaction to Obtain the Plastic Properties
!
! err set the initial error
       err = 1.d0
! maxit set the maximum steps allowed
       maxit = 100
! to record the number of steps in the procedure
       nit = 1
! to define the tolerance
       tol = 1.d-6*Y
       if (deplas==0.d0) deplas = 1.d-09/Y
! this procedure is to solve deplas, that is
       dta = 1.d-09/Y
       dta_2 = 1.d-09/Y*1.01
       deplas_2 = 1.d-09/Y
       do while (err>tol)
         
         deplas_2 = deplas*1.00001d0
           
         dta = (omega*dt/deplas-ta+ta*dexp(-deplas/omega))/
     1   (omega/deplas*dexp(-deplas/omega)+1)
         
         dta_2 = (omega*dt/deplas_2-ta+ta*dexp(-deplas_2/omega))/
     1   (omega/deplas_2*dexp(-deplas_2/omega)+1)
         
         ! The first part is to determine the function f
         c1 = Y/S*(1.d0+(eplas+deplas)/e0)**m
         c1_2 = Y/S*(1.d0+(eplas+deplas_2)/e0)**m
         
         c2 = H*(1.d0-dexp(-(((ta+dta)/td)**alfa)))
         c2_2 = H*(1.d0-dexp(-(((ta+dta_2)/td)**alfa)))
         
         c3 = dlog(deplas/dt/edot0)
         c3_2 = dlog(deplas_2/dt/edot0)

         f = sestar/S - 1.5d0*E*deplas/(1.d0+xnu)/S -
     1   (c1 + c2 + c3)
         f_2 = sestar/S - 1.5d0*E*deplas_2/(1.d0+xnu)/S - 
     1   (c1_2 + c2_2 + c3_2)
         
         ! the second part is to do the df/dde by numerical method
         dfde = (f_2-f)/(deplas_2-deplas)
         
! the obtain the new dedplas
         deplas_new = deplas - f/dfde
         if (deplas_new<0.d0) then
            deplas = deplas/10.d0
!            deplas = deplas_new
         else
            deplas = deplas_new
         endif

         nit = nit + 1
         ! dabs is to calculate the new f.
         err = dabs(f)
         if (nit>maxit) then
            write(6,*) ' Newton iterations in UMAT failed to converge '
            
            stop
         endif
         
       end do

!========================Store the new stress
!       
       
        stressNew(k,1:ntens) =
     1           ( 1.d0 - 1.5d0*deplas*E/(1.d0+xnu)/sestar )
     2                           *sdevstar(1:ntens)
        stressNew(k,1:ndir) =
     1              stressNew(k,1:ndir) + skkstar/3.d0

!========================Store the new state variables
!            
        
        stateNew(k,1) = eplas + deplas
        stateNew(k,2) = ta + dta
        stateNew(k,3) = deplas

        write(1,*) stressnew(1,1), ta,deplas
        
      end do


! Best luck~
!
      
      End Subroutine vumat
