!    ABAQUS format UEL subroutine_Written by Siyuan_Song in 31-Oct-2017
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 3D linear elastic continuum element, specifically for large deformation
!
!    The file also contains the following subrouines:
!          abq_UEL_3D_integrationpoints           - defines integration ponits for 3D continuum elements
!          abq_UEL_3D_shapefunctions              - defines shape functions for 3D continuum elements
!          abq_UEL_invert3D                       - computes the inverse and determinant of a 3x3 matrix
!          abq_facenodes_3D                       - returns list of nodes on the face of a 3D element
!
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL_hm7(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3     LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
    !
      INCLUDE 'ABA_PARAM.INC'
    !
    !
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1   SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2   DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3   JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4   PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

    !
    !       Variables that must be computed in this routine
    !       RHS(i)                     Right hand side vector.  In EN234_FEA the dimensions are always RHS(MLVARX,1)
    !       AMATRX(i,j)                Stiffness matrix d RHS(i)/ d DU(j)
    !       SVARS(1:NSVARS)            Element state variables.  Must be updated in this routine
    !       ENERGY(1:8)
    !                                  Energy(1) Kinetic Energy
    !                                  Energy(2) Elastic Strain Energy
    !                                  Energy(3) Creep Dissipation
    !                                  Energy(4) Plastic Dissipation
    !                                  Energy(5) Viscous Dissipation
    !                                  Energy(6) Artificial strain energy
    !                                  Energy(7) Electrostatic energy
    !                                  Energy(8) Incremental work done by loads applied to the element
    !       PNEWDT                     Allows user to control ABAQUS time increments.
    !                                  If PNEWDT<1 then time step is abandoned and computation is restarted with
    !                                  a time increment equal to PNEWDT*DTIME
    !                                  If PNEWDT>1 ABAQUS may increase the time increment by a factor PNEWDT
    !
    !       Variables provided for information
    !       NDOFEL                     Total # DOF for the element
    !       NRHS                       Dimension variable
    !       NSVARS                     Total # element state variables
    !       PROPS(1:NPROPS)            User-specified properties of the element
    !       NPROPS                     No. properties
    !       JPROPS(1:NJPROPS)          Integer valued user specified properties for the element
    !       NJPROPS                    No. integer valued properties
    !       COORDS(i,N)                ith coordinate of Nth node on element
    !       MCRD                       Maximum of (# coords,minimum of (3,#DOF)) on any node
    !       U                          Vector of DOF at the end of the increment
    !       DU                         Vector of DOF increments
    !       V                          Vector of velocities (defined only for implicit dynamics)
    !       A                          Vector of accelerations (defined only for implicit dynamics)
    !       TIME(1:2)                  TIME(1)   Current value of step time
    !                                  TIME(2)   Total time
    !       DTIME                      Time increment
    !       KSTEP                      Current step number (always 1 in EN234_FEA)
    !       KINC                       Increment number
    !       JELEM                      User assigned element number in ABAQUS (internally assigned in EN234_FEA)
    !       PARAMS(1:3)                Time increment parameters alpha, beta, gamma for implicit dynamics
    !       NDLOAD                     Number of user-defined distributed loads defined for this element
    !       JDLTYP(1:NDLOAD)           Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
    !       ADLMAG(1:NDLOAD)           Distributed load magnitudes
    !       DDLMAG(1:NDLOAD)           Increment in distributed load magnitudes
    !       PREDEF(1:2,1:NPREDF,1:NNODE)   Predefined fields.
    !       PREDEF(1,...)              Value of predefined field
    !       PREDEF(2,...)              Increment in predefined field
    !       PREDEF(1:2,1,k)            Value of temperature/temperature increment at kth node
    !       PREDEF(1:2,2:NPREDF,k)     Value of user defined field/field increment at kth node (not used in EN234FEA)
    !       NPREDF                     Number of predefined fields (1 for en234FEA)
    !       LFLAGS                     Control variable
    !       LFLAGS(1)                  Defines procedure type
    !       LFLAGS(2)                  0 => small displacement analysis  1 => Large displacement (NLGEOM option)
    !       LFLAGS(3)                   1 => Subroutine must return both RHS and AMATRX (always true in EN234FEA)
    !                                   2 => Subroutine must return stiffness AMATRX = -dF/du
    !                                   3 => Subroutine must return daming matrix AMATRX = -dF/dudot
    !                                   4 => Subroutine must return mass matrix AMATRX = -dF/duddot
    !                                   5 => Define the RHS only
    !                                   6 => Define the mass matrix for the initial acceleration calculation
    !                                   100 => Define perturbation quantities for output
    !       LFLAGS(4)                   0 => General step   1 => linear perturbation step
    !       LFLAGS(5)                   0 => current approximation to solution based on Newton correction; 1 => based on extrapolation
    !       MLVARX                      Dimension variable (equal to NDOFEL in EN234FEA)
    !       PERIOD                      Time period of the current step
    !
    !
    ! Local Variables
      integer      :: i,j,n_points,kint, nfacenodes, ipoin
      integer      :: face_node_list(8)                       ! List of nodes on an element face
    !
      double precision  ::  xi(3,64)                          ! Volumetric Integration points
      double precision  ::  w(64)                             ! Integration weights
      double precision  ::  N(20)                             ! 3D Shape functions
      double precision  ::  dNdxi(20,3)                       ! 3D Shape function derivatives
      double precision  ::  dxdxi(3,3)                        ! Derivative of position wrt normalized coords
      double precision  ::  dNdx(20,3)                        ! Derivative of shape functions wrt spatial coords
    !
    !   Variables below are for computing integrals over element faces
      double precision  ::  face_coords(3,8)                  ! Coords of nodes on an element face
      double precision  ::  xi2(2,9)                          ! Area integration points
      double precision  ::  N2(9)                             ! 2D shape functions
      double precision  ::  dNdxi2(9,2)                       ! 2D shape function derivatives
      double precision  ::  norm(3)                           ! Normal to an element face
      double precision  ::  dxdxi2(3,2)                       ! Derivative of spatial coord wrt normalized areal coord
    !
      double precision  ::  strain(6)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
      double precision  ::  stress_matrix(3,3)                ! Stress matrix contains nine terms
      double precision  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
      double precision  ::  Mmatrix_stress(3,3)               ! Material stress in matrix form
      double precision  ::  Mvector_stress(6)                 ! Material stress in vector form
      double precision  ::  D(6,6)                            ! stress = D*(strain)  (NOTE FACTOR OF 2 in shear strain)
      double precision  ::  B(6,60)                           ! strain = B*(dof_total)
      double precision  ::  Bstar(9,60)                       ! this is to store the (partial N/ partial x) in lagrangian coordinates
      double precision  ::  F(3,3)                            ! F is to restore deplacement gradient
      double precision  ::  F_inv(3,3)                        ! F_inv is to restore the inverse of displacement gradient
      double precision  ::  F_det                             ! F_det is to restore the determinate of F
      double precision  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
      double precision  ::  bulk_modulus, miu                 ! Material properties
      double precision  ::  G11, G22, G33, G44                ! To store g11, g22, g33, g44
      double precision  ::  G_matrix(6,6)                     ! To store G_matrix of the material properties
      double precision  ::  Cmatrix(3,3)                      ! Is the right Cauchy-Green Deformation tensor as matrix
      double precision  ::  Cmatrix_inv(3,3)                  ! To store the inverse of C as matrix
      double precision  ::  Cmatrix_ave(3,3)                  ! To store the average of C as matrix
      double precision  ::  Cmatrix_det                       ! To store the determinate of C as matrix
      double precision  ::  Cvector(6)                        ! Is the right Cauchy-Green Deformation tensor as vector
      double precision  ::  Cvector_inv(6)                    ! To store the inverse of C as vector
      double precision  ::  Cvector_ave(6)                    ! To store the inverse of C as vector
      double precision  ::  Cvector_det                       ! To store the determinate of C as vector
      double precision  ::  Cvector_star(6)                   ! To store the Cvector with the shear termed doubled
      double precision  ::  Cvector_ave_star(6)               ! To store the Cvector_ave with the shear termed doubled
      double precision  ::  C_I(6)                            ! To store [Cvector_ave_star - I]
      double precision  ::  Q                                 ! To store the temporary term Q
      double precision  ::  P(6)                            ! To store the temporary term P(3,3)
      double precision  ::  H(6,9)                            ! To store the temporary term H(6,9)
      double precision  ::  qmatrix(3,3)                      ! To store the q in matrix form
      double precision  ::  qvector(9)                        ! To store the q in vector form 
      double precision  ::  Omega(6,6)                        ! To store the Omega in matrix form
      double precision  ::  kmatrix(60,60)                   ! To store the little_k in the calculation of geometric stiffness
      double precision  ::  Ymatrix(60,60)                   ! To store the total geometric stiffness
    !
    !     Example ABAQUS UEL implementing 3D linear elastic elements
    !     El props are: miu, bulk_modulus, G11, G22, G33, G44

    !     PROPS(1)         miu
    !     PROPS(2)         bulk_modulus
    !     PROPS(3)         G11
    !     PROPS(4)         G22
    !     PROPS(5)         G33
    !     PROPS(6)         G44
      ! =======================================n_points Calculation===============================
      ! This part is to calculate the number of point in the system
      
      if (NNODE == 4) n_points = 1               ! Linear tet
      if (NNODE == 10) n_points = 4              ! Quadratic tet
      if (NNODE == 8) n_points = 8               ! Linear Hex
      if (NNODE == 20) n_points = 27             ! Quadratic hex

      ! End calculate the number of point
      ! =======================================xi, w, Calculation===============================
      
      call abq_UEL_3D_integrationpoints(n_points, NNODE, xi, w)

      ! =======================================Check the dimensional variable===============================
      if (MLVARX<3*NNODE) then
        write(6,*) ' Error in abaqus UEL '
        write(6,*) ' Variable MLVARX must exceed 3*NNODE'
        write(6,*) ' MLVARX = ',MLVARX,' NNODE = ',NNODE
        stop
      endif

      ! =======================================Initialize R and K matrix===============================
      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
      miu = PROPS(1)
      bulk_modulus = PROPS(2)
      G11 = PROPS(3)
      G22 = PROPS(4)
      G33 = PROPS(5)
      G44 = PROPS(6)
      
      ! The G_matrix should be stored here
      G_matrix(1:6,1:6) = 0.d0
      G_matrix(1,1) = G11
      G_matrix(2,2) = G22
      G_matrix(3,3) = G33
      G_matrix(4,4) = G44
      G_matrix(5,5) = G44
      G_matrix(6,6) = G44
      
      ! =======================================Initialize Energy matrix===============================
      ENERGY(1:8) = 0.d0
      
      ! =================================================================================================
      ! =======================================Main Part to Obtain R and K===============================

    !     --  Loop over integration points
      do kint = 1, n_points
          
        ! ============================To obtain dNdxi, dxdxi, dxidx and dNdx============================== 
        call abq_UEL_3D_shapefunctions(xi(1:3,kint),NNODE,N,dNdxi)
        dxdxi = matmul(coords(1:3,1:NNODE),dNdxi(1:NNODE,1:3))
        call abq_UEL_invert3d(dxdxi,dxidx,determinant)
        dNdx(1:NNODE,1:3) = matmul(dNdxi(1:NNODE,1:3),dxidx)
        
        ! ============================Bstar calculation==============================
        ! This part is to calculate Bstar
        
        Bstar = 0.d0
        
        Bstar(1,1:3*NNODE-2:3) = dNdx(1:NNODE,1)
        Bstar(2,2:3*NNODE-1:3) = dNdx(1:NNODE,2)
        Bstar(3,3:3*NNODE:3)   = dNdx(1:NNODE,3)
        Bstar(4,1:3*NNODE-2:3) = dNdx(1:NNODE,2)
        Bstar(5,2:3*NNODE-1:3) = dNdx(1:NNODE,1)
        Bstar(6,1:3*NNODE-2:3) = dNdx(1:NNODE,3)
        Bstar(7,3:3*NNODE:3)   = dNdx(1:NNODE,1)
        Bstar(8,2:3*NNODE-1:3) = dNdx(1:NNODE,3)
        Bstar(9,3:3*NNODE:3)   = dNdx(1:NNODE,2)       
        
        ! End Calculate Bstar
        ! ===========================F calculation=================================
        ! This part is to calculate F
        
        F(1:3,1:3) = 0.d0
        F_inv(1:3,1:3) = 0.d0
        F_det = 0.d0
        
        do i = 1,3
            do j = 1,3
                F(i,j) = dot_product(dNdx(1:NNODE,j),
     1          U(i:(3*NNODE-3+i):3))
            end do
        end do
        F(1,1) = F(1,1) + 1.d0
        F(2,2) = F(2,2) + 1.d0
        F(3,3) = F(3,3) + 1.d0
        call abq_UEL_invert3d(F(1:3,1:3),F_inv(1:3,1:3),F_det)
        
        ! End Calculate F
        ! ===========================Cauchy Green tensor calculation=================================

        ! Calculate the Right Cauchy Green tensor Cmatrix
        Cmatrix(1:3,1:3) = 0.d0
        Cmatrix_inv(1:3,1:3) = 0.d0
        Cmatrix_ave(1:3,1:3) = 0.d0
        Cmatrix_det = 0.d0
        
        Cmatrix(1:3,1:3) = matmul(transpose(F),F)
        call abq_UEL_invert3d(Cmatrix(1:3,1:3),
     1   Cmatrix_inv(1:3,1:3),Cmatrix_det)
        Cmatrix_ave(1:3,1:3) = Cmatrix(1:3,1:3)/(F_det**(2.d0/3.d0))
        
        ! Calculate the Right Cauchy Green tensor Cvector
        Cvector(1:6) = 0.d0
        Cvector_inv(1:6) = 0.d0
        Cvector_ave(1:6) = 0.d0
        Cvector_det = 0.d0
        
        Cvector(1) = Cmatrix(1,1)
        Cvector(2) = Cmatrix(2,2)
        Cvector(3) = Cmatrix(3,3)
        Cvector(4) = Cmatrix(1,2)
        Cvector(5) = Cmatrix(1,3)
        Cvector(6) = Cmatrix(2,3)
        
        Cvector_inv(1) = Cmatrix_inv(1,1)
        Cvector_inv(2) = Cmatrix_inv(2,2)
        Cvector_inv(3) = Cmatrix_inv(3,3)
        Cvector_inv(4) = Cmatrix_inv(1,2)
        Cvector_inv(5) = Cmatrix_inv(1,3)
        Cvector_inv(6) = Cmatrix_inv(2,3)    
        
        Cvector_ave(1:6) = Cvector(1:6)/(F_det**(2.d0/3.d0))
        
        Cvector_det = Cmatrix_det
        
        ! Calculate the Right Cauchy Green tensor Cvector_star
        Cvector_star(1:6) = 0.d0
        Cvector_ave_star(1:6) = 0.d0
        
        Cvector_star(1) = Cmatrix(1,1)
        Cvector_star(2) = Cmatrix(2,2)
        Cvector_star(3) = Cmatrix(3,3)
        Cvector_star(4) = Cmatrix(1,2)*2.d0
        Cvector_star(5) = Cmatrix(1,3)*2.d0
        Cvector_star(6) = Cmatrix(2,3)*2.d0
        
        Cvector_ave_star(1:6) = Cvector_star(1:6)/(F_det**(2.d0/3.d0))
        
        ! Calculate a special term C_I which calculate Cvector_ave_star - I
        C_I(1:6) = 0.d0
        C_I(1:6) = Cvector_ave_star(1:6)
        C_I(1:3) = C_I(1:3) - 1.d0 

        ! End calculate Right Cauchy Green tensor
        ! ===========================P and Q calculation=================================
        
        ! Calculate the temporary term P
        P = 0.d0
        P = 1.d0/2.d0/(F_det**(2.d0/3.d0))*
     1  (matmul(G_matrix,C_I)-1.d0/3.d0*
     2  dot_product(Cvector_star,matmul(G_matrix,C_I))*Cvector_inv) 
        
        ! Calculate the temporary term Q
        Q = 0.d0
        Q = 1.d0/4.d0*dot_product(C_I,matmul(G_matrix,C_I))

        ! End calculate P and Q
        ! ===========================Material Stress calculation=================================
        
        ! Calculate the material stress
        Mvector_stress = 0.d0
        Mvector_stress = miu*exp(Q)*P+
     1  bulk_modulus*F_det*(F_det-1)*Cvector_inv
        
        Mmatrix_stress = 0.d0
        Mmatrix_stress(1,1) = Mvector_stress(1)
        Mmatrix_stress(2,2) = Mvector_stress(2)
        Mmatrix_stress(3,3) = Mvector_stress(3)
        Mmatrix_stress(1,2) = Mvector_stress(4)
        Mmatrix_stress(1,3) = Mvector_stress(5)
        Mmatrix_stress(2,3) = Mvector_stress(6)
        Mmatrix_stress(2,1) = Mvector_stress(4)
        Mmatrix_stress(3,1) = Mvector_stress(5)
        Mmatrix_stress(3,2) = Mvector_stress(6)

        ! End calculate Material Stress
        ! ===========================qmatrix and q vector calculation=================================
        
        ! Calculate the qmatrix and qvector
        qmatrix = 0.d0
        qvector = 0.d0
        
        qmatrix = matmul(Mmatrix_stress,transpose(F))
        qvector(1) = qmatrix(1,1)
        qvector(2) = qmatrix(2,2)
        qvector(3) = qmatrix(3,3)
        qvector(4) = qmatrix(2,1)
        qvector(5) = qmatrix(1,2)
        qvector(6) = qmatrix(3,1)
        qvector(7) = qmatrix(1,3)
        qvector(8) = qmatrix(3,2)
        qvector(9) = qmatrix(2,3)

        ! End Calculate qvector
        ! ===========================Omega calculation=================================        
        ! Calculate the Omega
        Omega = 0.d0
        Omega(1,1) = Cmatrix_inv(1,1)*Cmatrix_inv(1,1)
        Omega(1,2) = Cmatrix_inv(1,2)*Cmatrix_inv(1,2)
        Omega(1,3) = Cmatrix_inv(1,3)*Cmatrix_inv(1,3)
        Omega(1,4) = Cmatrix_inv(1,1)*Cmatrix_inv(1,2)
        Omega(1,5) = Cmatrix_inv(1,1)*Cmatrix_inv(1,3)
        Omega(1,6) = Cmatrix_inv(1,2)*Cmatrix_inv(1,3)
        Omega(2,2) = Cmatrix_inv(2,2)*Cmatrix_inv(2,2)
        Omega(2,3) = Cmatrix_inv(2,3)*Cmatrix_inv(2,3)
        Omega(2,4) = Cmatrix_inv(2,1)*Cmatrix_inv(2,2)
        Omega(2,5) = Cmatrix_inv(2,1)*Cmatrix_inv(2,3)
        Omega(2,6) = Cmatrix_inv(2,2)*Cmatrix_inv(2,3)
        Omega(3,3) = Cmatrix_inv(3,3)*Cmatrix_inv(3,3)
        Omega(3,4) = Cmatrix_inv(3,1)*Cmatrix_inv(3,2)
        Omega(3,5) = Cmatrix_inv(3,1)*Cmatrix_inv(3,3)
        Omega(3,6) = Cmatrix_inv(3,2)*Cmatrix_inv(3,3)
        Omega(4,4) = (Cmatrix_inv(1,1)*Cmatrix_inv(2,2)+
     1  Cmatrix_inv(1,2)*Cmatrix_inv(1,2))/2.d0
        Omega(4,5) = (Cmatrix_inv(1,1)*Cmatrix_inv(2,3)+
     1  Cmatrix_inv(1,3)*Cmatrix_inv(1,2))/2.d0
        Omega(4,6) = (Cmatrix_inv(1,2)*Cmatrix_inv(2,3)+
     1  Cmatrix_inv(1,3)*Cmatrix_inv(2,2))/2.d0
        Omega(5,5) = (Cmatrix_inv(1,1)*Cmatrix_inv(3,3)+
     1  Cmatrix_inv(1,3)*Cmatrix_inv(1,3))/2.d0
        Omega(5,6) = (Cmatrix_inv(1,2)*Cmatrix_inv(3,3)+
     1  Cmatrix_inv(1,3)*Cmatrix_inv(2,3))/2.d0
        Omega(6,6) = (Cmatrix_inv(2,2)*Cmatrix_inv(3,3)+
     1  Cmatrix_inv(2,3)*Cmatrix_inv(2,3))/2.d0
        ! Symmetry Property
        do i=2,6
            do j=1,(i-1)
                Omega(i,j)=Omega(j,i)
            end do
        end do    
        
        ! End Calculate Omega 
        ! ===========================D_matrix calculation=================================   
        ! Calculate the D_matrix
        D = miu*exp(Q)*1.d0/(F_det**(4.d0/3.d0))*(
     1  G_matrix-1.d0/3.d0*(
     2  matmul(G_matrix,(spread(Cvector_star,dim=2,ncopies=6)*
     3  spread(Cvector_inv,dim=1,ncopies=6)))+
     4  spread(Cvector_inv,dim=2,ncopies=6)*
     5  spread(matmul(G_matrix,Cvector_star),dim=1,ncopies=6))
     6  -1.d0/3.d0*(F_det**(2.d0/3.d0))*
     7  dot_product(Cvector_star,matmul(G_matrix,C_I))*Omega+
     8  1.d0/9.d0*dot_product(Cvector_star,matmul(G_matrix,Cvector_star
     9  ))*spread(Cvector_inv,dim=2,ncopies=6)*
     1  spread(Cvector_inv,dim=1,ncopies=6))
        D = D+miu*exp(Q)*(2.d0*spread(P,dim=2,ncopies=6)*
     1  spread((P-1.d0/3.d0*Cvector_inv),dim=1,ncopies=6)-
     2  1.d0/3.d0/(F_det**(2.d0/3.d0))*
     3  spread(Cvector_inv,dim=2,ncopies=6)*
     4  spread((matmul(G_matrix,C_I)),dim=1,ncopies=6))
        D = D+bulk_modulus*F_det*(
     1  (2.d0*F_det-1.d0)*spread(Cvector_inv,dim=2,ncopies=6)*
     2  spread(Cvector_inv,dim=1,ncopies=6)+2.d0*(F_det-1.d0)*Omega)
        
        ! End Calculate the D_matrix
        ! ===========================D_matrix calculation=================================   
        ! Calculate H matrix
        H = 0.d0
        H(1,1) = F(1,1)
        H(1,5) = F(2,1)
        H(1,7) = F(3,1)
        H(2,2) = F(2,2)
        H(2,4) = F(1,2)
        H(2,9) = F(3,2)
        H(3,3) = F(3,3)
        H(3,6) = F(1,3)
        H(3,8) = F(2,3)
        H(4,1) = F(1,2)
        H(4,2) = F(2,1)
        H(4,4) = F(1,1)
        H(4,5) = F(2,2)
        H(4,7) = F(3,2)
        H(4,9) = F(3,1)
        H(5,1) = F(1,3)
        H(5,3) = F(3,1)
        H(5,5) = F(2,3)
        H(5,6) = F(1,1)
        H(5,7) = F(3,3)
        H(5,8) = F(2,1)
        H(6,2) = F(2,3)
        H(6,3) = F(3,2)
        H(6,4) = F(1,3)
        H(6,6) = F(1,2)
        H(6,8) = F(2,2)
        H(6,9) = F(3,3)
        
        ! End calculate H matrix
        ! ===========================Y_matrix calculation=================================
        ! Calculate kmatrix first
        kmatrix = 0.d0
        kmatrix(1:NNODE,1:NNODE) = matmul(matmul(dNdx(1:NNODE,1:3),
     1  Mmatrix_stress(1:3,1:3)),transpose(dNdx(1:NNODE,1:3)))
        
        ! Calculate the Ymatrix as following
        Ymatrix = 0.d0
        do i = 1,NNODE
            do j = 1,NNODE
                Ymatrix(3*i-2,3*j-2) = kmatrix(i,j)
                Ymatrix(3*i-1,3*j-1) = kmatrix(i,j)
                Ymatrix(3*i-0,3*j-0) = kmatrix(i,j)
            end do
        end do
        
        ! End calculate Ymatrix 
        ! ===========================Residual Calculation=================================
        ! Calculate the Residual
        RHS(1:3*NNODE,1) = RHS(1:3*NNODE,1)-
     1   matmul(transpose(Bstar(1:9,1:3*NNODE)),qvector(1:9))*
     2  w(kint)*determinant
        
        ! End Calculate the Residual
        ! ===========================Stiffness Calculation=================================    
        !Calculate the Stiffness
        AMATRX(1:3*NNODE,1:3*NNODE) = AMATRX(1:3*NNODE,1:3*NNODE)
     1  + (matmul(matmul(matmul(matmul(transpose(Bstar(1:9,1:3*NNODE))
     2  ,transpose(H)),D),H),Bstar(1:9,1:3*NNODE))+
     3  Ymatrix(1:3*NNODE,1:3*NNODE))*w(kint)*determinant
        
        ! End Calculate the Stiffness
        ! ===========================Cauchy stress calculation=================================
        ! This part is to calculate Cauchy stress
        stress_matrix = 0.d0
        stress = 0.d0
        stress_matrix = matmul(matmul(F,Mmatrix_stress),
     1  transpose(F))/F_det
        
        stress(1) = stress_matrix(1,1)
        stress(2) = stress_matrix(2,2)
        stress(3) = stress_matrix(3,3)
        stress(4) = stress_matrix(1,2)
        stress(5) = stress_matrix(1,3)
        stress(6) = stress_matrix(2,3)
        
        ! End calculate Cauchy stress
        ! ===========================Energy Calculation calculation=================================        
        ENERGY(2) = ENERGY(2)
     1  +(miu/2.d0*(exp(Q)-1)+bulk_modulus/2.d0*((F_det-1.d0)**2))
     2  *w(kint)*determinant   ! Store the elastic strain energy        
        
        ! ===========================Stress Storing================================= 
        if (NSVARS>=n_points*6) then   ! Store stress at each integration point (if space was allocated to do so)
            SVARS(6*kint-5:6*kint) = stress(1:6)
        endif
      end do
      
    ! ===========================All parts finished=================================
    ! ===========================I hope every thing runs well=================================

      PNEWDT = 1.d0          ! This leaves the timestep unchanged (ABAQUS will use its own algorithm to determine DTIME)
    !
    !   Apply distributed loads
    !
    !   Distributed loads are specified in the input file using the Un option in the input file.
    !   n specifies the face number, following the ABAQUS convention
    !
    !
      do j = 1,NDLOAD

        call abq_facenodes_3D(NNODE,iabs(JDLTYP(j,1)),
     1                                     face_node_list,nfacenodes)

        do i = 1,nfacenodes
            face_coords(1:3,i) = coords(1:3,face_node_list(i))
        end do

        if (nfacenodes == 3) n_points = 3
        if (nfacenodes == 6) n_points = 4
        if (nfacenodes == 4) n_points = 4
        if (nfacenodes == 8) n_points = 9

        call abq_UEL_2D_integrationpoints(n_points, nfacenodes, xi2, w)

        do kint = 1,n_points
            call abq_UEL_2D_shapefunctions(xi2(1:2,kint),
     1                        nfacenodes,N2,dNdxi2)
            dxdxi2 = matmul(face_coords(1:3,1:nfacenodes),
     1                           dNdxi2(1:nfacenodes,1:2))
            norm(1)=(dxdxi2(2,1)*dxdxi2(3,2))-(dxdxi2(2,2)*dxdxi2(3,1))
            norm(2)=(dxdxi2(1,1)*dxdxi2(3,2))-(dxdxi2(1,2)*dxdxi2(3,1))
            norm(3)=(dxdxi2(1,1)*dxdxi2(2,2))-(dxdxi2(1,2)*dxdxi2(2,1))

            do i = 1,nfacenodes
                ipoin = 3*face_node_list(i)-2
                RHS(ipoin:ipoin+2,1) = RHS(ipoin:ipoin+2,1)
     1                 - N2(1:nfacenodes)*adlmag(j,1)*norm(1:3)*w(kint)! Note determinant is already in normal
            end do
        end do
      end do

      return

      END SUBROUTINE UEL_hm7

      subroutine abq_UEL_3D_integrationpoints(n_points, n_nodes, xi, w)

      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes

      double precision, intent(out) :: xi(3,*)
      double precision, intent(out) :: w(*)

      integer :: i,j,k,n

      double precision x1D(4), w1D(4)

    !         Defines integration points and weights for 3D continuum elements

      if (n_nodes  == 4.or.n_nodes ==10 ) then   ! Tetrahedral elements
        if (n_points == 1) then
            xi(1,1) = 0.25D0
            xi(2,1) = 0.25D0
            xi(3,1) = 0.25D0
            w(1) = 1.D0/6.D0
        else if (n_points == 4) then
            xi(1,1) = 0.58541020
            xi(2,1) = 0.13819660
            xi(3,1) = xi(2,1)
            xi(1,2) = xi(2,1)
            xi(2,2) = xi(1,1)
            xi(3,2) = xi(2,1)
            xi(1,3) = xi(2,1)
            xi(2,3) = xi(2,1)
            xi(3,3) = xi(1,1)
            xi(1,4) = xi(2,1)
            xi(2,4) = xi(2,1)
            xi(3,4) = xi(2,1)
            w(1:4) = 1.D0/24.D0
        else if (n_points == 5) then
            xi(1,1) = 0.25d0
            xi(2,1) = 0.25d0
            xi(3,1) = 0.25d0
            xi(1,2) = 0.5d0
            xi(2,2) = 1.d0/6.d0
            xi(3,2) = 1.d0/6.d0
            xi(1,3) = 1.d0/6.d0
            xi(2,3) = 0.5d0
            xi(3,3) = 1.d0/6.d0
            xi(1,4) = 1.d0/6.d0
            xi(2,4) = 1.d0/6.d0
            xi(3,4) = 0.5d0
            xi(1,5) = 1.d0/6.d0
            xi(2,5) = 1.d0/6.d0
            xi(3,5) = 1.d0/6.d0
            w(1) = -4.d0/30.d0
            w(2:5) = 3.d0/40.d0
        else
            write(6,*) 'Incorrect # of int pts for tetrahedral element '
            write(6, *) ' called with ',n_points
            stop
        endif
      else if ( n_nodes == 8 .or. n_nodes == 20 ) then   ! 8 or 20 noded hexahedral elements
        if (n_points == 1) then
            xi(1,1) = 0.D0
            xi(2,1) = 0.D0
            xi(3,1) = 0.D0
            w(1) = 8.D0
        else if (n_points == 8) then
            x1D(1) = -0.5773502692
            x1D(2) =  0.5773502692
            do k = 1,2
                do j = 1,2
                    do i = 1,2
                        n = 4*(k-1) + 2*(j-1) + i
                        xi(1,n) = x1D(i)
                        xi(2,n) = x1D(j)
                        xi(3,n) = x1D(k)
                    end do
                end do
            end do
            w(1:8) = 1.D0
        else if (n_points == 27) then
            x1D(1) = -0.7745966692
            x1D(2) = 0.
            x1D(3) = 0.7745966692
            w1D(1) = 0.5555555555D0
            w1D(2) = 0.888888888D0
            w1D(3) = 0.55555555555D0
            do k = 1,3
                do j = 1,3
                    do i = 1,3
                        n = 9*(k-1) + 3*(j-1) + i
                        xi(1,n) = x1D(i)
                        xi(2,n) = x1D(j)
                        xi(3,n) = x1D(k)
                        w(n) = w1D(i)*w1D(j)*w1D(k)
                    end do
                end do
            end do
        else if (n_points == 64) then
            x1D(1) = .8611363115940526D+00
            x1D(2) = .3399810435848563D+00
            x1D(3) = -.3399810435848563D+00
            x1D(4) = -.8611363115940526D+00
            w1D(1) = .3478548451374538D+00
            w1D(2) = .6521451548625461D+00
            w1D(3) = .6521451548625461D+00
            w1D(4) = .3478548451374538D+00
            do k = 1,4
                do j = 1,4
                    do i = 1,4
                        n = 16*(k-1) + 4*(j-1) + i
                        xi(1,n) = x1D(i)
                        xi(2,n) = x1D(j)
                        xi(3,n) = x1D(k)
                        w(n) = w1D(i)*w1D(j)*w1D(k)
                    end do
                end do
            end do
        endif
      endif

      return

      end subroutine abq_UEL_3D_integrationpoints

      subroutine abq_UEL_3D_shapefunctions(xi,n_nodes,f,df)

      implicit none
      integer, intent(in) :: n_nodes

      double precision, intent(in) :: xi(3)
      double precision, intent(out) :: f(20)
      double precision, intent(out) :: df(20,3)
      double precision xi4

!   Defines shape functions for 3D continuum elements

      if (n_nodes == 4) then
        f(1) = xi(1)
        f(2) = xi(2)
        f(3) = xi(3)
        f(4) = 1.-xi(1)-xi(2)-xi(3)
        df(1,1) = 1.
        df(2,2) = 1.
        df(3,3) = 1.
        df(4,1) = -1.
        df(4,2) = -1.
        df(4,3) = -1.
      else if (n_nodes == 10) then
        xi4 = 1.D0-xi(1)-xi(2)-xi(3)
        f(1) = (2.*xi(1)-1.)*xi(1)
        f(2) = (2.*xi(2)-1.)*xi(2)
        f(3) = (2.*xi(3)-1.)*xi(3)
        f(4) = (2.*xi4-1.)*xi4
        f(5) = 4.*xi(1)*xi(2)
        f(6) = 4.*xi(2)*xi(3)
        f(7) = 4.*xi(3)*xi(1)
        f(8) = 4.*xi(1)*xi4
        f(9) = 4.*xi(2)*xi4
        f(10) = 4.*xi(3)*xi4
        df(1,1) = (4.*xi(1)-1.)
        df(2,2) = (4.*xi(2)-1.)
        df(3,3) = (4.*xi(3)-1.)
        df(4,1) = -(4.*xi4-1.)
        df(4,2) = -(4.*xi4-1.)
        df(4,3) = -(4.*xi4-1.)
        df(5,1) = 4.*xi(2)
        df(5,2) = 4.*xi(1)
        df(6,2) = 4.*xi(3)
        df(6,3) = 4.*xi(2)
        df(7,1) = 4.*xi(3)
        df(7,3) = 4.*xi(1)
        df(8,1) = 4.*(xi4-xi(1))
        df(8,2) = -4.*xi(1)
        df(8,3) = -4.*xi(1)
        df(9,1) = -4.*xi(2)
        df(9,2) = 4.*(xi4-xi(2))
        df(9,3) = -4.*xi(2)
        df(10,1) = -4.*xi(3)*xi4
        df(10,2) = -4.*xi(3)
        df(10,3) = 4.*(xi4-xi(3))
      else if (n_nodes == 8) then
        f(1) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3))/8.
        f(2) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3))/8.
        f(3) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3))/8.
        f(4) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3))/8.
        f(5) = (1.-xi(1))*(1.-xi(2))*(1.+xi(3))/8.
        f(6) = (1.+xi(1))*(1.-xi(2))*(1.+xi(3))/8.
        f(7) = (1.+xi(1))*(1.+xi(2))*(1.+xi(3))/8.
        f(8) = (1.-xi(1))*(1.+xi(2))*(1.+xi(3))/8.
        df(1,1) = -(1.-xi(2))*(1.-xi(3))/8.
        df(1,2) = -(1.-xi(1))*(1.-xi(3))/8.
        df(1,3) = -(1.-xi(1))*(1.-xi(2))/8.
        df(2,1) = (1.-xi(2))*(1.-xi(3))/8.
        df(2,2) = -(1.+xi(1))*(1.-xi(3))/8.
        df(2,3) = -(1.+xi(1))*(1.-xi(2))/8.
        df(3,1) = (1.+xi(2))*(1.-xi(3))/8.
        df(3,2) = (1.+xi(1))*(1.-xi(3))/8.
        df(3,3) = -(1.+xi(1))*(1.+xi(2))/8.
        df(4,1) = -(1.+xi(2))*(1.-xi(3))/8.
        df(4,2) = (1.-xi(1))*(1.-xi(3))/8.
        df(4,3) = -(1.-xi(1))*(1.+xi(2))/8.
        df(5,1) = -(1.-xi(2))*(1.+xi(3))/8.
        df(5,2) = -(1.-xi(1))*(1.+xi(3))/8.
        df(5,3) = (1.-xi(1))*(1.-xi(2))/8.
        df(6,1) = (1.-xi(2))*(1.+xi(3))/8.
        df(6,2) = -(1.+xi(1))*(1.+xi(3))/8.
        df(6,3) = (1.+xi(1))*(1.-xi(2))/8.
        df(7,1) = (1.+xi(2))*(1.+xi(3))/8.
        df(7,2) = (1.+xi(1))*(1.+xi(3))/8.
        df(7,3) = (1.+xi(1))*(1.+xi(2))/8.
        df(8,1) = -(1.+xi(2))*(1.+xi(3))/8.
        df(8,2) = (1.-xi(1))*(1.+xi(3))/8.
        df(8,3) = (1.-xi(1))*(1.+xi(2))/8.
      else if (n_nodes == 20) then
        f(1)=(1.-xi(1))*(1.-xi(2))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)/8.
        f(2)=(1.+xi(1))*(1.-xi(2))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)/8.
        f(3)=(1.+xi(1))*(1.+xi(2))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)/8.
        f(4)=(1.-xi(1))*(1.+xi(2))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)/8.
        f(5)=(1.-xi(1))*(1.-xi(2))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)/8.
        f(6)=(1.+xi(1))*(1.-xi(2))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)/8.
        f(7)=(1.+xi(1))*(1.+xi(2))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)/8.
        f(8)=(1.-xi(1))*(1.+xi(2))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)/8.
        f(9) = (1.-xi(1)**2.)*(1.-xi(2))*(1.-xi(3))/4.
        f(10) = (1.+xi(1))*(1.-xi(2)**2.)*(1.-xi(3))/4.
        f(11) = (1.-xi(1)**2.)*(1.+xi(2))*(1.-xi(3))/4.
        f(12) = (1.-xi(1))*(1.-xi(2)**2.)*(1.-xi(3))/4.
        f(13) = (1.-xi(1)**2.)*(1.-xi(2))*(1.+xi(3))/4.
        f(14) = (1.+xi(1))*(1.-xi(2)**2.)*(1.+xi(3))/4.
        f(15) = (1.-xi(1)**2.)*(1.+xi(2))*(1.+xi(3))/4.
        f(16) = (1.-xi(1))*(1.-xi(2)**2.)*(1.+xi(3))/4.
        f(17) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3)**2.)/4.
        f(18) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3)**2.)/4.
        f(19) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3)**2.)/4.
        f(20) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3)**2.)/4.
        df(1,1) = (-(1.-xi(2))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
        df(1,2) = (-(1.-xi(1))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
        df(1,3) = (-(1.-xi(1))*(1.-xi(2))*(-xi(1)-xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.

        df(2,1) = ((1.-xi(2))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)
     1           +(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
        df(2,2) = (-(1.+xi(1))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)
     1          -(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
        df(2,3) = (-(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-xi(3)-2.)
     1           -(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.

        df(3,1) = ((1.+xi(2))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(3,2) = ((1.+xi(1))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(3,3) = (-(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-xi(3)-2.)
     1           -(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.

        df(4,1) = (-(1.+xi(2))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(4,2) = ((1.-xi(1))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)
     1            +(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(4,3) = (-(1.-xi(1))*(1.+xi(2))*(-xi(1)+xi(2)-xi(3)-2.)
     1           -(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
        df(5,1) = (-(1.-xi(2))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(5,2) = (-(1.-xi(1))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)
     1           -(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(5,3) = ((1.-xi(1))*(1.-xi(2))*(-xi(1)-xi(2)+xi(3)-2.)
     1           +(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(6,1) = ((1.-xi(2))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(6,2) = (-(1.+xi(1))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)
     1           -(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(6,3) = ((1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
        df(7,1) = ((1.+xi(2))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(7,2) = ((1.+xi(1))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(7,3) = ((1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)+xi(3)-2.)
     1           +(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(8,1) = (-(1.+xi(2))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)
     1           -(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(8,2) = ((1.-xi(1))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)
     1           +(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(8,3) = ((1.-xi(1))*(1.+xi(2))*(-xi(1)+xi(2)+xi(3)-2.)
     1           +(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
        df(9,1)  = -2.*xi(1)*(1.-xi(2))*(1.-xi(3))/4.
        df(9,2)  = -(1.-xi(1)**2.)*(1.-xi(3))/4.
        df(9,3)  = -(1.-xi(1)**2.)*(1.-xi(2))/4.
        df(10,1)  = (1.-xi(2)**2.)*(1.-xi(3))/4.
        df(10,2)  = -2.*xi(2)*(1.+xi(1))*(1.-xi(3))/4.
        df(10,3)  = -(1.-xi(2)**2.)*(1.+xi(1))/4.
        df(11,1)  = -2.*xi(1)*(1.-xi(2))*(1.-xi(3))/4.
        df(11,2)  = -(1.-xi(1)**2.)*(1.-xi(3))/4.
        df(11,3)  = -(1.-xi(1)**2.)*(1.-xi(2))/4.
        df(12,1)  = -(1.-xi(2)**2.)*(1.-xi(3))/4.
        df(12,2)  = -2.*xi(2)*(1.-xi(1))*(1.-xi(3))/4.
        df(12,3)  = -(1.-xi(2)**2.)*(1.-xi(1))/4.
        df(13,1)  = -2.*xi(1)*(1.-xi(2))*(1.+xi(3))/4.
        df(13,2)  = -(1.-xi(1)**2.)*(1.+xi(3))/4.
        df(13,3)  = (1.-xi(1)**2.)*(1.-xi(2))/4.
        df(14,1)  = (1.-xi(2)**2.)*(1.+xi(3))/4.
        df(14,2)  = -2.*xi(2)*(1.+xi(1))*(1.+xi(3))/4.
        df(14,3)  = (1.-xi(2)**2.)*(1.+xi(1))/4.
        df(15,1)  = 2.*xi(1)*(1.+xi(2))*(1.+xi(3))/4.
        df(15,2)  = (1.-xi(1)**2.)*(1.+xi(3))/4.
        df(15,3)  = (1.-xi(1)**2.)*(1.+xi(2))/4.
        df(16,1)  = -(1.-xi(2)**2.)*(1.+xi(3))/4.
        df(16,2)  = -2.*xi(2)*(1.-xi(1))*(1.+xi(3))/4.
        df(16,3)  = (1.-xi(2)**2.)*(1.-xi(1))/4.
        df(17,1) = -(1.-xi(2))*(1.-xi(3)**2.)/4.
        df(17,2) = -(1.-xi(1))*(1.-xi(3)**2.)/4.
        df(17,3) = -xi(3)*(1.-xi(1))*(1.-xi(2))/2.
        df(18,1) = (1.-xi(2))*(1.-xi(3)**2.)/4.
        df(18,2) = -(1.+xi(1))*(1.-xi(3)**2.)/4.
        df(18,3) = -xi(3)*(1.+xi(1))*(1.-xi(2))/2.
        df(19,1) = (1.+xi(2))*(1.-xi(3)**2.)/4.
        df(19,2) = (1.+xi(1))*(1.-xi(3)**2.)/4.
        df(19,3) = -xi(3)*(1.+xi(1))*(1.+xi(2))/2.
        df(20,1) = -(1.+xi(2))*(1.-xi(3)**2.)/4.
        df(20,2) = (1.-xi(1))*(1.-xi(3)**2.)/4.
        df(20,3) = -xi(3)*(1.-xi(1))*(1.+xi(2))/2.
      endif


      end subroutine abq_UEL_3D_shapefunctions

      subroutine abq_UEL_invert3d(A,A_inverse,determinant)

      double precision, intent(in) :: A(3,3)
      double precision, intent(out) :: A_inverse(3,3)
      double precision, intent(out) :: determinant

      double precision COFACTOR(3,3)

!   Compute inverse and determinant of 3x3 matrix

      determinant =   A(1,1)*A(2,2)*A(3,3)
     1   - A(1,1)*A(2,3)*A(3,2)
     2   - A(1,2)*A(2,1)*A(3,3)
     3   + A(1,2)*A(2,3)*A(3,1)
     4   + A(1,3)*A(2,1)*A(3,2)
     5   - A(1,3)*A(2,2)*A(3,1)

      IF (determinant==0.d0) THEN
        write(6,*) ' Error in subroutine abq_UEL_inver3d'
        write(6,*) ' A 3x3 matrix has a zero determinant'
        stop
      endif
      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      A_inverse = transpose(COFACTOR) / determinant


      end subroutine abq_UEL_invert3d

      subroutine abq_facenodes_3D(nelnodes,face,list,nfacenodes)

      implicit none

      integer, intent (in)      :: nelnodes
      integer, intent (in)      :: face
      integer, intent (out)     :: list(*)
      integer, intent (out)     :: nfacenodes

    !
    !        Subroutine to return list of nodes on an element face for standard 3D solid elements
    !

      if (nelnodes == 4) then
        nfacenodes = 3
        if   (face == 1) list(1:3) = [1,2,3]
        if (face == 2) list(1:3) = [1,4,2]
        if (face == 3) list(1:3) = [2,4,3]
        if (face == 4) list(1:3) = [3,4,1]
      else if (nelnodes ==6) then
        nfacenodes = 3
        if (face==1) list(1:3) = [1,2,3]
        if (face==2) list(1:3) = [6,5,4]
        if (face==3) list(1:4) = [1,2,5,4]
        if (face==4) list(1:4) = [2,3,6,5]
        if (face==5) list(1:4) = [4,6,3,1]
        if (face>2) nfacenodes = 4
      else if (nelnodes == 10) then
        nfacenodes = 6
        if   (face == 1) list(1:6) = [1,2,3,5,6,7]
        if (face == 2) list(1:6) = [1,4,2,8,9,5]
        if (face == 3) list(1:6) = [2,4,3,9,10,6]
        if (face == 4) list(1:6) = [3,4,1,10,8,7]
      else if (nelnodes == 8) then
        nfacenodes = 4
        if (face==1) list(1:4) = [1,2,3,4]
        if (face==2) list(1:4) = [5,8,7,6]
        if (face==3) list(1:4) = [1,5,6,2]
        if (face==4) list(1:4) = [2,6,7,3]
        if (face==5) list(1:4) = [3,7,8,4]
        if (face==6) list(1:4) = [4,8,5,1]
      else if (nelnodes ==15) then
        nfacenodes = 6
        if (face==1) list(1:6) = [1,2,3,7,8,9]
        if (face==2) list(1:6) = [6,5,4,11,10,12]
        if (face==3) list(1:8) = [1,2,5,4,7,14,10,13]
        if (face==4) list(1:8) = [2,3,6,5,8,15,11,14]
        if (face==5) list(1:8) = [4,6,3,1,12,15,9,13]
        if (face>2) nfacenodes = 8
      else  if (nelnodes == 20) then
        nfacenodes = 8
        if (face == 1) list(1:8) = [1,2,3,4,9,10,11,12]
        if (face == 2) list(1:8) = [5,8,7,6,16,15,14,13]
        if (face == 3) list(1:8) = [1,5,6,2,17,13,18,9]
        if (face == 4) list(1:8) = [2,6,7,3,18,14,19,10]
        if (face == 5) list(1:8) = [3,7,8,4,19,15,6,11]
        if (face == 6) list(1:8) = [4,8,5,1,20,16,17,12]
      endif

      end subroutine abq_facenodes_3D