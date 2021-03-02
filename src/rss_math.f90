!---------------------------------------------------------------------------
!!
!>  @anchor     rss_functions
!!
!!  @brief      Math helper subroutines and functions
!!
!!  @author     Christopher Kebschull
!!
!!  @date       <ul>
!!                <li>24.08.2015 (initial implementation in SMART)</li>
!!                <li>29.02.2016 (using libslam now)</li>
!!                <li>11.07.2017 (migrating to librss now)</li>
!!              </ul>
!!
!---------------------------------------------------------------------------
module rss_math
    use slam_error_handling,    only: isControlled,hasToReturn,hasFailed,       &
                                        checkIn,checkOut,                       &
                                        setError,WARNING,FATAL,E_SPECIAL,       &
                                        E_MATRIX_NOT_POSITIVE_DEFINITE,         &
                                        E_UNKNOWN_PARAMETER
    use slam_strings,           only: toString
    use rss_types,              only: dp

    contains

    !===================================================
    !
    !> @anchor    argpi
    !!
    !! @brief     Subroutine for determining the angle between -PI (-180) and 
    !!              +PI (180)
    !!
    !! @author    Joerg Boettcher (03/86)
    !! @author    Bernhard Schlarmann (07/92)
    !! @aithor    Christopher Kebschull
    !!
    !! @date      <ul>
    !!              <li> 11.07.2017 (migrated to librss)</li>
    !!            </ul>
    !!
    !!---------------------------------------------------------------------
    real(dp) function argpi(x)
        use slam_types,     only: dp
        use slam_math,      only: pi,twopi

        character(len=*), parameter                 :: csubid = "argpi"         ! Function id
        real(dp),intent(in) :: x
        real(dp)            :: xvalue

        if(isControlled()) then
            if(hasToReturn()) return
            call checkIn(csubid)
        end if

        xvalue = dmod (x,twopi)
        if (dabs(xvalue).gt.pi) then
            xvalue = xvalue - dsign(twopi,x)
        end if

        argpi = xvalue

        if(isControlled()) then
            call checkOut(csubid)
        end if

    end function argpi
    
    !!--------------------------------------------------------------------------
    !!
    !> @anchor    unscented_transformation
    !!
    !> @brief     Subroutine for sigma point decomposition
    !!
    !> @author    Eduard Gamper
    !!
    !> @date      <ul>
    !!              <li> 08.01.2014 (initial design)</li>
    !!            </ul>
    !!
    !> @details   Based on equations from "Kalman Filtering and Neural Networks"
    !!              p. 232 table 7.2 The sigma points of the state vectors are
    !!              dervied.
    !!
    !> @param[in] weighting - Kind of weighting (1->Haykin, 2->Wan, 3->Julier)
    !!
    !< @param[in] alpha - Spread of sigma points around state (1>=alpha>=0.0001)
    !!
    !> @param[in] ephemeris - contains state vector and covariance
    !!
    !> @param[out] x_chi - Sigma point matrix
    !!
    !> @result lsuccess - .true. when the process was successful
    !!
    !!--------------------------------------------------------------------------
    ! function unscented_transformation(weighting,alpha,ephemeris,x_chi) result(lsuccess)
    !     use slam_linAlgebra,    only: invertMatrix
    !     use rss_types,          only: t_ephemeris,dp

    !     integer,parameter                           :: n = 6                    ! State vector elements

    !     !Input
    !     integer,intent(in)                          :: weighting                ! Kind of weighting (1->Haykin, 2->Wan, 3->Julier)
    !     real(dp),intent(in)                         :: alpha                    ! Spread of sigma points around state (1>=alpha>=0.0001)
    !     class(t_ephemeris),intent(in)               :: ephemeris                ! Ephemeris contains state vector and covariance
    !     !Output
    !     real(dp),dimension(n,2*n+1),intent(out)     :: x_chi                    ! Sigmapoint-matrix
    !     logical lsuccess

    !     !Locals
    !     integer                                     :: i                        ! Loop-counter
    !     integer                                     :: j                        ! Loop-counter
    !     real(dp)                                    :: beta                     ! Knowledge of distribution of state-vector (Gaussian distribution: beta=2)
    !     real(dp)                                    :: lambda                   ! Scaling parameter
    !     real(dp)                                    :: kappa                    ! Secondary scaling parameter
    !     real(dp),dimension(n,n)                     :: P_sqrt                   ! Root of Covariance-matrix
    !     real(dp),dimension(n,n)                     :: sqroot                   ! Root (dummy)

    !     character(len=*), parameter                 :: csubid = 'unscented_transformation'

    !     if(isControlled()) then
    !       if(hasToReturn()) return
    !       call checkIn(csubid)
    !     end if

    !     lsuccess = .false.
    !     sqroot = 0.0d0

    !     if(weighting == 1) then
    !         !!Nach Haykin.2001
    !         beta=2.0d0
    !         kappa=3.d0-n
    !         lambda=alpha**2*(n+kappa)-n
    !     else if(weighting == 2) then
    !         !!Nach wan.2000a
    !         beta=2.0d0
    !         kappa=0.d0
    !         lambda=alpha**2*(n+kappa)-n
    !     else if(weighting == 3) then
    !         !!Nach Julier.1997
    !         kappa=3-n
    !     else
    !         call setError(E_UNKNOWN_PARAMETER, FATAL, (/toString(weighting)/))
    !         return
    !     end if

    !     ! Cholesky decomposition
    !     !call invertMatrix(ephemeris%covariance_matrix%elem,P_sqrt,"CHOLESKY")
    !     if (.not. cholesky(n, ephemeris%covariance_matrix%elem, P_sqrt)) then
    !         return
    !     end if

    !     !Determine sigma point based on 7.30 (x_chi(1)=x)
    !     if(weighting==1 .or. weighting==2) then
    !         sqroot=sqrt(n+lambda)*P_sqrt    !Based on Haykin
    !     else if(weighting==3) then
    !         sqroot=sqrt(n+kappa)*P_sqrt     !Based on Julier 1997
    !     else
    !         call setError(E_UNKNOWN_PARAMETER, FATAL, (/toString(weighting)/))
    !         return
    !     end if

    !     x_chi(:,1)=ephemeris%state_vector(:)

    !     do j=2, n+1
    !         do i=1, n
    !             x_chi(i,j)=ephemeris%state_vector(i)+sqroot(i,j-1)
    !             x_chi(i,j+n)=ephemeris%state_vector(i)-sqroot(i,j-1)
    !         end do
    !     end do

    !     lsuccess = .true.

    !     if(isControlled()) then
    !         call checkOut(csubid)
    !     end if

    ! end function unscented_transformation

    !!--------------------------------------------------------------------------
    !!
    !> @anchor    cholesky
    !!
    !> @brief     Subroutine for Cholesky Decomposition
    !!
    !> @author    Eduard Gamper
    !!
    !> @date      <ul>
    !!              <li> 03.01.2014 (initial design)</li>
    !!            </ul>
    !!
    !> @details  Cholesky-Decomposition nxn-matrix. Matrix HAS TO BE symmetric
    !!              and positive definite. Equations from Tapley B., Schutz B.,
    !!              Born G. "Statistical orbit determination" (AP, 2004)
    !!              (ISBN 0126836302), p.287
    !!
    !> @param[in] n - Number of elements in state vector
    !!
    !> @param[in] A - quadratic, symmetric nxn matrix A
    !!
    !> @param[out] L - lower triangular matrix
    !!
    !> @result lsuccess - .true. when the process was successful
    !!
    !!--------------------------------------------------------------------------
    function cholesky(n, A, L) result(lsuccess)
        use IEEE_ARITHMETIC
        integer, intent(in)                     ::  n                           ! Number of elements in state vector
        real(dp), dimension(n,n), intent(in)    ::  A                           ! Random symmetric matrix
        real(dp), dimension(n,n), intent(out)   ::  L                           ! Root of matrix A (lower triangular matrix)

        !Locals
        integer                                 :: i,j,k                        ! Loop-counter
        real(dp)                                :: summe1 = 0                   ! Random sum for caching
        real(dp)                                :: summe2 = 0                   ! Random sum for caching
        logical                                 :: p                            ! Flag symmetry
        logical                                 :: s                            ! Flag positive definite quadratic form
        logical                                 :: lsuccess
        character(len=255)                      :: cmess                        ! Error message

        character(len=*), parameter             :: csubid = 'cholesky'

        if(isControlled()) then
          if(hasToReturn()) return
          call checkIn(csubid)
        end if

        lsuccess = .false.
        s = symmetrie_func(n,A)
        if(s .eqv. .true.) then
            p=positiv_func(n,A)             !Wenn Berechnung der Determinante funktioniert, wieder aktivieren!
            if (p .eqv. .true.) then
                do i=1, n
                    summe1 = 0
                    do k=1,i-1
                        summe1 = summe1+L(i,k)**2
                    end do
                    L(i,i)=sqrt(A(i,i)-summe1)
                    do j=i+1, n
                        summe2 = 0
                        do k=1, i-1
                            summe2 = summe2+L(i,k)*L(j,k)
                        end do
                        !if(abs(A(j,i)-summe2)<=0.0000000001) then   !Sicherheit gegen NaN bei Abschnittsfehlern
                        !    L(j,i)=0
                        !else
                            L(j,i)=(A(j,i)-summe2)/L(i,i)
                        !end if
                        if (ieee_is_nan( L(j,i))) then
                             L(j,i) = 0.0d0
                        end if
                    end do
                    do j=i+1,n !! obere Dreiecksmatrix zu 0 setzen
                        L(i,j) = 0.0d0
                    end do
                end do
            else
                !write (*,*) "Cholesky decomposition not successful. Returning"
                cmess = "Cholesky decomposition not successful"
                call setError(E_SPECIAL, FATAL, (/cmess/))
                return
            end if
        else
            L = 0.0d0
        end if
        lsuccess = .true.

        if(isControlled()) then
            call checkOut(csubid)
        end if

    end function cholesky

    !!--------------------------------------------------------------------------
    !!
    !> @anchor    symmetrie_func
    !!
    !> @brief     Function for testing the symmetry of quadratic matrix
    !!
    !> @author    Eduard Gamper
    !!
    !> @date      <ul>
    !!              <li> 01.01.2014 (initial design)</li>
    !!            </ul>
    !!
    !> @param[in] n - Number of elements in state vector
    !!
    !> @param[in] A - quadratic, symmetric nxn matrix A
    !!
    !> result .true. when the matrix is symmetric
    !!
    !!--------------------------------------------------------------------------
    logical function symmetrie_func(n, A)

        !Input
        integer, intent(in)                 ::  n                               ! Number of elements in state vector
        real(dp), dimension(n,n), intent(in)::  A                               ! Random matrix

        !Locals
        integer                             :: i,j                              ! Loop-counter
        logical                             :: symm = .true.                    ! Cache
        character(len=255)                  :: cmess                            ! Error message

        character(len=*), parameter         :: csubid = 'symmetrie_func'

        if(isControlled()) then
          if(hasToReturn()) return
          call checkIn(csubid)
        end if

        outer:do i=1,n
            do j=i+1, n
                if(abs(abs(A(i,j))-abs(A(j,i))) >= 1e-7) then
        !           if(abs(abs(A(i,j))-abs(A(j,i))) >= 0.00000001) then
                    symm=.false.
                    !write(*,*) 'Covariance matrix not symetric'
                    !write(*,*) 'i: ', i, 'j: ', j
                    !write (*,*) 'delta: ', abs(A(i,j))-abs(A(j,i))
                    !write(*,*) 'symmetrie A(2,1), A(1,2): ', A(j,i), '|', A(i,j)
                    cmess = "Covariance matrix not symetric"//new_line('A') &
                    //"i: "//toString(i)//"j: "//toString(j)//new_line('A') &
                    //"delta: "//toString(abs(A(i,j))-abs(A(j,i)))//new_line('A') &
                    //"symmetrie A(2,1), A(1,2): "//toString(A(j,i))//" | "//toString(A(i,j))
                    call setError(E_SPECIAL, FATAL, (/cmess/))
                    return
                    exit outer                ! Exit program on first .false.
                end if
            end do
            if (i == n) then
                symm=.true.
            end if
        end do outer

        symmetrie_func = symm

        if(isControlled()) then
            call checkOut(csubid)
        end if

        return
    end function symmetrie_func


    !!--------------------------------------------------------------------------
    !!
    !> @anchor    Positive definite quadratic form
    !!
    !> @brief     Function for testing nxn matrix for positive
    !!              definite quadratic form.
    !!
    !> @author    Eduard Gamper
    !!
    !> @date      <ul>
    !!              <li> 01.01.2014 (initial design)</li>
    !!            </ul>
    !!
    !> @param[in] n - Number of elements in state vector
    !!
    !> @param[in] A - quadratic, symmetric nxn matrix A
    !!
    !> result .true. when the matrix is symmetric
    !!
    !!--------------------------------------------------------------------------
    logical function positiv_func(n,A)

        !Input
        integer, intent(in)                 ::            n                     ! Number of elements in state vector
        real(dp), dimension(n,n), intent(in)::            A                     ! Random matrix

        !Locals
        integer     ::  i,j,k,l                                                 ! Loop-counter
        real(dp)    ::  det                                                     ! Determinant
        real(dp), allocatable , dimension(:,:)::matrix                          ! Matrix for caching

        character(len=*), parameter     :: csubid = 'positiv_func'

        if(isControlled()) then
          if(hasToReturn()) return
          call checkIn(csubid)
        end if

        if(n==2) then
          if(A(1,1)*A(2,2)-A(1,2)*A(2,1)>=0) then
            positiv_func=.true.
          else
            positiv_func=.false.
          end if
          return
        end if

        do i=3,n                     !Bestimmen der Determinante f¸r Matrizen ab 3x3 bis nxn
          allocate(matrix(i,i))
          do j=1, i
            do k=1, i
              matrix(k,j)=A(k,j)
            end do
          end do
          det=determinant(i, matrix)
          if (det<=0) then
            call setError(E_MATRIX_NOT_POSITIVE_DEFINITE, FATAL)
            !write(*,*) 'Covariance in ', i, '. minor not positiv defined.'
            positiv_func=.false.
            return
          end if
          deallocate(matrix)
        end do

        positiv_func=.true.

        if(isControlled()) then
            call checkOut(csubid)
        end if

        return
        201 format(6f18.10)
    end function positiv_func

    !--------------------------------------------------------------------------
    !
    !> @anchor    determinant
    !!
    !> @brief     Function for determing the determinant of nxn-matrix
    !!              for dimensions greater then 3x3
    !!
    !> @author    Eduard Gamper
    !!
    !> @details   Using a recursive calculation with Laplace's formula.
    !!
    !> @date      <ul>
    !!              <li> 01.01.2014 (initial design)</li>
    !!            </ul>
    !!
    !> @param[in] n - Number of elements in state vector
    !!
    !> @param[in] A - quadratic, symmetric nxn matrix A
    !!
    !> result .true. when the matrix is symmetric
    !!
    !!--------------------------------------------------------------------------
    recursive function determinant(n, A) result(summe)

        !Input
        integer, intent(in)                 ::  n                               ! Number of elements in state vector
        real(dp), dimension(n,n), intent(in)::  A                               ! Random matrix

        !Locals
        integer     ::  posi                                                    ! Loop-counter for column that is developed
        integer     ::  i,k,l                                                   ! Loop-counter for row
        integer     ::  j                                                       ! Loop-counter for column
        real(dp)    ::  summe                                                   ! Random sum for caching
        real(dp)    ::  zw_summe                                                ! Random sum for caching
        real(dp), dimension(n-1,n-1)::    S                                     ! Sub matrix

        character(len=*), parameter     :: csubid = 'determinant'

        if(isControlled()) then
          if(hasToReturn()) return
          call checkIn(csubid)
        end if

        if(n==2) then
          summe=A(1,1)*A(2,2)-A(1,2)*A(2,1)
        else
          do posi=1,n
            do j=1,n
            inner:     do i=2,n
              if (j .ne. posi) then
                if(j<posi) then
                  S(i-1,j)=A(i,j)
                else
                  S(i-1,j-1)=A(i,j)
                end if
              end if
            end do inner
            end do
          if(posi .eq. 1) then
          zw_summe=0
          end if
          zw_summe=zw_summe+A(1,posi)*((-1)**(posi+1))*determinant(n-1, S)
          end do
          summe=zw_summe
        end if

        if(isControlled()) then
            call checkOut(csubid)
        end if

        return
        201 format(6f18.10)
    end function determinant

end module rss_math
