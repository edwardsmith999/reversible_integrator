module stuff

contains

subroutine integrate_normally(QMM, QM, Q0, QP, iters, dt)
    implicit none

    integer, intent(in) :: iters
    double precision, intent(inout) :: QMM, QM, Q0, QP

    integer :: it
    double precision ::  QPP, AP, A0, AM, dt, dt2

    dt2 = dt**2

    !Integrate a number of iterations
    DO IT = 1,iters

        ! Acceletaton is just momentum as this example is
        ! a harmonic oscillator \ddot{q} = −q
        AP = QP
        A0 = Q0
        AM = QM

        !Compute Ingredients Of Three Accelerations
        AP = 0.25d0*dt2*(5.D0*AP)
        A0 = 0.25d0*dt2*(2.D0*A0)
        AM = 0.25d0*dt2*(5.D0*AM)

        ! This is the "step" (DOI:10.12921/cmst.2017.0000031)
        QPP = QP + QM - QMM - (AP + A0 + AM)

        ! ======== Print stuff =========
        !Get momenta (4th order) (http://arxiv.org/abs/1801.09899v3)
!        IP0 =   (4.d0/3.d0)*(QP -QM )/(2.d0*dt) & 
!              - (1.d0/3.d0)*(QPP-QMM)/(4.d0*dt)
!        print*, "time, pos, momentum", it*dt, Q0, P0
        ! ======== Print stuff =========

        !Coordinate Updates For Five Successive Times
        QMM = QM
        QM  = Q0
        Q0  = QP
        QP  = QPP

    END DO

end subroutine integrate_normally


subroutine integrate_reversably(QMM, QM, Q0, QP, iters, dt)
    implicit none

    integer, intent(in) :: iters
    double precision, intent(inout) :: QMM, QM, Q0, QP

    integer :: it
    integer, parameter :: intkind = selected_int_kind(8)
    integer(kind=intkind) :: IQMM,IQM,IQ0,IQP,IQPP,IP0 !Contiguous integer coordinates
    integer(kind=intkind) :: IAP,IA0,IAM           !Ingredients of the accelerations
    double precision ::dt, dt2

    dt2 = dt**2

    !Convert Coordinates To 15-Digit Integers
    IQMM = QMM*(10.D0**(intkind-1))
    IQM  = QM *(10.D0**(intkind-1))
    IQ0  = Q0 *(10.D0**(intkind-1))
    IQP  = QP *(10.D0**(intkind-1))

    !Integrate a number of iterations
    DO IT = 1,iters

        ! Acceletaton is just momentum as this example is
        ! a harmonic oscillator \ddot{q} = −q
        IAP = IQP
        IA0 = IQ0
        IAM = IQM

        !Compute Ingredients Of Three Accelerations
        IAP = 0.25d0*dt2*(5.D0*IAP)
        IA0 = 0.25d0*dt2*(2.D0*IA0)
        IAM = 0.25d0*dt2*(5.D0*IAM)

        ! This is the "step" (DOI:10.12921/cmst.2017.0000031)
        IQPP = IQP + IQM - IQMM - (IAP + IA0 + IAM)

        ! ======== Print stuff =========
        !Get momenta (4th order) (http://arxiv.org/abs/1801.09899v3)
!        IP0 =   (4.d0/3.d0)*(IQP -IQM )/(2.d0*dt) & 
!              - (1.d0/3.d0)*(IQPP-IQMM)/(4.d0*dt)
!        print*, "time, pos, momentum", it*dt, IQ0/10.D0**15, IP0/10.D0**15
        ! ======== Print stuff =========

        !Coordinate Updates For Five Successive Times
        IQMM = IQM
        IQM  = IQ0
        IQ0  = IQP
        IQP  = IQPP

    END DO

    !Convert back to double precision
    QMM = IQMM/10.D0**(intkind-1)
    QM = IQM/10.D0**(intkind-1)
    Q0 = IQ0/10.D0**(intkind-1)
    QP = IQP/10.D0**(intkind-1)

end subroutine integrate_reversably


end module


program revsersable_integrator
    use stuff
    implicit none

    integer :: it, itmax

    ! QMM is q_{t-2dt}, QM is q_{t+dt}, 
    !     Q0 is q_{t}, P0 is p_{t}
    ! QP is q_{t-dt} and QPP is q_{t-2dt}
    double precision ::dt, dt2, QMM, QM, Q0, QP, P0, time, t1, t2
    double precision ::initP, init0, initM, initMM
    double precision ::tempP, temp0, tempM, tempMM

    !Timestep and number of steps
    ITMAX = 100000000
    dt = 0.005
    dt2 = dt**2

    !Setup initial condition
    QMM = DCOS(-2.D0*DT)
    QM = DCOS(-1.D0*DT)
    Q0 = DCOS( 0.D0*DT)
    QP = DCOS(+1.D0*DT)

    !==============================
    !    Test 1 - Integer
    !==============================

    !Store initial conditions
    initMM = QMM; initM = QM; init0 = Q0; initP = QP

    call cpu_time(t1)

    !Integrate
    call integrate_reversably(QMM, QM, Q0, QP, ITMAX, dt)

    !Integerate backwards
    call integrate_reversably(QP, Q0, QM, QMM, ITMAX, dt)

    call cpu_time(t2)

    !Check divergence from initial condition
    print'(a, 2f10.5, f24.17, a, f10.5 )', "Error Integer = ", Q0, init0, Q0-init0, " Cpu time = ", t2-t1

    !==============================
    !    Test 1 - Double precision
    !==============================

    !Store initial conditions
    initMM = QMM; initM = QM; init0 = Q0; initP = QP

    call cpu_time(t1)

    !Integrate
    call integrate_normally(QMM, QM, Q0, QP, ITMAX, dt)

    !Integerate backwards
    call integrate_normally(QP, Q0, QM, QMM, ITMAX, dt)

    call cpu_time(t2)

    !Check divergence from initial condition
    print'(a, 2f10.5, f24.17, a, f10.5 )', "Error double = ", Q0, init0, Q0-init0, " Cpu time = ", t2-t1

    !==============================
    !   Test 2 - stepwise
    !==============================

!    integer*16 :: IQMM,IQM,IQ0,IQP,IQPP,IP0 !Contiguous integer coordinates
!    integer*16 :: IAP,IA0,IAM           !Ingredients of the accelerations

!    !Store initial conditions
!    initMM = QMM; initM = QM; init0 = Q0; initP = QP

!    !Integrate
!    DO IT = 1,ITMAX
!        call integrate_reversably(QMM, QM, Q0, QP, 1, dt)
!    enddo
!    !Integerate backwards
!    DO IT = 1,ITMAX
!        call integrate_reversably(QP, Q0, QM, QMM, 1, dt)
!    enddo
!    !Check divergence from initial condition
!    print'(a, 2f10.5, f24.17)', "Error stepwise = ", Q0, init0, Q0-init0
!  
    !Convert Coordinates To 15-Digit Integers
!    IQMM = QMM*(10.D0**15)
!    IQM  = QM *(10.D0**15)
!    IQ0  = Q0 *(10.D0**15)
!    IQP  = QP *(10.D0**15)

!    DO IT = 1,ITMAX
!        TIME = IT*DT

!        ! Acceletaton is just momentum as this example is
!        ! a harmonic oscillator \ddot{q} = −q
!        IAP = IQP
!        IA0 = IQ0
!        IAM = IQM

!        !Compute Ingredients Of Three Accelerations
!        IAP = 0.25d0*dt2*(5.D0*IAP)
!        IA0 = 0.25d0*dt2*(2.D0*IA0)
!        IAM = 0.25d0*dt2*(5.D0*IAM)

!        ! This is the "step" (DOI:10.12921/cmst.2017.0000031)
!        IQPP = IQP + IQM - IQMM - (IAP + IA0 + IAM)

!        !Get momenta (4th order) (http://arxiv.org/abs/1801.09899v3)
!        IP0 =   (4.d0/3.d0)*(IQP -IQM )/(2.d0*dt) & 
!              - (1.d0/3.d0)*(IQPP-IQMM)/(4.d0*dt)

!        !Coordinate Updates For Five Successive Times
!        IQMM = IQM
!        IQM = IQ0
!        IQ0 = IQP
!        IQP = IQPP

!        !Convert back to double precision
!        Q0 = IQ0/10.D0**15
!        P0 = IP0/10.D0**15

!        print*, it, time, Q0, P0
!    END DO

!    !Swap order of initial conditions to reverse trajectory
!    tempMM = QMM; tempM = QM; temp0 = Q0; tempP = QP
!    QMM = tempP; QM = temp0; Q0 = tempM; QP = tempMM

end program revsersable_integrator
