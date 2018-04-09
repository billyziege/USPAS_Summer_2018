program MDelectron
!unit setting: using electron@LAMMPS
!atomic mass units - Bohr radius - fs - Hartree
!!!velocity = [Bohr radius/atomic time unit] ~ 1.03275 fs
!Force = [Hartrees/Bohr radius]


  implicit none
  !system parameters
  real(8), parameter :: Pi=4*atan(1._8)
  real(8), parameter :: tt = 1.03275  !time-factor due to units
  real(8), parameter :: vc =  5.85E3  !speed of light
  real(8), parameter :: m = 5.4858E-4
  
  real(8), parameter :: cutoff2 = 400.0 !proximity limit
  integer, parameter :: N=8000, Ntime=100
  real(8) :: dt=1.0, realt = 0.0
  integer :: plotstride = 10
  !experiment: lz = 0.4 um, rx = ry = 100.0 um with 10^6 electrons
  real(8), parameter :: L = 3.779E5  , Lz = 1511.8
  
  !using rest mass of electrondd
  !define the Electronic Field force coefficient in [Hartree/(a_0*electron_charge)]
  !to make F=eE agree with the Force unit
  real(8), parameter :: E_coeff=1.94E-12
  !Extraction Field in [Volts/meter]
  real(8), parameter :: EE = 5.0E6  ! for 5MV/m
  ! z_anode = 2cm to make 100keV electrons
  real(8), parameter :: z_anode_on = 1.89E8, z_anode_off = 3*z_anode_on
  real(8), parameter :: z_anode = 2*z_anode_on
  !Relativistic effect (c in velocity unit)
  
  real(8) :: R(3,N), V(3,N), F(3,N), rel(3)
  real(8) :: r1, r2, r3,r4, PE, KE, rellength, s1,s2,s3
  !R is for position of the atoms, V for velocity, F for Force
  !real(8) :: KEnergy(Ntime),PE,PE_0, Etot(Ntime), Eofel, Eofel_0
  !KEnergy = total Kinetic energy of atomic system
  !PE and PE_0 is the total Potential energy of atoms
  !Eofel is energy of electronic system
  integer :: Time,i,j,k,numb,check
   
  open(UNIT=13, file="RandV_3D_Uniform.xyz", status="replace")
   
  !-----Initialization-----
  call init_R_Uniform(R,N)
  V = 0.0
  F = 0.0
  
  !Warm-up for the simulation
  do Time=1,Ntime
    call verlet_init(V,R,dt,N,F,m)
    realt = realt + dt
    if (modulo(Time,plotstride) == 1) then
      call getPE(R,N,PE)
      KE = 0.5*m*sum(V**2)
      write(*,'(4F15.5)') realt, PE, KE, PE+KE
      write(13,'(i5)') N
      write(13,'(F15.5)') realt
      do i =1, N
        write(13,'(i5,3F15.3,3F15.7)') 1,R(:,i),V(:,i)
      end do
    end if  
    !change time step size for better resolution
    if (Time == 100) then
      dt = 1.0
      print *, 'dt = ', dt
    endif
    if (Time == 200) then
      dt = 5.0
      print *, 'dt = ', dt
    endif
    if (Time == 500) then
      dt = 15.0
      print *, 'dt = ', dt
    endif
    if (Time == 700) then
      dt = 50.0
      print *, 'dt = ', dt
    endif
    if (Time == 2000) then
      plotstride = 100
    endif
  end do

  close(13)

contains
 
  subroutine Init_R_Uniform(R,N)
    real(8), intent(inout) :: R(:,:)
    integer, intent(in) :: N
    integer :: numb,check,i
    real(8) :: r1,r2,r3,s1,s2,s3,rel(3)
    call random_seed()
    numb = 0
    do while (numb < N)
      call random_number(r1)
      call random_number(r2)
      call random_number(r3)
      !Uniform distribution
      s1 = 2.0*r1-1.0
      s2 = 2.0*r2-1.0
      s3 = 2.0*r3-1.0
      !check for proximity
      check = 0
      if (sqrt(s1*s1+s2*s2+s3*s3) < 1.0) then
        s1 = L*s1
        s2 = L*s2
        s3 = Lz*s3
        do i = 1, numb
          rel = R(:,i) - (/s1,s2,s3/)
          if (sum(rel**2) < cutoff2 ) then
            check = 1
            exit
          endif
        enddo
        if (check == 0) then
          numb = numb + 1
          R(:,numb) = (/s1,s2,s3/)
        endif
      endif
    enddo
  end subroutine Init_R_Uniform

  subroutine verlet_init(V,R, dt, N, F, m)
    real(8), intent(inout) :: V(:,:), R(:,:), F(:,:)
    real(8), intent(in) :: dt, m
    integer , intent(in) :: N
    real(8) :: dttt,hdtttm
    dttt = dt*tt
    hdtttm = 0.5*dt*tt/m
      V=V+F*hdtttm
      R=R+V*dttt
      call Force(F,R,N)
      V=V+F*hdtttm
  end subroutine verlet_init
 
  !Force Calculation
  subroutine Force(F,R,N)
    real(8), intent(inout) :: F(:,:)
    real(8), intent(in)   :: R(:,:)
    integer , intent(in) :: N
    integer :: i,j 
    real(8) :: REL(3), RELlength, Fij(3), dz_ee
    real(8) :: zforee !this create a register for R(3,i) to save time
    F=0.0
    !$omp parallel do private(fij,rel,rellength)
    do j=1,N-1
      do i=j+1,N
        REL=R(:,i)-R(:,j)
        !REL=REL-nint(Rel/L)*L
        RELlength=sqrt(sum(REL**2))
        !FC = 1.0 so we can skip it
        Fij=(1/(RELlength)**3)*REL
        !$omp atomic
        F(1,j)=F(1,j)-Fij(1)
        !$omp atomic
        F(2,j)=F(2,j)-Fij(2)
        !$omp atomic
        F(3,j)=F(3,j)-Fij(3)
        !$omp atomic
        F(1,i)=F(1,i)+Fij(1)
        !$omp atomic
        F(2,i)=F(2,i)+Fij(2)
        !$omp atomic
        F(3,i)=F(3,i)+Fij(3)
      end do 
    end do 
    !$omp end parallel do

    !the Extraction Field on z-dirction
!    do i = 1, N
!      zforee = R(3,i)
!      if ( zforee < z_anode_off) then 
!        if (zforee < z_anode_on) then
!          F(3,i) = F(3,i) + E_coeff*EE 
!        else
!          dz_ee = (zforee - z_anode)/z_anode_on
!          F(3,i) = F(3,i) + E_coeff*0.5*EE*(1.0-tanh(5.0*dz_ee))
!        endif
!      endif
!    enddo

  end subroutine Force

  !get the potential energy of Atomic system
  subroutine getPE(R,N,PE)
  real(8), intent(in) :: R(:,:)
  integer, intent(in) :: N
  real(8), intent(out):: PE
  integer :: i,j
  real(8) :: rel_pe(3),rel_2,vij
  PE = 0.0
  do j=1,N-1
    do i=j+1,N
      rel_pe=R(:,i)-R(:,j)
      rel_2=sqrt(sum(rel_pe**2))
      vij = 1.0/(rel_2)
      PE = PE + vij
    end do
  end do
  end subroutine getPE

end program MDelectron
