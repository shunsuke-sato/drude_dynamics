! Drude model
! m d/dt v = q*E(t) - m*gamma*v

module global_variables
  implicit none
! math constants
  complex(8),parameter :: zi = (0d0, 1d0)
  real(8),parameter :: pi = 4d0*atan(1d0)

! physics constants
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: clight = 137.035999084d0
  real(8),parameter :: cm = 1d0/0.529177210903e-8
  real(8),parameter :: nm = 1d-9/0.529177210903d-10
  real(8),parameter :: fs = 1d0/0.024189d0

! Drude parameters
  real(8) :: ne_dns, mass, gamma, delta_gamma

! Physical system
  real(8) :: velocity_t

! time propagation
  integer :: nt
  real(8) :: Tprop, dt
  real(8), allocatable :: jt(:)

! frequency domain
  integer,parameter :: nw = 128
  real(8),parameter :: wi = 30d0*ev, wf = 50d0*ev, dw = (wf-wi)/nw
  complex(8) :: zsigma_w(0:nw), zsigma_gs_w(0:nw), zEw_store(0:nw)

! laser
  real(8),allocatable :: Eprobe(:)
  real(8),allocatable :: Ipump(:), Ipump_int(:)
  real(8) :: Tpump, Tprobe_train, omega_probe, omega_pump, tdelay, t_offset
  
  
end module global_variables
!------------------------------------------------------------------------------!
program main
  use global_variables
  implicit none
  

  call set_parameter
  call calc_stationary
!  call set_laser
  call pump_probe_dynamics

end program main
!------------------------------------------------------------------------------!
subroutine set_parameter
  use global_variables
  implicit none


  mass = 1d0
!  ne_dns = (15d0*ev)**2*mass/(4d0*pi)
!  ne_dns = (7.5d0*ev)**2*mass/(4d0*pi)
  ne_dns = (15.3d0*ev)**2*mass/(4d0*pi)


!  gamma = 7.5d0*ev
!  gamma = 1d0/(30d0/0.024189d3)

!  gamma = 0.5984d0*ev
  gamma = 1d0/(3d0/0.024189d0)
  delta_gamma = gamma*0.01d0


  Tprop = 200d0*fs
  dt = 0.04d0
  nt = aint(Tprop/dt)+1
  dt = Tprop/nt

end subroutine set_parameter
!------------------------------------------------------------------------------!
subroutine calc_stationary
  use global_variables
  implicit none
  real(8),parameter :: L_thick = 200d0*nm
  integer :: iw
  real(8) :: ww
  complex(8) :: zsigma, zeps, zn, zk, zt
  real(8) :: abs_coeff, k0, transmission
  real(8) :: n0_ox
  
!  n0_ox = 1.4d0
  n0_ox = 1d0

  open(20,file='absorption_Drude.out')
  do iw = 0, nw
    ww = wi + dw*iw

    zsigma = zi*ne_dns/mass/(ww+zi*gamma)
    zeps = 1d0 + 4d0*pi*zi*zsigma/ww
    zn = sqrt(zeps)

    k0 = ww/clight
    zk = k0*sqrt(zeps)

!exact
    zt = 2d0*zk*k0/&
      (2d0*zk*k0*cos(zk*L_thick)-zi*(zk**2+k0**2)*sin(zk*L_thick))
    transmission = abs(zt)**2

! approx
!    transmission = (1d0-abs(n0_ox-sqrt(zeps))**2/abs(n0_ox+sqrt(zeps))**2)**2&
!                  *(1d0-abs(1d0-n0_ox)**2/abs(1d0+n0_ox)**2)**2&
!                  *exp(-2d0*k0*aimag(sqrt(zeps))*L_thick)

    abs_coeff = 2d0*aimag(zn)*ww/clight*cm


    write(20,"(999e26.16e3)")ww/ev,transmission,abs_coeff,zsigma,zeps,-aimag(1d0/zeps)
  end do
  close(20)
  


end subroutine calc_stationary
!------------------------------------------------------------------------------!
subroutine set_laser
  use global_variables
  implicit none
  integer :: it
  real(8) :: tt, xx, ss

  allocate(Eprobe(-1:nt+1), Ipump(-1:nt+1), Ipump_int(-1:nt+1))
  omega_pump = 1.5d0*ev
  omega_probe = omega_pump*27d0
  Tpump = 20d0*fs
  Tprobe_train = 10d0*fs
!  tdelay = 0d0*fs

  t_offset = min(-0.5d0*Tpump, -0.5d0*tprobe_train+tdelay)

  Ipump = 0d0
  Ipump_int = 0d0
  Eprobe = 0d0
! pump
  do it = 0, nt+1
    tt = dt*it
    xx = tt + t_offset
    if(abs(xx) <= 0.5d0*Tpump)then
      Ipump(it) = cos(pi*xx/Tpump)**2
    end if
  end do

  ss = 0d0
  do it = 0, nt+1
    ss = ss + 0.5d0*(Ipump(it)+Ipump(it-1))*dt
    Ipump_int(it) = ss
  end do
  ss = Ipump_int(nt+1)
  Ipump_int = Ipump_int/ss
  Ipump = Ipump/ss
  

! probe
  do it = 0, nt+1
    tt = dt*it
    xx = tt + t_offset-tdelay
    if(abs(xx) <= 0.5d0*Tprobe_train)then
      Eprobe(it) = cos(pi*xx/Tprobe_train)**2*cos(omega_pump*xx)**6*sin(omega_probe*xx)
    end if
  end do




end subroutine set_laser
!------------------------------------------------------------------------------!
subroutine unset_laser
  use global_variables
  implicit none

  deallocate(Eprobe, Ipump, Ipump_int)
end subroutine unset_laser
!------------------------------------------------------------------------------!
subroutine pump_probe_dynamics
  use global_variables
  implicit none
  integer :: ndelay, idelay
  real(8) :: tdelay_ini, tdelay_fin
  integer :: iw, it
  real(8) :: jule_heat,jule_heat0

  allocate(jt(0:nt))

  tdelay_ini = -20d0*fs
  tdelay_fin =  20d0*fs
  ndelay = 100

! no-pump calculation
  tdelay = 0d0
  call set_laser
  open(20,file='laser.dat')
  do it = 0, nt
    write(20,"(999e26.16e3)")dt*it,Ipump(it),Ipump_int(it),Eprobe(it)
  end do
  close(20)
  Ipump = 0d0; Ipump_int = 0d0
  call calc_dynamics
  jule_heat0 = sum(Eprobe(0:nt)*jt(0:nt))*dt
  call calc_conductivity
  open(30,file='sigma_gs.out')
  do iw = 0, nw
    write(30,"(999e26.16e3)")(wi+dw*iw)/ev,zsigma_w(iw),abs(zEw_store(iw))**2
  end do
  close(30)
  call unset_laser
  zsigma_gs_w = zsigma_w

  open(50,file='tr_jule_heat.out')
  open(40,file='tr_sigma.out')
  do idelay = 0, ndelay
    tdelay = tdelay_ini + idelay*(tdelay_fin-tdelay_ini) /ndelay
    call set_laser
    call calc_dynamics
    jule_heat = sum(Eprobe(0:nt)*jt(0:nt))*dt
    call calc_conductivity
    do iw = 0, nw
      write(40,"(999e26.16e3)")tdelay/fs,(wi+dw*iw)/ev,zsigma_w(iw)-zsigma_gs_w(iw)
    end do
    write(40,*)

    write(50,"(999e26.16e3)")tdelay/fs,jule_heat,jule_heat-jule_heat0

    call unset_laser
  end do
  close(40)
  close(50)

end subroutine pump_probe_dynamics
!------------------------------------------------------------------------------!
subroutine calc_dynamics
  use global_variables
  implicit none
  integer :: it

  jt = 0d0
  velocity_t = 0d0
  jt(0) = velocity_t*ne_dns/mass

  do it = 0, nt
    call dt_evolve(it)
    jt(it+1) = velocity_t*ne_dns/mass
  end do


end subroutine calc_dynamics
!------------------------------------------------------------------------------!
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer, intent(in) :: it
  real(8) :: gamma_eff, gamma_n, gamma_o

  gamma_o = gamma + delta_gamma*Ipump_int(it)
  gamma_n = gamma + delta_gamma*Ipump_int(it+1)
  
  gamma_eff = 0.5d0*(gamma_n + gamma_o)
  
  velocity_t = velocity_t*exp(-dt*gamma_eff) &
    +0.5d0*dt*(Eprobe(it+1) + Eprobe(it) )/mass

end subroutine dt_evolve
!------------------------------------------------------------------------------!
subroutine calc_conductivity
  use global_variables
  implicit none
  integer :: iw, it
  real(8) :: ww,tt
  complex(8) :: zjw, zEw, zfact

  zsigma_w = 0d0

  do iw = 0, nw
    ww = wi + dw*iw

    zjw = 0d0
    zEw = 0d0
    do it = -1,nt+1
      tt = dt*it
      zfact = exp(zI*ww*tt)
      zjw = zjw + jt(it)    *zfact
      zEw = zEw + Eprobe(it)*zfact
    end do
    zsigma_w(iw) = zjw/zEw
    zEw_store(iw) = zEw*dt

  end do



end subroutine calc_conductivity
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
