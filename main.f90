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
  real(8) :: ne_dns, mass, gamma

! Physical system
  real(8) :: velocity_t

! time propagation
  integer :: nt
  real(8) :: Tprop, dt

! laser
  real(8),allocatable :: Eprobe(:)
  real(8),allocatable :: Ipump(:)
  real(8) :: Tpump, Tprobe_train, omega_probe, omega_pump, tdelay, t_offset
  
  
end module global_variables
!------------------------------------------------------------------------------!
program main
  use global_variables
  implicit none
  

  call set_parameter
  call calc_stationary
  call set_laser

end program main
!------------------------------------------------------------------------------!
subroutine set_parameter
  use global_variables
  implicit none


  mass = 1d0
!  ne_dns = (15d0*ev)**2*mass/(4d0*pi)
  ne_dns = (7.5d0*ev)**2*mass/(4d0*pi)


!  gamma = 7.5d0*ev
  gamma = 1d0/(30d0/0.024189d3)


  Tprop = 100d0*fs
  dt = 0.1d0
  nt = aint(Tprop/dt)+1
  dt = Tprop/nt

end subroutine set_parameter
!------------------------------------------------------------------------------!
subroutine calc_stationary
  use global_variables
  implicit none
  integer,parameter :: nw = 512
  real(8),parameter :: wi = 0.1d0*ev, wf = 90d0*ev
  real(8),parameter :: dw = (wf-wi)/nw
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


    write(20,"(999e26.16e3)")ww/ev,transmission,abs_coeff,zeps,-aimag(1d0/zeps)
  end do
  close(20)
  



end subroutine calc_stationary
!------------------------------------------------------------------------------!
subroutine set_laser
  use global_variables
  implicit none
  integer :: it
  real(8) :: tt, xx, ss

  allocate(Eprobe(-1:nt+1), Ipump(-1:nt+1))
  omega_pump = 1.5d0*ev
  omega_probe = omega_pump*27d0
  Tpump = 40d0*fs
  Tprobe_train = 10d0*fs
  tdelay = 0d0*fs

  t_offset = min(-0.5d0*Tpump, -0.5d0*tprobe_train+tdelay)

  Ipump = 0d0
  Eprobe = 0d0
! pump
  do it = 0, nt
    tt = dt*it
    xx = tt + t_offset
    if(abs(xx) <= 0.5d0*Tpump)then
      Ipump(it) = cos(pi*xx/Tpump)**2
    end if
  end do
  ss = sum(Ipump)*dt
  Ipump = Ipump/ss

! probe
  do it = 0, nt
    tt = dt*it
    xx = tt + t_offset-tdelay
    if(abs(xx) <= 0.5d0*Tprobe_train)then
      Eprobe(it) = cos(pi*xx/Tprobe_train)**2*cos(omega_pump*xx)**6*sin(omega_probe*xx)
    end if
  end do


  open(20,file='laser.dat')
  do it = 0, nt
    write(20,"(999e26.16e3)")dt*it,Ipump(it),Eprobe(it)
  end do
  close(20)




end subroutine set_laser
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
