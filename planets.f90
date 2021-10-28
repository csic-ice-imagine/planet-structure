program main
  implicit none

  integer, parameter :: nr = 1000, np = 50           ! radial resolution
  real*8, parameter :: plogmin = 5., plogmax = 10.    ! Min and Max log(P[bar]) to be explored
  real*8 :: pc                                       ! Pressure in the center in bar
  character, dimension(4), parameter :: comp = ["f","r","w","g"]
  integer :: i, c

  do c = 4, 4

    open(unit = 14, file = comp(c) // '_mass-radius.d', status='replace')
    write(14, *) "Pc [Mbar], R [10^7 m], Mearth: "
    close(14)
  
    do i = 0, np
      pc = 10.**(plogmin + i*(plogmax - plogmin)/np)
      call structure(pc, comp(c), nr)
    end do

  end do


  return

end program



subroutine structure(pc, comp, nr)
  implicit none
  integer, intent(in) :: nr           ! radial resolution
  real*8, intent(in) :: pc            ! Pressure in the center in bar
  character, intent(in) :: comp
  real*8, parameter :: pi = acos(-1d0)
  real*8, parameter :: UNIT_R = 1d7        ! in m
  real*8, parameter :: UNIT_RHO = 1d3      ! Density g/cm^3 in terms of kg/m^3
  real*8, parameter :: UNIT_M = UNIT_RHO*UNIT_R**3   ! Unit mass in kg
  real*8, parameter :: EARTH_M = 5.97e24   ! in kg
  real*8, parameter :: UNIT_P = 1d5        ! Pressure bar in terms of Pascal
  real*8, parameter :: G = 6.673d-11/(UNIT_P*UNIT_R/(UNIT_M*UNIT_RHO))    ! gravitational constant in terms of the other quantities
  real*8, parameter :: ps = 1d0            ! Pressure at the surface in bar
  real*8, parameter :: eps = 1.0         ! First point: dr(1->2)/r(1) it has to be of order 1
  real*8, parameter :: gamma = 2.5         ! Adiabatic index in P(T)
  real*8, parameter :: logk_ad = 8       ! Constant of P = k_ad T^gamma
  real*8, dimension(nr) :: r, rho, p, temp, mass       ! P is in Mbar
  real*8 :: dlnp, mdm, rho_eos
  real*8 :: pint, rhoint, rint, mint
  real*8 :: dm_1, dm_2, dm_3, dm_4, dr_1, dr_2, dr_3, dr_4
  integer i

  dlnp = dlog(ps/pc)/(nr-1.)
  do i = 1, nr
   p(i) = pc*dexp((i-1)*dlnp)
   temp(i) = (p(i)/10**logk_ad)**(1/gamma)
   rho(i) = rho_eos(p(i),temp(i),comp)
  enddo
   ! dr(1->2) = -dlnp p(1) r(1)**2 / (m(1) G rho(1)) = -dlnp p(1) 3/(4pi rho(1)^2 G r(1))
  ! Fxi r(1->2) so that dr(1)/r(1) simeq eps < 1
  r(1) = eps*dsqrt(-3d0*dlnp*pc/(4d0*pi*rho(1)**2*G))
  mass(1) = (4d0/3d0)*pi*rho(1)*(r(1))**3

! EULER   
  do i = 2, nr
    mdm = - dlnp*4*pi*(r(i-1))**4*p(i)/G
    mass(i) = dsqrt( mass(i-1)**2 + 2d0*mdm)
    r(i) = ( r(i-1)**3 + 3d0*(mdm/mass(i))/(4d0*pi*rho(i)) )**(1d0/3d0)
  enddo

  open(unit = 12, file = 'eul.d', status='replace')
  write(12, *) "r, rho, p, mass"
  do i = 1, nr
    write(12, '(1pe11.3,1pe11.3,1pe11.3,1pe11.3)') r(i)*UNIT_R, rho(i), p(i), mass(i)*UNIT_M/EARTH_M
  end do    
  close(12)

  print*, "Euler:   Pc, R, M: ", pc/1e6, r(nr), mass(nr)*UNIT_M/EARTH_M

! RK4 WITH VOLUME INTEGRATION
  do i = 2, nr
   ! Runge-Kutta 4th order to obtain m and r
   ! 1st sub-step: previous point
   dm_1 = - dlnp*4*pi*(r(i-1))**4*p(i-1)/(mass(i-1)*G)
   dr_1 = dm_1/(4d0*pi*rho(i-1)) / r(i-1)**2

   ! 2nd substep: updated intermediate values of p and rho
   mint = sqrt(mass(i-1)**2 + mass(i-1)*dm_1)
   rint = ( r(i-1)**3 + 0.5d0*3d0*dr_1*r(i-1)**2 )**(1d0/3d0)
   pint = pc*dexp((i-0.5d0)*dlnp)
   rhoint = rho_eos(pint,0.5d0*(temp(i)+temp(i-1)),comp)
   dm_2 = -dlnp*4*pi*rint**4*pint/(mint*G)
   dr_2 = dm_2/(4d0*pi*rhoint) / rint**2

   ! 3rd substep: same values of p and rho, update m and r
   mint = sqrt(mass(i-1)**2 + mint*dm_2)
   rint = ( r(i-1)**3 + 0.5d0*3d0*dr_2*rint**2 )**(1d0/3d0)
   dm_3 = -dlnp*4*pi*rint**4*pint/(mint*G)
   dr_3 = dm_3/(4d0*pi*rhoint) / rint**2

   ! 4th substep: new values of p and rho, update m and r
   mint = sqrt(mass(i-1)**2 + 2d0*mint*dm_1)
   rint = ( r(i-1)**3 + 3d0*dr_1*rint**2 )**(1d0/3d0)
   dm_4 = -dlnp*4*pi*rint**4*p(i)/(mint*G)
   dr_4 = dm_4/(4d0*pi*rho(i)) / rint**2

   mass(i) = mass(i-1) + (dm_1 + 2d0*dm_2 + 2d0*dm_3 + dm_4)/6d0
   r(i) = r(i-1) + (dr_1 + 2d0*dr_2 + 2d0*dr_3 + dr_4)/6d0

 enddo

 open(unit = 12, file = 'rk4vol.d', status='replace')
 write(12, *) "r, rho, p, mass"
 do i = 1, nr
   write(12, '(1pe11.3,1pe11.3,1pe11.3,1pe11.3)') r(i)*UNIT_R, rho(i), p(i), mass(i)*UNIT_M/EARTH_M
 end do    
 close(12)

 print*, comp, "RK4 vol: Pc, R, M: ", pc/1e6, r(nr), mass(nr)*UNIT_M/EARTH_M
 open(unit = 14, file = comp // '_mass-radius.d', access='append')
 write(14,*) pc/1e6, r(nr), mass(nr)*UNIT_M/EARTH_M
 close(14)

! RK4 WITH DIFFERENCES
 do i = 2, nr
   ! Runge-Kutta 4th order to obtain m and r
   ! 1st sub-step: previous point
   dm_1 = -dlnp*4*pi*(r(i-1))**4*p(i-1)/(mass(i-1)*G)
   dr_1 = dm_1/(4d0*pi*rho(i-1)) / r(i-1)**2

   ! 2nd substep: updated intermediate values of p and rho
   mint = mass(i-1) + 0.5d0*dm_1
   rint = r(i-1) + 0.5d0*dr_1
   pint = pc*dexp((i-0.5d0)*dlnp)
   rhoint = rho_eos(pint,0.5d0*(temp(i)+temp(i-1)),comp)
   dm_2 = -dlnp*4*pi*rint**4*pint/(mint*G)
   dr_2 = dm_2/(4d0*pi*rhoint) / rint**2

   ! 3rd substep: same values of p and rho, update m and r
   mint = mass(i-1) + 0.5d0*dm_2
   rint = r(i-1) + 0.5d0*dr_2
   dm_3 = -dlnp*4*pi*rint**4*pint/(mint*G)
   dr_3 = dm_3/(4d0*pi*rhoint) / rint**2

   ! 4th substep: new values of p and rho, update m and r
   mint = mass(i-1) + dm_3
   rint = r(i-1) + dr_3
   dm_4 = -dlnp*4*pi*rint**4*p(i)/(mint*G)
   dr_4 = dm_4/(4d0*pi*rho(i)) / rint**2

   mass(i) = mass(i-1) + (dm_1 + 2d0*dm_2 + 2d0*dm_3 + dm_4)/6d0
   r(i) = r(i-1) + (dr_1 + 2d0*dr_2 + 2d0*dr_3 + dr_4)/6d0

 enddo

 open(unit = 12, file = 'rk4.d', status='replace')
 write(12, *) "r, rho, p, mass"
 do i = 1, nr
   write(12, '(1pe11.3,1pe11.3,1pe11.3,1pe11.3)') r(i)*UNIT_R, rho(i), p(i), mass(i)*UNIT_M/EARTH_M
 end do    
 close(12)
 print*, "RK4 fd:  Pc, R, M: ", pc/1e6, r(nr), mass(nr)*UNIT_M/EARTH_M


end subroutine structure



! p[Mbar] = 3.59d-5*rho*temp          ! eq. (2) Fortney et al. 2007 & Fig. 1 H20 very rough
real*8 function rho_eos(pres_bar, temperature, comp)
  implicit None
  real*8, intent(in) :: pres_bar, temperature
  character, intent(in) :: comp
  real*8, parameter :: B = 3.59e-5 ! Term hot ideal EoS, eq. (2) Fortney

  select case(comp)

    case("f")
      rho_eos = 11d0*(pres_bar*1e-6)**0.28      ! Iron EoS from Fortney, manual fit
    case("r")
      rho_eos = 5d0*(pres_bar*1e-6)**0.32      ! Rock EoS from Fortney, manual fit
    case("w")
      rho_eos = 2.45d0*(pres_bar*1e-6)**0.38    ! H20 EoS from Fortney, manual fit
    case("g")
      rho_eos = pres_bar*1e-6/(B*temperature)
  end select

end function rho_eos
