module mod_usr
  use mod_hd
  use mod_eos, only: eos
  implicit none
  ! Massaglia et al. 2016, Eqs. (5)-(8), Case B, 3D Cartesian jet.
  ! Derived from local jet_cloud; original routine names/order retained.
  double precision :: eta_jet, mach, rc, alpha, rhoj, pj, pa, vj
  logical, parameter :: inject_jet = .true.

contains

  subroutine usr_init()
  
    call set_coordinate_system("Cartesian_3D")
    
    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_aux_output      => extra_var_output
    usr_add_aux_names   => extra_var_names_output

    call hd_activate()

  end subroutine usr_init
  
  subroutine initglobaldata_usr
    if (ndim /= 3 .or. ndir /= 3) call mpistop('Case B requires -d=3 -v=3')
    if (.not. hd_energy .or. hd_n_tracer /= 1) call mpistop('Require energy and one tracer')
    if (abs(eos%gamma-5.d0/3.d0)>1.d-12) call mpistop('Require gamma=5/3')
    if (abs(xprobmin2)>1.d-12) call mpistop('Require bottom boundary at y=0')
    if (nwauxio /= 1) call mpistop('Require nwauxio=1 for T output')
    if (nghostcells /= 2 .and. nghostcells /= 3) &
      call mpistop('This setup supports two or three ghost layers')
    if (mype == 0) print *, 'FRJETS3D: ndim, ndir, nghostcells = ', ndim, ndir, nghostcells
    ! Units: jet radius, jet density, jet sound speed. No gravity/cooling.
    eta_jet=0.01d0
    mach=4.d0
    rc=40.d0 ! Ambient CORE radius, no longer the old cloud radius.
    alpha=2.d0
    rhoj=1.d0
    pj=1.d0/eos%gamma
    pa=pj
    vj=mach*sqrt(eos%gamma*pj/rhoj)
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: radius(ixI^S)
    ! Initially only the static atmosphere: no cloud or pre-filled jet.
    radius(ixO^S)=sqrt(x(ixO^S,1)**2+x(ixO^S,2)**2+x(ixO^S,3)**2)
    w(ixO^S,1:nw)=0.d0
    w(ixO^S,rho_)=(rhoj/eta_jet)/(1.d0+(radius(ixO^S)/rc)**alpha)
    w(ixO^S,e_)=pa
    w(ixO^S,tracer(1))=0.d0
    ! Here mom slots contain velocity and e_ pressure; convert exactly once.
    call eos%to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initonegrid_usr

  subroutine specialbound_usr(qdt,qt, ixI^L,ixO^L, iB, w, x)
    integer, intent(in) :: ixI^L,ixO^L, iB
    double precision, intent(in) :: qdt,qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    integer :: i,j,k,jmirror
    double precision :: rad, sr, sq, rr, vxr, vzr, myr, pr, fr
    double precision :: den,vx,vz,my,press,frac
    select case(iB)
    case(3)
      ! Mirror EACH ghost layer. Normal momentum odd; all other U even.
      do j=ixOmin2,ixOmax2
        jmirror=2*ixOmax2+1-j
        w(ixOmin1:ixOmax1,j,ixOmin3:ixOmax3,1:nw)=w(ixOmin1:ixOmax1,jmirror,ixOmin3:ixOmax3,1:nw)
        w(ixOmin1:ixOmax1,j,ixOmin3:ixOmax3,mom(2))=-w(ixOmin1:ixOmax1,jmirror,ixOmin3:ixOmax3,mom(2))
      end do
      if (inject_jet) then
        do k=ixOmin3,ixOmax3
        do j=ixOmin2,ixOmax2
          do i=ixOmin1,ixOmax1
            rad=hypot(x(i,j,k,1),x(i,j,k,3)) ! Circular nozzle; hypot avoids squaring overflow.
            ! Cap BEFORE exponentiation: avoids both power and cosh overflow.
            sr=1.d0/cosh(min(rad/1.4d0,80.d0**(1.d0/8.d0))**8)
            sq=1.d0/cosh(min(rad,80.d0**(1.d0/6.d0))**6)
            rr=w(i,j,k,rho_)
            vxr=w(i,j,k,mom(1))/rr
            vzr=w(i,j,k,mom(3))/rr
            myr=w(i,j,k,mom(2))
            pr=(eos%gamma-1.d0)*(w(i,j,k,e_)- &
                 0.5d0*(w(i,j,k,mom(1))**2+myr**2+w(i,j,k,mom(3))**2)/rr)
            fr=w(i,j,k,tracer(1))/rr
            ! Eq.(8): mix density, transverse velocity, AXIAL MOMENTUM,
            ! pressure and tracer FRACTION, using two different weights.
            den=(1.d0-sr)*rr+sr*rhoj
            vx=(1.d0-sq)*vxr
            vz=(1.d0-sq)*vzr
            my=(1.d0-sq)*myr+sq*rhoj*vj
            press=(1.d0-sq)*pr+sq*pj
            frac=(1.d0-sq)*fr+sq
            w(i,j,k,rho_)=den
            w(i,j,k,mom(1))=den*vx
            w(i,j,k,mom(2))=my
            w(i,j,k,mom(3))=den*vz
            w(i,j,k,e_)=press/(eos%gamma-1.d0)+0.5d0*(den*(vx**2+vz**2)+my**2/den)
            w(i,j,k,tracer(1))=den*frac
          end do
        end do
        end do
      end if
    case default
      call mpistop('Only bottom y boundary may be special')
    end select
  end subroutine specialbound_usr
  
  subroutine extra_var_output(ixI^L,ixO^L, w, x, normconv)
    use mod_physics
    use mod_global_parameters
    integer, intent(in)              :: ixI^L,ixO^L
    double precision, intent(in)     :: x(ixI^S, 1:ndim)
    double precision                 :: w(ixI^S, nw+nwauxio), wlocal(ixI^S, 1:nw)
    double precision                 :: normconv(0:nw+nwauxio)

    double precision                 :: pth(ixI^S)
    
    wlocal(ixI^S, 1:nw) = w(ixI^S, 1:nw)

    ! Dimensionless temperature proxy p/rho, NOT Kelvin
    call phys_get_pthermal(wlocal, x, ixI^L,ixO^L, pth)
    w(ixO^S, nw+1) = pth(ixO^S) / w(ixO^S, rho_)
    
  end subroutine extra_var_output
  
  subroutine extra_var_names_output(varnames)
    character(len=*)  :: varnames
    
    varnames = "T"
    
  end subroutine extra_var_names_output
    
end module mod_usr
