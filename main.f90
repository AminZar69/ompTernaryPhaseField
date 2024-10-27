!**********************************************************************

module global_var
    implicit none
    integer :: t
    integer,parameter :: tf   = 3000000
    integer,parameter :: step = 20000
    integer,parameter :: nx = 399
    integer,parameter :: ny = 99
    integer,parameter :: num_threads = 38
    integer,parameter :: ex(0:8) = [0, 1, 0,-1, 0, 1,-1,-1, 1]
    integer,parameter :: ey(0:8) = [0, 0, 1, 0,-1, 1, 1,-1,-1]

    real(8),parameter :: wa(0:8) = [16,4, 4, 4, 4, 1, 1, 1, 1] / 36.d0
    real(8),parameter :: r1 = 40.0d0!l0/8.d0
    real(8),parameter :: r2 = 40.0d0!l0/8.d0
    real(8),parameter :: w  = 4.0d0
    real(8),parameter :: drop_cent_y= 1.
    real(8),parameter :: density(3) = [1000., 500., 1.]
    real(8),parameter :: tau1  = 0.8d0
    real(8),parameter :: tau2  = 0.8d0
    real(8),parameter :: tau3  = 0.8d0
    real(8),parameter :: sigma(3,3) = reshape((/1.5, 1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5/), (/3,3/))
    real(8),parameter :: lambda(3) = [(sigma(2,1) +  &
                                    sigma(3,1) - &
                                    sigma(3,2)), &
                                    (sigma(2,1) + sigma(3,2) - sigma(3,1)), &
                                    (sigma(3,1) + sigma(3,2) - sigma(2,1)) ]
    real(8),parameter :: lambda_t = 3.0d0 / ((1.0d0 / lambda(1)) + &
                                    (1.0d0 / lambda(2)) + (1.0d0 / lambda(3)))
    real(8),parameter :: pi = acos(-1.0d0)
    real(8),parameter :: m = 0.01d0
    real(8),parameter :: w_phi = 1.d0/(0.5d0 + 3.d0*m)
    real(8),parameter :: tetad13 = 140.0d0
    real(8),parameter :: tetad23 = 90.0d0
    real(8),parameter :: teta13 = tetad13 * pi / 180.0d0
    real(8),parameter :: teta23 = tetad23 * pi / 180.0d0
    real(8),parameter :: teta12 = dacos((sigma(3,1) * dcos(teta13) - sigma(3,2) * dcos(teta23))/sigma(2,1))
    real(8),parameter :: gx(0:8)= [0., 2., 0.,-2., 0., 1.,-1.,-1., 1.]
    real(8),parameter :: gy(0:8)= [0., 0., -2., 0., 2., -1.,-1., 1., 1.]

    real(8) :: h(2,0:8,0:nx+1,0:ny+1), g(0:8,0:nx+1,0:ny+1)
    real(8) :: gamma(0:8), heq(2,0:8), geq(0:8), hlph(2,0:8), efh(2,0:8), hlpg(0:8), ef(0:8)
    real(8) :: phi(3,0:nx+1,0:ny+1), p(nx,ny), mu(3,nx,ny), dphi_dx(3,nx,ny), dphi_dy(3,nx,ny), mio(1:3,nx,ny)
    real(8) :: rho(0:nx+1,0:ny+1), ux(nx,ny), uy(nx,ny), ni(3,nx,ny), nj(3,nx,ny)
    real(8) :: ds_dx(nx,ny), ds_dy(nx,ny), abs_ds(nx,ny), ds_dx8(nx,ny), ds_dy8(nx,ny)
    integer :: is_solid_node(0:nx+1,0:ny+1),solid_plt(nx,ny), angle_s8(nx,ny)

end module global_var

!**********************************************************************

program TernaryContactAngleOnFlat
    use omp_lib
    use global_var
    implicit none
    integer :: conjunction_x, i
    real(8) :: t1, t2
    real(8) :: angle_l, angle_r, angle_m
    real(8) :: grad_inside, grad_conj, grad_outside

    do i = 1, 3
        if(lambda(i)<=(lambda_t/2.).or.lambda(i)<0.or.lambda_t<0) then
            print*, "please change the values of surface tension!"
            stop
        end if
    end do

!    call OMP_SET_NUM_THREADS(num_threads)

    t1 = OMP_GET_WTIME()

    call Set_Solid_Geometry
    call Solid_Show
    call Cal_Solid_Gradient
    call Init_Populations


    !--------------------------------------------------------------------------------------------------
    print '(/a7,5a12,7a12/)', 't', 'phi_min', 'phi_max', 'ux_max', 'uy_max', '|u_max|', 'mass_phi', &
    'left_ca', 'middle_ca', 'right_ca', 'grad_inside', 'grad_conj', 'grad_outisde'
    !--------------------------------------------------------------------------------------------------


    do t = 0, tf


        if( mod(t,5000)==0 )then

            call Get_Left_Contact_Angle(angle_l, conjunction_x)
            call Get_Right_Contact_Angle(angle_r)
            call Get_Middle_Contact_Angle(angle_m)
            call Get_Solid_Delta_Phi(conjunction_x, grad_inside, grad_conj, grad_outside)
            call Write_Phi_Profile

            call Write_Contact_Angle_Height_Length

            print '(i7,2f12.6,3e12.3,4f12.2, 3e12.3)', &

              t, minval(phi(3,:,:)), &
                maxval(phi(3,:,:)), &
                maxval(abs(ux)), &
                maxval(abs(uy)), &
                dsqrt( maxval(ux**2+uy**2) ), &
                sum(rho(1:nx,1:ny)), angle_l, &
                angle_m, angle_r, &
                grad_inside, grad_conj, grad_outside

        end if

        if( mod(t,step)==0 )then

            call Write_Results

        end if


        call Collision

        call Set_Periodic_h( h )
        call Set_Periodic_g( g )

        call Propagation_h( h )
        call Propagation_g( g )

        call Interface_Capturing_Properties
        call Hydrodynamics_Properties

    end do
    !
    t2 = omp_get_wtime()
    print '(/a)',' *****************************************'
    print '(a,f12.1)', ' time elapsed:', t2-t1
    print '( a)',' *****************************************'

end

!**********************************************************************

subroutine Init_Populations
    use global_var
    implicit none
    integer :: i
    integer :: x, y
    real(8) :: ga_wa(0:8), ri(2)

    do y = 0, ny+1
        do x = 0, nx+1

            if(is_solid_node(x,y)==0) then

                ri(1) = dsqrt( (x-((nx/2.0d0-2.*r1)))**2.d0 + (y-(1.))**2.d0 )


                ri(2) = dsqrt( (x-(nx/2.0d0+2.*r1))**2.d0 + (y-(1.))**2.d0 )


                if(x<=nint(nx/2.0d0)) then

                    phi(1,x,y) = 0.5d0 + 0.5d0 * tanh(2.0d0*(r1-ri(1))/w)

                else

                    phi(2,x,y) = 0.5d0 + 0.5d0 * tanh(2.0d0*(r2-ri(2))/w)


                end if

                phi(3,x,y) = 1.0d0 - phi(1,x,y) - phi(2,x,y)

            end if

            rho(x,y) = sum(phi(:,x,y) * density(:), 1)

        end do

    end do

    call Set_Periodic_Phi( phi )

    call Get_Chemical_Potential

    call Cal_Phi_Gradient(phi, dphi_dx, dphi_dy)

    call Get_Fluid_Fluid_Normal

    p=0.0
    ux=0.0d0
    uy=0.0d0

    !--- Initializing the interface-capturing distribution functions and the hydrodynamics distribution function with their equilibrium ---!
    do y = 1, ny
        do x = 1, nx

            call Cal_Equilibrium( ux(x,y), uy(x,y), ga_wa )

            gamma(:) = ga_wa(:) + wa(:)


            efh(1,:)  = 2.0d0 * ((4.0d0/w) * dabs(phi(1,x,y) * (1.0d0 - phi(1,x,y))) * (ex(:) * ni(1,x,y) + ey(:) * nj(1,x,y))) &
                - ((4.0d0/w) * dabs(phi(2,x,y) * (1.0d0 - phi(2,x,y))) * (ex(:) * ni(2,x,y) + ey(:) * nj(2,x,y)) &
                + (4.0d0/w) * dabs(phi(3,x,y) * (1.0d0 - phi(3,x,y))) * (ex(:) * ni(3,x,y) + ey(:) * nj(3,x,y)))

            efh(2,:)  = 2.0d0 * ((4.0d0/w) * dabs(phi(2,x,y) * (1.0d0 - phi(2,x,y))) * (ex(:) * ni(2,x,y) + ey(:) * nj(2,x,y))) &
                - ((4.0d0/w) * dabs(phi(1,x,y) * (1.0d0 - phi(1,x,y))) * (ex(:) * ni(1,x,y) + ey(:) * nj(1,x,y)) &
                + (4.0d0/w) * dabs(phi(3,x,y) * (1.0d0 - phi(3,x,y))) * (ex(:) * ni(3,x,y) + ey(:) * nj(3,x,y)))


            do i=1,2


                hlph(i,:) = wa(:) * efh(i,:) / 3.0d0

                h(i,:, x,y) = phi(i,x,y) * gamma(:) - 0.5d0 * hlph(i,:)


            end do


            g(:, x,y) = p(x,y) * wa(:) + ga_wa(:)

        end do
    end do

end

!**********************************************************************

subroutine Cal_Equilibrium( u, v, ga_wa )
    use global_var, only: ex, ey, wa
    implicit none
    real(8), intent(in) :: u, v
    real(8), intent(out) :: ga_wa(0:8)

    real(8) :: u2, eu(0:8)

    u2 = u*u + v*v

    eu(:) = ex(:) * u  + ey(:) * v

    ga_wa(:) = wa(:) * ( eu(:)*(3.d0 + 4.5d0*eu(:)) - 1.5d0*u2 )

end

!**********************************************************************

subroutine Collision
    use global_var
    implicit none
    integer :: i
    integer :: x, y
    real(8) :: fp_x, fp_y, fm_x, fm_y,heav1,heav2, fs_x, fs_y, fx, fy, ga_wa(0:8), tauu


    call Get_Fluid_Fluid_Normal

    !$omp parallel do &
    !$omp private(fp_x, fp_y, fm_x, fm_y,fs_x,fs_y,tauu,geq, ga_wa, heav1,heav2, fx, fy, efh, hlph,heq,ef,hlpg, gamma)&
    !$omp schedule(runtime)
    do y = 1, ny
        do x = 1, nx

            !--- Interface-capturing distribution functions---!

            call Cal_Equilibrium( ux(x,y), uy(x,y), ga_wa )

            gamma(:) = ga_wa(:) + wa(:)


            efh(1,:)  = 2.0d0 * ((4.0d0/w) * dabs(phi(1,x,y) * (1.0d0 - phi(1,x,y))) * (ex(:) * ni(1,x,y) + ey(:) * nj(1,x,y))) &
                - ((4.0d0/w) * dabs(phi(2,x,y) * (1.0d0 - phi(2,x,y))) * (ex(:) * ni(2,x,y) + ey(:) * nj(2,x,y)) &
                + (4.0d0/w) * dabs(phi(3,x,y) * (1.0d0 - phi(3,x,y))) * (ex(:) * ni(3,x,y) + ey(:) * nj(3,x,y)))

            efh(2,:)  = 2.0d0 * ((4.0d0/w) * dabs(phi(2,x,y) * (1.0d0 - phi(2,x,y))) * (ex(:) * ni(2,x,y) + ey(:) * nj(2,x,y))) &
                - ((4.0d0/w) * dabs(phi(1,x,y) * (1.0d0 - phi(1,x,y))) * (ex(:) * ni(1,x,y) + ey(:) * nj(1,x,y)) &
                + (4.0d0/w) * dabs(phi(3,x,y) * (1.0d0 - phi(3,x,y))) * (ex(:) * ni(3,x,y) + ey(:) * nj(3,x,y)))


            do i=1,2

                hlph(i,:) = wa(:) * efh(i,:) / 3.0d0

                heq(i,:) = phi(i,x,y)*gamma(:) - 0.5d0 * hlph(i,:)

                if (is_solid_node(x,y)==0) then

                    h(i,:, x,y) = h(i,:, x,y) * (1.d0-w_phi) + heq(i,:) * w_phi + hlph(i,:)

                end if

            end do


            !--- Hydrodynamics distribution function---!

            fp_x = - p(x,y) * sum(density(:) * dphi_dx(:,x,y), 1) / 3.0d0
            fp_y = - p(x,y) * sum(density(:) * dphi_dy(:,x,y), 1) / 3.0d0

            call Cal_Heaviside(phi(1,x,y)-0.5d0, heav1)
            call Cal_Heaviside(phi(2,x,y)-0.5d0, heav2)

            tauu = (tau1 - tau3) * heav1 + (tau2 - tau3) * heav2 + tau3

            geq(:) = p(x,y) * wa(:) + ga_wa(:)

            call Cal_Viscous_Force(tauu, dphi_dx(:,x,y), dphi_dy(:,x,y), g(:,x,y)-geq(:), fm_x, fm_y)

            fs_x = sum(mu(:,x,y)*dphi_dx(:,x,y))
            fs_y = sum(mu(:,x,y)*dphi_dy(:,x,y))


            if (is_solid_node(x,y)==0) then
                fx =fs_x + fp_x + fm_x
                fy =fs_y + fp_y + fm_y
            end if

            ef(:) = ex(:) * fx + ey(:) * fy

            hlpg(:) = 3.d0 * wa(:) * ef(:) / rho(x,y)

            geq(:) = p(x,y) * wa(:) + ga_wa(:) - 0.5d0 * hlpg(:)

            if (is_solid_node(x,y)==0) then


                g(:, x,y) = g(:, x,y) * (1-(1.0d0/tauu)) + geq(:) * (1.0d0/tauu) +  hlpg(:)


            end if

        end do
    end do
    !$omp end parallel do

end

!**********************************************************************

subroutine Set_Periodic_g( f )
    use global_var, only: nx, ny
    implicit none
    real(8), intent(inout) :: f(0:8, 0:nx+1,0:ny+1)

    !$omp parallel workshare
    !-- left and right boundaries
    f(:,  0  ,:)=f(:, nx,:) !periodic
    f(:, nx+1,:)=f(:, 1 ,:) !periodic

    !-- bottom and top boundaries
    f(:, :, 0  )=f(:, :,ny) !periodic
    f(:, :,ny+1)=f(:, :,1 ) !periodic
    !$omp end parallel workshare
end

!**********************************************************************

subroutine Set_Periodic_h( f )
    use global_var, only: nx, ny
    implicit none
    integer :: o
    real(8), intent(inout) :: f(2,0:8, 0:nx+1,0:ny+1)

    do o = 1,2

        !-- left and right boundaries
        f(o,:,  0  ,:)=f(o,:, nx,:) !periodic
        f(o,:, nx+1,:)=f(o,:, 1 ,:) !periodic

        !-- bottom and top boundaries
        f(o,:, :, 0  )=f(o,:, :,ny) !periodic
        f(o,:, :,ny+1)=f(o,:, :,1 ) !periodic

    end do

end

!**********************************************************************

subroutine Set_Periodic_Phi( a )
    use global_var, only: nx, ny
    implicit none
    real(8),intent(inout):: a(3,0:nx+1,0:ny+1)

    !$omp parallel workshare

    !-- bottom and top boundaries
    a(:,:,  0 ) = a(:,:,ny  )   ! periodic
    a(:,:,ny+1) = a(:,:, 1  )   ! periodic

    !-- left and right boundaries
    a(:,  0 ,:) = a(:,nx  ,:)   ! periodic
    a(:,nx+1,:) = a(:, 1  ,:)   ! periodic

    !$omp end parallel workshare

end

!**********************************************************************

subroutine Set_Periodic_Rho( a )
    use global_var, only: nx, ny
    implicit none
    real(8),intent(inout):: a(0:nx+1,0:ny+1)

    !$omp parallel workshare
    !-- bottom and top boundaries
    a(:,  0 ) = a(:,ny  )   ! periodic
    a(:,ny+1) = a(:, 1  )   ! periodic

    !-- left and right boundaries
    a(  0 ,:) = a(nx  ,:)   ! periodic
    a(nx+1,:) = a( 1  ,:)   ! periodic
    !$omp end parallel workshare

end

!**********************************************************************

subroutine Propagation_g( f )
    use global_var
    implicit none
    real(8), intent(inout) :: f(0:8, 0:nx+1,0:ny+1)
    integer :: i
    integer :: x, y
    real(8) :: fnew(8,nx,ny)

    !$omp parallel
    !$omp do schedule(runtime)
    do y=1,ny
        do x=1,nx
            do i = 1, 8

                fnew(i,x,y) = f(i,x-ex(i),y-ey(i))

                !--- Half-way bounce-back ---!
                if (is_solid_node(x-1,y)==1.or. &
                    is_solid_node(x+1,y)==1.or. &
                    is_solid_node(x,y+1)==1.or. &
                    is_solid_node(x,y-1)==1.or. &
                    is_solid_node(x-1,y+1)==1.or. &
                    is_solid_node(x+1,y+1)==1.or. &
                    is_solid_node(x+1,y-1)==1.or. &
                    is_solid_node(x-1,y-1)==1) then

                    if (is_solid_node(x-1,y)==1) then
                        fnew(1,x,y) = f(3,x,y)
                    end if
                    if (is_solid_node(x+1,y)==1) then
                        fnew(3,x,y) = f(1,x,y)
                    end if
                    if (is_solid_node(x,y-1)==1) then
                        fnew(2,x,y) = f(4,x,y)
                    end if
                    if (is_solid_node(x,y+1)==1) then
                        fnew(4,x,y) = f(2,x,y)
                    end if

                    if (is_solid_node(x+1,y+1)==1) then
                        fnew(7,x,y) = f(5,x,y)
                    end if

                    if (is_solid_node(x-1,y-1)==1) then
                        fnew(5,x,y) = f(7,x,y)
                    end if

                    if (is_solid_node(x-1,y+1)==1) then
                        fnew(8,x,y) = f(6,x,y)
                    end if
                    if (is_solid_node(x+1,y-1)==1) then
                        fnew(6,x,y) = f(8,x,y)
                    end if
                end if

            end do
        end do
    end do
    !$omp end do
    !$omp workshare
    f(1:8,1:nx,1:ny) = fnew
    !$omp end workshare
    !$omp end parallel

end

!**********************************************************************

subroutine Propagation_h( f )
    use global_var
    implicit none
    real(8), intent(inout) :: f(2,0:8, 0:nx+1,0:ny+1)

    integer :: i
    integer :: x, y
    integer :: o
    real(8) :: fnew(2,8,nx,ny)

    !$omp parallel
    !$omp do schedule(runtime)
    do y=1,ny
        do x=1,nx
            do i = 1, 8
                do o=1,2

                    fnew(o,i,x,y) = f(o,i,x-ex(i),y-ey(i))

                    !--- Half-way bounce-back ---!
                    if (is_solid_node(x-1,y)==1.or. &
                        is_solid_node(x+1,y)==1.or. &
                        is_solid_node(x,y+1)==1.or. &
                        is_solid_node(x,y-1)==1.or. &
                        is_solid_node(x-1,y+1)==1.or. &
                        is_solid_node(x+1,y+1)==1.or. &
                        is_solid_node(x+1,y-1)==1.or. &
                        is_solid_node(x-1,y-1)==1) then

                        if (is_solid_node(x-1,y)==1) then
                            fnew(o,1,x,y) = f(o,3,x,y)
                        end if
                        if (is_solid_node(x+1,y)==1) then
                            fnew(o,3,x,y) = f(o,1,x,y)
                        end if
                        if (is_solid_node(x,y-1)==1) then
                            fnew(o,2,x,y) = f(o,4,x,y)
                        end if
                        if (is_solid_node(x,y+1)==1) then
                            fnew(o,4,x,y) = f(o,2,x,y)
                        end if

                        if (is_solid_node(x+1,y+1)==1) then
                            fnew(o,7,x,y) = f(o,5,x,y)
                        end if

                        if (is_solid_node(x-1,y-1)==1) then
                            fnew(o,5,x,y) = f(o,7,x,y)
                        end if

                        if (is_solid_node(x-1,y+1)==1) then
                            fnew(o,8,x,y) = f(o,6,x,y)
                        end if
                        if (is_solid_node(x+1,y-1)==1) then
                            fnew(o,6,x,y) = f(o,8,x,y)
                        end if
                    end if

                end do
            end do

        end do

    end do
    !$omp end do

    !$omp workshare

    f(:,1:8,1:nx,1:ny) = fnew

    !$omp end workshare

    !$omp end parallel

end

!**********************************************************************

subroutine Interface_Capturing_Properties
    use global_var
    implicit none
    integer :: x, y
    integer :: o
    real(8) :: phi_f(3), h1, phi_temp(3,0:nx+1,0:ny+1), rho_temp(0:nx+1,0:ny+1)

    phi_temp=0
    rho_temp=0
    !$omp parallel
    !$omp do private(h1, phi_f)  &
    !$omp schedule(runtime)
    do y=1, ny
        do x=1, nx
            if (is_solid_node(x,y)==0) then

                do o=1,2

                    phi_temp(o,x,y) = sum(h(o,:,x,y))

                end do

                phi_temp(3,x,y) = 1.0d0 - phi_temp(1,x,y) - phi_temp(2,x,y)

            else if(abs_ds(x,y)>0) then


                call Get_Neighbouring_Phi(x,y,phi_f,h1)

!                call Explicit_Wetting_Condition(x,y,phi_f,h1,phi_temp)

                call Implicit_Wetting_Condition(x,y,phi_f,h1,phi_temp)

            end if

            rho_temp(x,y) = sum(phi_temp(:,x,y) * density(:),1)

        end do
    end do
    !$omp end  do

    !$omp workshare
        phi=phi_temp
        rho=rho_temp
    !$omp end  workshare
    !$omp end parallel

end

!**********************************************************************

subroutine Hydrodynamics_Properties
    use global_var
    implicit none
    integer :: x, y
    real(8) :: fp_x, fp_y, fm_x, fm_y, heav1, heav2, fs_x, fs_y, fx, fy, ga_wa(0:8), tauu


    call Set_Periodic_Phi( phi )
    call Get_Chemical_Potential
    call Cal_Phi_Gradient(phi, dphi_dx, dphi_dy)

    !$omp parallel do private(fp_x, fp_y, fm_x, fm_y,fs_x,fs_y,tauu,geq, ga_wa, heav1,heav2, fx, fy) schedule(runtime)
    do y = 1, ny

        do x = 1, nx

            p(x,y) = sum( g(:, x,y), 1)

            fp_x = - p(x,y) * sum(density(:) * dphi_dx(:,x,y), 1) / 3.0d0
            fp_y = - p(x,y) * sum(density(:) * dphi_dy(:,x,y), 1) / 3.0d0

            call Cal_Equilibrium( ux(x,y), uy(x,y), ga_wa )

            geq(:) = p(x,y) * wa(:) + ga_wa(:)

            call Cal_Heaviside(phi(1,x,y)-0.5d0, heav1)
            call Cal_Heaviside(phi(2,x,y)-0.5d0, heav2)

            tauu = (tau1 - tau3) * heav1 + (tau2 - tau3) * heav2 + tau3

            call Cal_Viscous_Force( tauu, dphi_dx(:,x,y), dphi_dy(:,x,y), g(:,x,y)-geq(:), fm_x, fm_y )

            fs_x = sum(mu(:,x,y)*dphi_dx(:,x,y))
            fs_y = sum(mu(:,x,y)*dphi_dy(:,x,y))

            if (is_solid_node(x,y)==0) then
                fx = fs_x + fp_x + fm_x
                fy = fs_y + fp_y + fm_y
            end if

            if(is_solid_node(x,y)==0) then
                ux(x,y) = (g(1,x,y)-g(3,x,y)+g(5,x,y)-g(6,x,y)-g(7,x,y)+g(8,x,y)) + 0.5d0*fx/rho(x,y)
                uy(x,y) = (g(2,x,y)-g(4,x,y)+g(5,x,y)+g(6,x,y)-g(7,x,y)-g(8,x,y)) + 0.5d0*fy/rho(x,y)
            else
                ux(x,y) = 0
                uy(x,y) = 0
            end if

        end do

    end do
    !$omp end parallel do

end

!**********************************************************************

subroutine Get_Chemical_Potential
    use global_var, only:  nx, ny, phi, mu, lambda, lambda_t, w, wa, ex, ey, sigma
    implicit none
    integer :: x, y
    real(8) :: d2_phi(3)

    !$omp parallel do private(d2_phi) schedule(runtime)

    do y = 1, ny
        do x = 1, nx

            d2_phi(:) = ( phi(:,x-1,y-1)+phi(:,x+1,y-1)+phi(:,x-1,y+1)+phi(:,x+1,y+1) &
                +4.0d0*(phi(:,x  ,y-1)+phi(:,x-1,y  )+phi(:,x+1,y  )+phi(:,x  ,y+1)) - 20.0d0*phi(:,x,y) )/6.0d0


            mu(:,x,y) = (12.0d0 / w) * (lambda(:) * phi(:,x,y) * (1.0d0 - phi(:,x,y)) * (1.0d0 - 2.0d0 * phi(:,x,y))  &
                - 2.0d0 * lambda_t * phi(1,x,y) * phi(2,x,y) * (1.0d0 - phi(1,x,y) - phi(2,x,y))) &
                - (3.0d0 / 4.0d0) * w * lambda(:) * d2_phi(:) !- 2.0d0 * 0.03d0*phi(:,x,y)**2.0d0

        end do
    end do

    !$omp end parallel do

end

!**********************************************************************

subroutine Get_Fluid_Fluid_Normal
    use global_var
    implicit none
    integer :: x, y
    integer :: o
    real(8) :: temp
    !$omp parallel do private(temp) schedule(runtime)
    do y = 1, ny

        do x = 1, nx

            do o = 1, 3


                temp = dsqrt( dphi_dx(o,x,y)**2.0d0 + dphi_dy(o,x,y)**2.0d0 )

                if(temp>0) then

                    ni(o,x,y) = dphi_dx(o,x,y) / temp
                    nj(o,x,y) = dphi_dy(o,x,y) / temp

                else
                    temp = temp + 1e-32
                    ni(o,x,y) = dphi_dx(o,x,y) / temp
                    nj(o,x,y) = dphi_dy(o,x,y) / temp

                end if


            end do
        end do
    end do
    !$omp end parallel do

end

!**********************************************************************

subroutine Cal_Phi_Gradient( phi, dphi_dx, dphi_dy)
    use global_var, only: nx, ny, wa, ex, ey
    implicit none

    integer :: x, y
    real(8),intent(in) :: phi(1:3,0:nx+1,0:ny+1)
    real(8),intent(out):: dphi_dx(1:3,nx,ny), dphi_dy(1:3,nx,ny)
    !$omp parallel do schedule(runtime)
    do y = 1, ny
        do x = 1, nx

            dphi_dx(:,x,y) = (phi(:,x+1,y) - phi(:,x-1,y))/3.d0 + &
                            ( phi(:,x+1,y-1) + phi(:,x+1,y+1) - phi(:,x-1,y-1) - phi(:,x-1,y+1))/12.d0
            dphi_dy(:,x,y) = (phi(:,x,y+1) - phi(:,x,y-1))/3.d0 + &
                             ( phi(:,x-1,y+1) + phi(:,x+1,y+1) - phi(:,x-1,y-1) - phi(:,x+1,y-1))/12.d0

        end do
    end do
    !$omp end parallel do

end

!**********************************************************************

subroutine Cal_Viscous_Force( tau, dphi_dx, dphi_dy, gneq, fm_x, fm_y )
    implicit none
    real(8),intent(in)  :: tau, dphi_dx(3), dphi_dy(3), gneq(0:8)
    real(8),intent(out) :: fm_x, fm_y

    call Cal_Viscous_Force_BGK( tau, dphi_dx, dphi_dy, gneq, fm_x, fm_y )

end

!**********************************************************************

subroutine Cal_Viscous_Force_BGK( tau, dphi_dx, dphi_dy, gneq, fm_x, fm_y )
    use global_var, only: density
    implicit none
    real(8),intent(in)  :: tau, dphi_dx(3), dphi_dy(3), gneq(0:8)
    real(8),intent(out) :: fm_x, fm_y

    real(8) :: sxx, sxy, syy

    call Cal_Stress_Tensor_BGK( gneq(1:), sxx, sxy, syy )

    fm_x = (0.5d0-tau)/tau * (sxx * sum(density(:)*dphi_dx(:)) + sxy*sum(density(:)*dphi_dy(:)))
    fm_y = (0.5d0-tau)/tau * (sxy*sum(density(:)*dphi_dx(:))+syy*sum(density(:)*dphi_dy(:)))

end

!**********************************************************************

subroutine Cal_Stress_Tensor_BGK( gneq, sxx, sxy, syy )
    use global_var, only: ex, ey
    implicit none
    real(8),intent(in) :: gneq(1:8)
    real(8),intent(out):: sxx, sxy, syy

    sxx = sum( gneq(1:) * ex(1:) * ex(1:) )
    sxy = sum( gneq(1:) * ex(1:) * ey(1:) )
    syy = sum( gneq(1:) * ey(1:) * ey(1:) )

end

!**********************************************************************

subroutine Write_Results
    use global_var

    implicit none
    integer :: x, y

    character(len=18):: filename,form
    integer :: ierror_result

    select case (int(t/step))
        case(0:9)
            write(form,'(a9,i1)') 'xml_000000',int(t/step)
        case(10:99)
            write(form,'(a8,i2)') 'xml_00000',int(t/step)
        case(100:999)
            write(form,'(a7,i3)') 'xml_0000',int(t/step)
        case(1000:9999)
            write(form,'(a6,i4)') 'xml_000',int(t/step)
        case(10000:99999)
            write(form,'(a5,i5)') 'xml_00',int(t/step)
        case(100000:999999)
            write(form,'(a4,i6)') 'xml_0',int(t/step)
        case default
            write(form,'(a9,i1)') 'xml_000000',int(-1)
    end select

    write(filename,'(a10,a4)')trim(form),".vtr"

    open(unit=1,file=filename,status='unknown',iostat=ierror_result)

    write (1,'(a21)') '<?xml version="1.0"?>'
    write (1,'(a)') '<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">'
    write (1,*) '<RectilinearGrid WholeExtent="1 ',nx,' 1 ',ny,' 1 ','1','">'
    write (1,*) '<Piece Extent="1 ',nx,' 1 ',ny,' 1 ','1',' ">'

    write (1,*) '<Coordinates>'
    write (1,'(a)') '<DataArray type="Float32" Name="X" NumberOfComponents="1" format="ascii">'
    do x= 1,nx
        write (1,*) x
    end do
    write (1,'(a)') '</DataArray>'

    write (1,'(a)') '<DataArray type="Float32" name="Y" NumberOfComponents="1" format="ascii">'
    do y= 1,ny

        write (1,*) y

    end do
    write (1,'(a)') '</DataArray>'

    write (1,'(a)') '<DataArray type="Float32" name="Z" NumberOfComponents="1" format="ascii">'

    write (1,*) 1

    write (1,'(a)') '</DataArray>'
    write (1,*) '</Coordinates>'

    write (1,'(a)') '<PointData>'
    write (1,'(a)') '<DataArray type="Float32" Name="Density" NumberOfComponents="1" format="ascii">'
    do y = 1, ny
        do x = 1, nx
            write (1,*) rho(x,y)
        end do
    end do
    write (1,'(a)') '</DataArray>'

    write (1,'(a)') '<DataArray type="Float32" Name="Phi1" NumberOfComponents="1" format="ascii">'
    do y= 1, ny
        do x= 1, nx
            write (1,*) phi(1,x,y)
        end do
    end do
    write (1,'(a)') '</DataArray>'

    write (1,'(a)') '<DataArray type="Float32" Name="Phi2" NumberOfComponents="1" format="ascii">'
    do y= 1, ny
        do x= 1, nx
            write (1,*) phi(2,x,y)
        end do
    end do
    write (1,'(a)') '</DataArray>'

    write (1,'(a)') '<DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">'
    do y= 1, ny
        do x= 1, nx
            write (1,*) ux(x,y), uy(x,y),0
        end do
    end do
    write (1,'(a)') '</DataArray>'
    write (1,'(a)') '<DataArray type="Float32" Name="Pressure" NumberOfComponents="1" format="ascii">'
    do y= 1, ny
        do x= 1, nx
            write (1,*) p(x,y)*rho(x,y)/3.0d0
        end do
    end do
    write (1,'(a)') '</DataArray>'

    write (1,'(a)') '<DataArray type="Float32" Name="Solid" NumberOfComponents="1" format="ascii">'
    do y= 1, ny
        do x= 1, nx
            write (1,*) solid_plt(x,y)
        end do
    end do
    write (1,'(a)') '</DataArray>'

    write (1,'(a)') '</PointData>'
    write (1,'(a)') '</Piece>'
    write (1,'(a)') '</RectilinearGrid>'
    write (1,'(a)') '</VTKFile>'
    close(1)

end

!**********************************************************************

subroutine Set_Solid_Geometry
    use global_var
    implicit none
    integer :: x

    do x=0,nx+1
        is_solid_node(x,1)=1
        is_solid_node(x,0)=1
        is_solid_node(x,ny+1)=1
        is_solid_node(x,ny)=1
    end do


end subroutine Set_Solid_Geometry

!**********************************************************************

subroutine Solid_Show
    use global_var
    implicit none
    integer :: x, y

    do y=1,ny
        do x=1,nx

            if(is_solid_node(x,y)==1.and.abs_ds(x,y)==0) then

                solid_plt(x,y) = is_solid_node(x,y)

            end if

        end do
    end do
end

!**********************************************************************

subroutine Cal_Heaviside(dum, heav)
    use global_var
    implicit none
    real(8), intent(in) :: dum
    real(8), intent(out) :: heav

    heav = 0.5d0 + sign(0.5d0,dum)

end

!**********************************************************************

subroutine Get_Neighbouring_Phi(x, y, phi_f, h1)
    use global_var, only : phi, angle_s8, pi
    implicit none
    integer, intent(in) :: x, y
    real(8), intent(out) :: phi_f(3), h1
    real(8) :: wei1, wei2


    if (angle_s8(x,y)==45) then ! angle=45
        h1=dsqrt(2.0d0)
        phi_f(:) = phi(:,x-1,y+1)
    end if

    if (angle_s8(x,y)==135) then ! angle=135
        h1=dsqrt(2.0d0)
        phi_f(:) = phi(:,x+1,y+1)
    end if

    if (angle_s8(x,y)==(-135)) then ! angle=-135
        h1=dsqrt(2.0d0)
        phi_f(:) = phi(:,x+1,y-1)
    end if

    if (angle_s8(x,y)==(-45)) then ! angle=-45
        h1=dsqrt(2.0d0)
        phi_f(:) = phi(:,x-1,y-1)
    end if

    if (angle_s8(x,y)==0) then ! angle=0
        h1=1.0d0
        phi_f(:) = phi(:,x-1,y)
    end if

    if (angle_s8(x,y)==90) then ! angle=90
        h1=1.0d0
        phi_f(:) = phi(:,x,y+1)
    end if

    if (angle_s8(x,y)==-90) then ! angle=-90
        h1=1.0d0
        phi_f(:) = phi(:,x,y-1)
    end if

    if (angle_s8(x,y)==180.or.angle_s8(x,y)==(-180)) then ! angle=180
        h1=1.0d0
        phi_f(:) = phi(:,x+1,y)
    end if

    if (0<angle_s8(x,y).and.angle_s8(x,y)<45) then !0-45
        wei1=abs(angle_s8(x,y)-0.0d0)/45.0d0
        wei2=abs(angle_s8(x,y)-45.0d0)/45.0d0
        h1=1.0d0/cos(angle_s8(x,y)*pi/180.0d0)
        phi_f(:) = wei1 * phi(:,x-1,y) + wei2 * phi(:,x-1,y+1)
    end if

    if (45<angle_s8(x,y).and.angle_s8(x,y)<90) then !45-90
        wei1=abs(angle_s8(x,y)-45.0d0)/45.0d0
        wei2=abs(angle_s8(x,y)-90.0d0)/45.0d0
        h1=1.0d0/cos(abs(angle_s8(x,y)-90.0d0)*pi/180.0d0)
        phi_f(:) = wei1 * phi(:,x-1,y+1) + wei2 * phi(:,x,y+1)
    end if

    if (90<angle_s8(x,y).and.angle_s8(x,y)<135) then !90-135
        wei1=abs(angle_s8(x,y)-90.0d0)/45.0d0
        wei2=abs(angle_s8(x,y)-135.0d0)/45.0d0
        h1=1.0d0/cos(abs(90.0d0-angle_s8(x,y))*pi/180.0d0)
        phi_f(:) = wei1 * phi(:,x,y+1) + wei2 * phi(:,x+1,y+1)
    end if

    if (135<angle_s8(x,y).and.angle_s8(x,y)<180) then !135-180
        wei1=abs(angle_s8(x,y)-180.0d0)/45.0d0
        wei2=abs(angle_s8(x,y)-135.0d0)/45.0d0
        h1=1.0d0/cos(abs(180.0d0-angle_s8(x,y))*pi/180.0d0)
        phi_f(:) = wei1 * phi(:,x+1,y) + wei2 * phi(:,x+1,y+1)
    end if

    if ((-180)<angle_s8(x,y).and.angle_s8(x,y)<(-135)) then !(-180)-(-135)
        wei1=abs(angle_s8(x,y)-(-180.0d0))/45.0d0
        wei2=abs(angle_s8(x,y)-(-135.0d0))/45.0d0
        h1=1.0d0/cos(abs(-180.0d0-angle_s8(x,y))*pi/180.0d0)
        phi_f(:) = wei1 * phi(:,x+1,y) + wei2 * phi(:,x+1,y-1)
    end if

    if ((-135)<angle_s8(x,y).and.angle_s8(x,y)<(-90)) then !(-135)-(-90)
        wei1=abs(angle_s8(x,y)-(-90.0d0))/45.0d0
        wei2=abs(angle_s8(x,y)-(-135.0d0))/45.0d0
        h1=1.0d0/cos(abs(-90.0d0-angle_s8(x,y))*pi/180.0d0)
        phi_f(:) = wei1 * phi(:,x,y-1) + wei2 * phi(:,x+1,y-1)
    end if

    if ((-90)<angle_s8(x,y).and.angle_s8(x,y)<(-45)) then !(-90)-(-45)
        wei1=abs(angle_s8(x,y)-(-90.0d0))/45.0d0
        wei2=abs(angle_s8(x,y)-(-45.0d0))/45.0d0
        h1=1.0d0/cos(abs(-90.0d0-angle_s8(x,y))*pi/180.0d0)
        phi_f(:) = wei1 * phi(:,x,y-1) + wei2 * phi(:,x-1,y-1)
    end if

    if ((-45)<angle_s8(x,y).and.angle_s8(x,y)<0) then !(0)-(-45)
        wei1=abs(angle_s8(x,y)-0.0d0)/45.0d0
        wei2=abs(angle_s8(x,y)-(-45.0d0))/45.0d0
        h1=1.0d0/cos(abs(-angle_s8(x,y))*pi/180.0d0)
        phi_f(:) = wei1 * phi(:,x-1,y) + wei2 * phi(:,x-1,y-1)
    end if

end

!**********************************************************************
!Gaussian elimination with partial pivoting and scaling
subroutine Gauss_Elimination(a,b,x,n)
    implicit none
    integer :: n
    real(8) :: a(n,n), b(n), x(n)
    real(8) :: s(n)
    real(8) :: c, pivot, store
    integer i, j, k, l

    ! step 1: begin forward elimination
    do k=1, n-1

        ! step 2: "scaling"
        ! s(i) will have the largest element from row i
        do i=k,n                       ! loop over rows
            s(i) = 0.0
            do j=k,n                    ! loop over elements of row i
                s(i) = max(s(i),abs(a(i,j)))
            end do
        end do

        ! step 3: "pivoting 1"
        ! find a row with the largest pivoting element
        pivot = abs(a(k,k)/s(k))
        l = k
        do j=k+1,n
            if(abs(a(j,k)/s(j)) > pivot) then
                pivot = abs(a(j,k)/s(j))
                l = j
            end if
        end do

        ! check if the system has a sigular matrix
        if(abs(pivot) < 1.0e-15) then
            write(*,*) ' the matrix is singular '
            return
        end if

        ! step 4: "pivoting 2" interchange rows k and l (if needed)
        if (l /= k) then
            do j=k,n
                store = a(k,j)
                a(k,j) = a(l,j)
                a(l,j) = store
            end do
            store = b(k)
            b(k) = b(l)
            b(l) = store
        end if

        ! step 5: the elimination (after scaling and pivoting)
        do i=k+1,n
            c=a(i,k)/a(k,k)
            a(i,k) = 0.0
            b(i)=b(i)- c*b(k)
            do j=k+1,n
                a(i,j) = a(i,j)-c*a(k,j)
            end do
        end do
    end do

    ! step 6: back substitution
    x(n) = b(n)/a(n,n)
    do i=n-1,1,-1
        c=0.0
        do j=i+1,n
            c= c + a(i,j)*x(j)
        end do
        x(i) = (b(i)- c)/a(i,i)
    end do

end subroutine Gauss_Elimination

!**********************************************************************

subroutine Implicit_Wetting_Condition(x,y,phi_f,h1,phi_temp)
    use global_var, only : phi, teta12, teta13, teta23, w, nx, ny
    implicit none
    integer, intent(in) :: x, y
    real(8), intent(in) :: phi_f(3), h1
    real(8), intent(out) :: phi_temp(3,0:nx+1,0:ny+1)
    real(8) :: j(3,3), f(3), phi_s(3), deltax(3), phi_s_new(3), err

    phi_s(:) = phi(:,x,y) ! initial guess
    err=10.

    do while(err>1e-8)

        f(1) = (1.0d0/w)*(-cos(teta13)* &
            ((phi_f(1) + phi_s(1))*(phi_f(3)+phi_s(3))) - &
            cos(teta12)*((phi_f(1) + phi_s(1))*(phi_f(2)+phi_s(2)))) + &
            ((phi_s(1)-phi_f(1))/h1)
        f(2) = (1.0d0/w)*(cos(teta12)* &
            ((phi_f(1) + phi_s(1))*(phi_f(2)+phi_s(2))) - &
            cos(teta23)*((phi_f(2) + phi_s(2))*(phi_f(3)+phi_s(3)))) + &
            ((phi_s(2)-phi_f(2))/h1)
        f(3) = (1.0d0/w)*(cos(teta23)* &
            ((phi_f(2) + phi_s(2))*(phi_f(3)+phi_s(3))) + &
            cos(teta13)*((phi_f(1) + phi_s(1))*(phi_f(3)+phi_s(3)))) + &
            ((phi_s(3)-phi_f(3))/h1)

        j(1,1) = (1.0d0/h1) + (1.0d0/w) * (-cos(teta13) * (phi_f(3) + phi_s(3)) - cos(teta12) * (phi_f(2) + phi_s(2)))
        j(1,2) = (1.0d0/w) * (-cos(teta12) * (phi_f(1) + phi_s(1)))
        j(1,3) = (1.0d0/w) * (-cos(teta13) * (phi_f(1) + phi_s(1)))
        j(2,1) = (1.0d0/w) * (cos(teta12) * (phi_f(2) + phi_s(2)))
        j(2,2) = (1.0d0/h1) + (1.0d0/w) * (cos(teta12) * (phi_f(1) + phi_s(1)) - cos(teta23) * (phi_f(3) + phi_s(3)))
        j(2,3) = (1.0d0/w) * (-cos(teta23) * (phi_f(2) + phi_s(2)))
        j(3,1) = (1.0d0/w) * (cos(teta13) * (phi_f(3) + phi_s(3)))
        j(3,2) = (1.0d0/w) * (cos(teta23) * (phi_f(3) + phi_s(3)))
        j(3,3) = (1.0d0/h1) + (1.0d0/w) * (cos(teta23) * (phi_f(2) + phi_s(2)) + cos(teta13) * (phi_f(1) + phi_s(1)))

        call Gauss_Elimination(j,-f,deltax,3)

        phi_s_new = deltax + phi_s

        !err = abs(phi_s_new(3) + phi_s_new(2) + phi_s_new(1) - (phi_s(3) + phi_s(2) + phi_s(1)))
        !err = abs(phi_s_new(1) - phi_s(1))
        err = dsqrt(f(1)**2.0d0 + f(2)**2.0d0 + f(3)**2.0d0)

        phi_s = phi_s_new
    end do

    phi_temp(:,x,y) = phi_s(:)

end

!**********************************************************************

subroutine Cal_Solid_Gradient
    use global_var, only : is_solid_node, nx, ny, ds_dx, ds_dy, gx, gy, abs_ds, pi, ds_dx8, ds_dy8, angle_s8
    implicit none
    integer :: x, y
    do y=1,ny
        do x=1,nx
            if (is_solid_node(x,y)==1) then
                !---4th-order isotropy---!
                ds_dx(x,y) = is_solid_node(x,y)*gx(0) + &
                    is_solid_node(x+1,y)*gx(1) + &
                    is_solid_node(x,y+1)*gx(2) + &
                    is_solid_node(x-1,y)*gx(3) + &
                    is_solid_node(x,y-1)*gx(4) + &
                    is_solid_node(x+1,y+1)*gx(5) + &
                    is_solid_node(x-1,y+1)*gx(6) + &
                    is_solid_node(x-1,y-1)*gx(7) + &
                    is_solid_node(x+1,y-1)*gx(8)
                ds_dy(x,y) = is_solid_node(x,y)*gy(0) + &
                    is_solid_node(x+1,y)*gy(1) + &
                    is_solid_node(x,y+1)*gy(2) + &
                    is_solid_node(x-1,y)*gy(3) + &
                    is_solid_node(x,y-1)*gy(4) + &
                    is_solid_node(x+1,y+1)*gy(5) + &
                    is_solid_node(x-1,y+1)*gy(6) + &
                    is_solid_node(x-1,y-1)*gy(7) + &
                    is_solid_node(x+1,y-1)*gy(8)
                abs_ds(x,y)= dsqrt(ds_dx(x,y)**2+ds_dy(x,y)**2)

                !---8th-order isotropy---!
                ds_dx8(x,y) = (4.0d0/21.0d0)*(is_solid_node(x+1,y)-is_solid_node(x-1,y)) + &
                    (4.0d0/45.0d0)*(is_solid_node(x+1,y+1)-is_solid_node(x-1,y+1)) + &
                    (4.0d0/45.0d0)*(is_solid_node(x+1,y-1)-is_solid_node(x-1,y-1)) + &
                    (2.0d0/60.0d0)*(is_solid_node(x+2,y)-is_solid_node(x-2,y)) + &
                    (4.0d0/315.0d0)*(is_solid_node(x+2,y+1)-is_solid_node(x-2,y+1)) + &
                    (4.0d0/315.0d0)*(is_solid_node(x+2,y-1)-is_solid_node(x-2,y-1)) + &
                    (2.0d0/315.0d0)*(is_solid_node(x+1,y+2)-is_solid_node(x-1,y+2)) + &
                    (2.0d0/315.0d0)*(is_solid_node(x+1,y-2)-is_solid_node(x-1,y-2)) + &
                    (2.0d0/5040.0d0)*(is_solid_node(x+2,y+2)-is_solid_node(x-2,y+2))+ &
                    (2.0d0/5040.0d0)*(is_solid_node(x+2,y-2)-is_solid_node(x-2,y-2))

                ds_dy8(x,y) = (4.0d0/21.0d0)*(is_solid_node(x,y+1)-is_solid_node(x,y-1)) + &
                    (4.0d0/45.0d0)*(is_solid_node(x+1,y+1)-is_solid_node(x+1,y-1)) + &
                    (4.0d0/45.0d0)*(is_solid_node(x-1,y+1)-is_solid_node(x-1,y-1)) + &
                    (2.0d0/60.0d0)*(is_solid_node(x,y+2)-is_solid_node(x,y-2)) + &
                    (2.0d0/315.0d0)*(is_solid_node(x+2,y+1)-is_solid_node(x+2,y-1)) + &
                    (2.0d0/315.0d0)*(is_solid_node(x-2,y+1)-is_solid_node(x-2,y-1)) + &
                    (4.0d0/315.0d0)*(is_solid_node(x+1,y+2)-is_solid_node(x+1,y-2)) + &
                    (4.0d0/315.0d0)*(is_solid_node(x-1,y+2)-is_solid_node(x-1,y-2)) + &
                    (2.0d0/5040.0d0)*(is_solid_node(x+2,y+2)-is_solid_node(x+2,y-2))+ &
                    (2.0d0/5040.0d0)*(is_solid_node(x-2,y+2)-is_solid_node(x-2,y-2))

                angle_s8(x,y)=nint(atan2(-ds_dy8(x,y),ds_dx8(x,y))*180.0d0/pi)

            else

               angle_s8(x,y) = 90

            end if

        end do
    end do
end

!**********************************************************************

subroutine Explicit_Wetting_Condition(x,y,phi_f,h1,phi_temp)
    use global_var, only : teta12, teta13, teta23, w, ny, nx
    implicit none
    integer, intent(in) :: x, y
    real(8), intent(in) :: phi_f(3), h1
    real(8), intent(out) :: phi_temp(3,0:nx+1,0:ny+1)

    phi_temp(1,x,y) = phi_f(1) - h1 * (4.0d0/w)*(-cos(teta13)*(phi_f(1)*phi_f(3)) - cos(teta12)*(phi_f(1)*phi_f(2)))
    phi_temp(2,x,y) = phi_f(2) - h1 * (4.0d0/w)*(cos(teta12)*(phi_f(1)*phi_f(2)) - cos(teta23)*(phi_f(2)*phi_f(3)))
    phi_temp(3,x,y) = phi_f(3) - h1 * (4.0d0/w)*(cos(teta23)*(phi_f(2)*phi_f(3)) + cos(teta13)*(phi_f(1)*phi_f(3)))

end

!**********************************************************************

subroutine Get_Left_Contact_Angle(angle_l, conjunction_x)
    use global_var
    implicit none
    real(8), intent(out) :: angle_l
    integer, intent(out) :: conjunction_x
    integer :: p_i(4),p_j(4), x, y
    real (8) :: xp(2),yp(2)
    real (8) :: phi1,phi2,phi3,phi4

    logical :: is_x1_found,is_x2_found,is_inside

    xp=20
    yp=20
    p_i = 20
    p_j = 20

    ! is it inside drop?
    is_inside = .false.
    ! is point 1 or 3 found?
    is_x1_found = .false.
    is_x2_found = .false.

    y = 1

    is_inside = .true.
    do x=1, nx
        if (is_inside.and.0.5<phi(1,x,y).and.0.5>phi(1,x-1,y)) then
            xp(1) = real(0.5-phi(1,x-1,y))/(phi(1,x,y)-phi(1,x-1,y)) + (x-1)
            yp(1) = y
            p_i(1) = x
            conjunction_x = x
            p_j(1) = y
            is_x1_found = .true.
            exit
        endif
    end do



    p_i(2) = p_i(1) - 1
    p_j(2) = p_j(1)

    p_i(3) = p_i(2)
    p_j(3) = p_j(2) + 1

    p_i(4) = p_i(1)
    p_j(4) = p_j(3)


    phi1 = phi(1,p_i(1),p_j(1))
    phi2 = phi(1,p_i(2),p_j(2))
    phi3 = phi(1,p_i(3),p_j(3))
    phi4 = phi(1,p_i(4),p_j(4))
    !print*, p_i(4),p_j(4)

    is_x2_found = .true.
    !   p2 = vapor, p3 = liquid
    if (0.5.le.phi3.and.0.5.gt.phi2) then
        xp(2) = p_i(2)
        yp(2) = (0.5-phi2)/(phi3-phi2) + p_j(2)
    !   p4 = vapor, p1 = liquid
    elseif(0.5.le.phi1.and.0.5.gt.phi4) then
        xp(2) = p_i(1)
        yp(2) = (0.5-phi1)/(phi4-phi1) + p_j(1)
    !   p3 = vapor, p4 = liquid
    elseif(0.5.le.phi4.and.0.5.gt.phi3) then
        xp(2) =(0.5-phi3)/(phi4-phi3) + p_i(3)
        yp(2) = p_j(3)
    else
        is_x2_found = .false.
    endif

    ! all points are found.
    if (is_x1_found.and.is_x2_found) then
        angle_l = datan( (yp(2)-yp(1)) / (xp(2)-xp(1)) ) * 180/pi
        if(angle_l<0) angle_l = angle_l+180
    else
        angle_l = -1.d0
    endif
!print*, angle_l
end subroutine Get_Left_Contact_Angle

!**********************************************************************

subroutine Get_Right_Contact_Angle(angle_r)
    use global_var
    implicit none
    real(8), intent(out) :: angle_r
    integer :: p_i(4),p_j(4), x, y
    real (8) :: xp(2),yp(2)
    real (8) :: phi1,phi2,phi3,phi4

    logical :: is_x1_found,is_x2_found,is_inside
    x=20
    y=20

    p_i = 20
    p_j = 20

    ! is it inside drop?
    is_inside = .false.
    ! is point 1 or 3 found?
    is_x1_found = .false.
    is_x2_found = .false.

    ! sweeping right-bottom of the drop
    ! j=1 is the wall.
    y = 1
    is_inside = .true.
    do x=1, nx
        if (is_inside.and.0.5.le.phi(2,x,y).and.0.5.gt.phi(2,x+1,y)) then
            xp(1) = real(0.5-phi(2,x,y))/(phi(2,x+1,y)-phi(2,x,y)) + x
            yp(1) = y
            p_i(1) = x
            p_j(1) = y
            is_x1_found = .true.
            exit
        endif
    end do

    p_i(2) = p_i(1)+1
    p_j(2) = p_j(1)

    p_i(3) = p_i(2)
    p_j(3) = p_j(2) + 1

    p_i(4) = p_i(1)
    p_j(4) = p_j(3)


    phi1 = phi(2,p_i(1),p_j(1))
    phi2 = phi(2,p_i(2),p_j(2))
    phi3 = phi(2,p_i(3),p_j(3))
    phi4 = phi(2,p_i(4),p_j(4))

    is_x2_found = .true.
    !   p2 = vapor, p3 = liquid
    if (0.5.le.phi3.and.0.5.gt.phi2) then
        xp(2) = p_i(2)
        yp(2) = (0.5-phi2)/(phi3-phi2) + p_j(2)
    !   p4 = vapor, p1 = liquid
    elseif(0.5.le.phi1.and.0.5.gt.phi4) then
        xp(2) = p_i(1)
        yp(2) = (0.5-phi1)/(phi4-phi1) + p_j(1)
    !   p3 = vapor, p4 = liquid
    elseif(0.5.le.phi4.and.0.5.gt.phi3) then
        xp(2) =(0.5-phi4)/(phi3-phi4) + p_i(4)
        yp(2) = p_j(4)
    else
        is_x2_found = .false.
    endif

    ! all points are found.
    if (is_x1_found.and.is_x2_found) then
        angle_r = atan( (yp(2)-yp(1)) / (xp(1)-xp(2)) ) * 180/pi
        if(angle_r<0) angle_r = angle_r+180
    else
        angle_r = -1.d0
    endif

end subroutine Get_Right_Contact_Angle

!**********************************************************************

subroutine Get_Middle_Contact_Angle(angle_m)
    use global_var
    implicit none
    real(8), intent(out) :: angle_m
    integer :: p_i(4),p_j(4), x, y
    real (8) :: xp(2),yp(2)
    real (8) :: phi1,phi2,phi3,phi4

    logical :: is_x1_found,is_x2_found,is_inside
    x=20
    y=20
    p_i = 20
    p_j = 20

    ! is it inside drop?
    is_inside = .false.
    ! is point 1 or 3 found?
    is_x1_found = .false.
    is_x2_found = .false.

    ! sweeping right-bottom of the drop
    ! j=1 is the wall.
    y = 1
    is_inside = .true.

    do x=1, nx
        if (is_inside.and.0.5.le.phi(1,x,y).and.0.5.gt.phi(1,x+1,y)) then
            xp(1) = real(0.5-phi(1,x,y))/(phi(1,x+1,y)-phi(1,x,y)) + x
            yp(1) = y
            p_i(1) = x
            p_j(1) = y
            is_x1_found = .true.
            exit
        endif
    end do

    p_i(2) = p_i(1)+1
    p_j(2) = p_j(1)

    p_i(3) = p_i(2)
    p_j(3) = p_j(2) + 1

    p_i(4) = p_i(1)
    p_j(4) = p_j(3)


    phi1 = phi(1,p_i(1),p_j(1))
    phi2 = phi(1,p_i(2),p_j(2))
    phi3 = phi(1,p_i(3),p_j(3))
    phi4 = phi(1,p_i(4),p_j(4))

    is_x2_found = .true.
    !   p2 = vapor, p3 = liquid
    if (0.5.le.phi3.and.0.5.gt.phi2) then
        xp(2) = p_i(2)
        yp(2) = (0.5-phi2)/(phi3-phi2) + p_j(2)
    !   p4 = vapor, p1 = liquid
    elseif(0.5.le.phi1.and.0.5.gt.phi4) then
        xp(2) = p_i(1)
        yp(2) = (0.5-phi1)/(phi4-phi1) + p_j(1)
    !   p3 = vapor, p4 = liquid
    elseif(0.5.le.phi4.and.0.5.gt.phi3) then
        xp(2) =(0.5-phi4)/(phi3-phi4) + p_i(4)
        yp(2) = p_j(4)
    else
        is_x2_found = .false.
    endif

    ! all points are found.
    if (is_x1_found.and.is_x2_found) then
        angle_m = atan( (yp(2)-yp(1)) / (xp(1)-xp(2)) ) * 180/pi
        if(angle_m<0) angle_m = angle_m+180
    else
        angle_m = -1.d0
    endif

end subroutine Get_Middle_Contact_Angle

!**********************************************************************

subroutine Get_Solid_Delta_Phi(conjunction_x, grad_inside, grad_conj, grad_outside)

    use global_var
    implicit none
    integer, intent(in) :: conjunction_x
    real(8), intent(out) :: grad_inside, grad_conj, grad_outside
    !integer :: insidex, outside x


    grad_conj = phi(1,conjunction_x,2) - phi(1,conjunction_x,1)
    grad_inside = phi(1,int(nx/2.0d0-2.*r1),2) - phi(1,int(nx/2.0d0-2.*r1),1)
    grad_outside = phi(1,10,2) - phi(1,10,1)

end subroutine

!**********************************************************************

subroutine Write_Phi_Profile
    use global_var
    implicit none
    integer :: x

    open(2,file = 'x1andx2.csv')
    open(3,file = 'delta_phi.csv')

    do x = 1, nx

        write(2,'(i6.1,a2,f12.6,a2,f12.6)') x, ',', phi(1,x,1), ',', phi(1,x,2)
        write(3,'(i6.1,a2,e12.3)') x, ',', phi(1,x,2)-phi(1,x,1)

    end do

    close(2)
    close(3)

end

!**********************************************************************
! Calculates the contact angle based on the shape of droplet i.e., height and length
subroutine Write_Contact_Angle_Height_Length
    use global_var
    implicit none

    integer :: h_max1, length_max1, x, y
    integer :: ss1, h_temp1, length_temp1, west, east
    real(8) :: contact_angel

    open(5,file = 'length_height_contact_angle.csv')

    ss1 = 0
    h_temp1 = 0
    h_max1 = 0
    length_max1=0
    west=0
    east=0
    length_temp1=0

    x = int(nx/2.0-2.*r1)
    do y=ny,1,-1
        if (is_solid_node(x,y)==0) then
            if (phi(1,x,y)>=0.5.and.phi(1,x,y+1)<0.5) then
                ss1=y
            end if
        end if

        h_temp1 =  abs(ss1)

        if (h_max1<h_temp1) then
            h_max1=h_temp1
        end if
    end do

    y=1
    do x = 1, nx
        if(phi(1,x,y)>=0.5.and.phi(1,x-1,y)<0.5) then
            west = x
        end if
        if(phi(1,x,y)>=0.5.and.phi(1,x+1,y)<0.5) then
            east = x
        end if

         length_temp1 = abs(west-east)

        if (length_max1<length_temp1) then
            length_max1=length_temp1
        end if
    end do

    contact_angel = 2.0d0 * datan(2.0d0*h_max1/length_max1) * 180.0d0/pi

    write(5,'(i6.1,a2,i6.1,a2,f12.2)') length_max1, ',', h_max1, ',', contact_angel
    close(5)

end
!**********************************************************************
