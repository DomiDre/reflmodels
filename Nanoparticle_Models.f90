module nanoparticle_models
use math
implicit none
contains
    !Sphere CS NP
    double precision function cs_sphere_core_a_fraction(x, p, Np)
        double precision, intent(in) :: x
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np
        
        double precision :: x0
        double precision :: density, core_radius, shell_thickness
        
        double precision :: particle_size, xshift
        
        x0 = p(1)
        density = p(2)
        core_radius = p(3)
        shell_thickness = p(4)
        
        particle_size = core_radius + shell_thickness
        xshift = x - particle_size - x0
        
        if (xshift < -core_radius) then
            cs_sphere_core_a_fraction = 0d0
        else if (xshift  < core_radius) then
            cs_sphere_core_a_fraction = ((core_radius/particle_size)**2 -&
                                         (xshift/particle_size)**2)
        else
            cs_sphere_core_a_fraction = 0d0
        end if
    end function cs_sphere_core_a_fraction
    
    double precision function cs_sphere_shell_a_fraction(x, p, Np)
        double precision, intent(in) :: x
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np
        
        double precision :: x0
        double precision :: density, core_radius, shell_thickness
        
        double precision :: particle_size, xshift
        
        x0 = p(1)
        density = p(2)
        core_radius = p(3)
        shell_thickness = p(4)
        
        particle_size = core_radius + shell_thickness
        xshift = x - particle_size - x0
        
        if (xshift < -particle_size) then
            cs_sphere_shell_a_fraction = 0d0
        else if (xshift < -core_radius) then
            cs_sphere_shell_a_fraction =&
                                (1d0- (xshift/particle_size)**2 )
        else if (xshift  < core_radius) then
            cs_sphere_shell_a_fraction =&
                                (1d0 - (core_radius/particle_size)**2)
        else if (xshift < particle_size) then 
            cs_sphere_shell_a_fraction =&
                                (1d0- (xshift/particle_size)**2 )
        else
            cs_sphere_shell_a_fraction = 0d0
        end if
    end function cs_sphere_shell_a_fraction
    
    subroutine sld_coreshell_sphere(x, x0,&
                                    density, core_radius, shell_thickness,&
                                    sig_radius, sig_x0,&
                                    sld_core, sld_shell, sld_matrix,&
                                    Nx, sld)
    double precision, intent(in), dimension(Nx) :: x
    double precision, intent(in) :: x0, density, core_radius, shell_thickness
    double precision, intent(in) :: sig_radius, sig_x0
    complex*16, intent(in) :: sld_core, sld_shell, sld_matrix
    integer, intent(in) :: Nx
    complex*16, intent(out), dimension(Nx) :: sld
    
    integer, parameter :: Np=4
    double precision, dimension(Nx) :: area_core, area_shell
    double precision, dimension(Np) :: p
    double precision :: Rmin, Rmax
    double precision :: x0min, x0max

    integer :: ix
    
    p = (/x0, density, core_radius, shell_thickness/)
    call get_cutoff_gaussian(x0, sig_x0, x0min, x0max) 
    call get_cutoff_lognormal(core_radius, sig_radius, Rmin, Rmax) 

    !$omp parallel
    !$omp do
    do ix=1, Nx
        call integrate_two_size_distributions(x(ix), p, Np, &
                        1, x0min, x0max, sig_x0, &
                        3, Rmin, Rmax, sig_radius, &
                        cs_sphere_core_a_fraction, gaussian, lognormal, area_core(ix))
        call integrate_two_size_distributions(x(ix), p, Np, &
                        1, x0min, x0max, sig_x0, &
                        3, Rmin, Rmax, sig_radius, &
                        cs_sphere_shell_a_fraction, gaussian, lognormal, area_shell(ix))
    end do
    !$omp end do
    !$omp end parallel

    sld = sld_core*area_core*density + sld_shell*area_shell*density +&
          sld_matrix*(1d0 - area_core*density - area_shell*density)
    end subroutine sld_coreshell_sphere

    subroutine sld_immersed_coreshell_sphere(x, x0,&
                                    density, core_radius, shell_thickness,&
                                    immersion_depth,&
                                    sig_radius,&
                                    sld_core, sld_shell, sld_matrix,&
                                    Nx, sld)
    double precision, intent(in), dimension(Nx) :: x
    double precision, intent(in) :: x0, density, core_radius, shell_thickness
    double precision, intent(in) :: immersion_depth
    double precision, intent(in) :: sig_radius
    complex*16, intent(in) :: sld_core, sld_shell, sld_matrix
    integer, intent(in) :: Nx
    complex*16, intent(out), dimension(Nx) :: sld
    
    integer, parameter :: Np=4
    double precision, dimension(Nx) :: area_core, area_shell, sld_immersion
    double precision, dimension(Np) :: p
    double precision :: Rmin, Rmax
    double precision :: x0min, x0max

    integer :: ix
    
    p = (/x0, density, core_radius, shell_thickness/)
    call get_cutoff_lognormal(core_radius, sig_radius, Rmin, Rmax) 

    !$omp parallel
    !$omp do
    do ix=1, Nx
        call integrate_size_distribution(x(ix), p, Np, &
                        3, Rmin, Rmax, sig_radius, &
                        cs_sphere_core_a_fraction, lognormal, area_core(ix))
        call integrate_size_distribution(x(ix), p, Np, &
                        3, Rmin, Rmax, sig_radius, &
                        cs_sphere_shell_a_fraction, lognormal, area_shell(ix))
          
        if (x(ix) < immersion_depth) then
            sld_immersion(ix) = sld_matrix
        else
            sld_immersion(ix) = 0d0
        end if
    end do
    !$omp end do
    !$omp end parallel


    sld = sld_core*area_core*density + sld_shell*area_shell*density +&
          sld_immersion*(1d0 - area_core*density - area_shell*density)
    end subroutine sld_immersed_coreshell_sphere
    
    ! Cube CS NP
    double precision function cs_cube_core_a_fraction(x, p, Np)
        double precision, intent(in) :: x
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np
        
        double precision :: x0
        double precision ::density, core_edge_length, shell_thickness
        
        double precision :: particle_size, xshift
        
        x0 = p(1)
        density = p(2)
        core_edge_length = p(3)
        shell_thickness = p(4)
        
        particle_size = core_edge_length + 2*shell_thickness
        xshift = x - x0
        
        if (xshift < shell_thickness) then
            cs_cube_core_a_fraction = 0d0
        else if (xshift  < core_edge_length+shell_thickness) then
            cs_cube_core_a_fraction =&
                                 ((core_edge_length/particle_size)**2)*density
        else
            cs_cube_core_a_fraction = 0d0
        end if
    end function cs_cube_core_a_fraction
    
    double precision function cs_cube_shell_a_fraction(x, p, Np)
        double precision, intent(in) :: x
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np
        
        double precision :: x0
        double precision :: density, core_edge_length, shell_thickness
        
        double precision :: particle_size, xshift
        
        x0 = p(1)
        density = p(2)
        core_edge_length = p(3)
        shell_thickness = p(4)
        
        particle_size = core_edge_length + 2*shell_thickness
        xshift = x - x0
        
        if (xshift < 0d0) then
            cs_cube_shell_a_fraction = 0d0
        else if (xshift < shell_thickness) then
            cs_cube_shell_a_fraction = density
        else if (xshift  < shell_thickness+core_edge_length) then
            cs_cube_shell_a_fraction =&
                            (1d0 - (core_edge_length/particle_size)**2)*density
        else if (xshift < particle_size) then 
            cs_cube_shell_a_fraction = density
        else
            cs_cube_shell_a_fraction = 0d0
        end if
    end function cs_cube_shell_a_fraction
    
    subroutine sld_coreshell_cube(x, x0,&
                                    density, core_edgelength, shell_thickness,&
                                    sig_edgelength, sig_x0,&
                                    sld_core, sld_shell, sld_matrix,&
                                    Nx, sld)
    double precision, intent(in), dimension(Nx) :: x
    double precision, intent(in) :: x0, density, core_edgelength, shell_thickness
    double precision, intent(in) :: sig_edgelength, sig_x0
    complex*16, intent(in) :: sld_core, sld_shell, sld_matrix
    integer, intent(in) :: Nx
    complex*16, intent(out), dimension(Nx) :: sld
    
    integer, parameter :: Np=4
    double precision, dimension(Nx) :: area_core, area_shell
    double precision, dimension(Np) :: p
    double precision :: edgemin, edgemax
    double precision :: x0min, x0max

    integer :: ix
    
    p = (/x0, density, core_edgelength, shell_thickness/)
    call get_cutoff_gaussian(x0, sig_x0, x0min, x0max) 
    call get_cutoff_lognormal(core_edgelength, sig_edgelength, edgemin, edgemax) 

    !$omp parallel
    !$omp do
    do ix=1, Nx
        call integrate_two_size_distributions(x(ix), p, Np, &
                        1, x0min, x0max, sig_x0, &
                        3, edgemin, edgemax, sig_edgelength, &
                        cs_cube_core_a_fraction, gaussian, lognormal, area_core(ix))
        call integrate_two_size_distributions(x(ix), p, Np, &
                        1, x0min, x0max, sig_x0, &
                        3, edgemin, edgemax, sig_edgelength, &
                        cs_cube_shell_a_fraction, gaussian, lognormal, area_shell(ix))
    end do
    !$omp end do
    !$omp end parallel

    sld = sld_core*area_core + sld_shell*area_shell +&
          sld_matrix*(1d0 - area_core - area_shell)
    end subroutine sld_coreshell_cube


    subroutine sld_colloidal_coreshell_cubes_45deg_period(x, x0,&
                                a, d,&
                                sld_core, sld_shell, sld_matrix,&
                                N_periods, Nx, sld)

    double precision, intent(in), dimension(Nx) :: x
    double precision, intent(in) :: x0
    double precision, intent(in) :: a, d
    complex*16, intent(in) :: sld_core, sld_shell, sld_matrix
    integer, intent(in) :: N_periods, Nx
    complex*16, intent(out), dimension(Nx) :: sld
    
    double precision, dimension(Nx) :: area_core, area_shell
    integer :: ix, iper
    double precision :: z
    double precision :: A_ref,sq2d, sq2a, sq2a2, diagonal_center, diagonal_length
    double precision :: half_cube_diagonal
    !a: core edge length
    !d: shell thickness
    
    sq2a = sq2*a
    sq2a2 = sq2*a**2
    sq2d = sq2*d
    half_cube_diagonal = 5d-1*sq2a
    diagonal_center = half_cube_diagonal + sq2d
    diagonal_length = sq2*(a+2d0*d)
    area_shell = 0d0
    area_core = 0d0
    A_ref = 1d0/(2d0*(a+2*d)**2)

!    !$omp parallel
!    !$omp do
    do ix=1, Nx
        z = x(ix) - x0
        ! bottom cube layer
        if (z < 0 .or. z > (N_periods+1)*diagonal_length) continue
        
        if (z < sq2d) then
            continue
        else if (z < diagonal_center) then
            area_core(ix) = area_core(ix) + 2d0*a*(z - sq2d)*A_ref
        else
            continue
        end if
        z = z - diagonal_center
        do iper=0, N_periods-1
            ! outer 4 cubes area fraction:
            if (z < 0d0) then
                continue
            else if (z < half_cube_diagonal) then
                area_core(ix) = area_core(ix) + (sq2a2 - 2d0*a*z) * A_ref
            else if (z < diagonal_center+sq2d) then
                continue
            else if (z <= diagonal_length) then
                area_core(ix) = area_core(ix) + 2d0*a*(z - (diagonal_center+sq2d))*A_ref
            else
                continue
            end if
            
            ! center cube area fraction:
            if (z < sq2d) then
                continue
            else if (z < diagonal_center) then
                area_core(ix) = area_core(ix) + 2d0*a*(z - sq2d)*A_ref
            else if (z < diagonal_length - sq2d) then
                area_core(ix) = area_core(ix) + (sq2a2 - 2d0*a*(z - diagonal_center))*A_ref
            else
                continue
            end if
            z = z - diagonal_length
        end do
!         top cube layer
        if (z < 0d0) then
            continue
        else if (z < half_cube_diagonal) then
            area_core(ix) = area_core(ix) + (sq2a2 - 2d0*a*z)*A_ref
        else
            continue
        end if
    end do
!    !$omp end do
!    !$omp end parallel

    sld = sld_core*area_core + sld_shell*area_shell + sld_matrix*(1d0 - area_core - area_shell)
    end subroutine sld_colloidal_coreshell_cubes_45deg_period

end module nanoparticle_models
