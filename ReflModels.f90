module models
use math
implicit none
contains
    subroutine parrat_amplitude(q, sld, sigma, thickness, NLayers, Nq, refl)
        double precision, intent(in), dimension(Nq) :: q
        complex*16, intent(in), dimension(NLayers) :: sld
        double precision, intent(in), dimension(NLayers) :: sigma, thickness
        integer, intent(in):: NLayers, Nq
        complex*16, intent(out), dimension(Nq) :: refl

        integer :: iq, ir
        
        double precision :: k0, k02
        complex*16 :: ksub, kj, knext, rnext, p2, refl_val

        !$omp parallel
        !$omp do
        do iq=1,Nq
            k0 = q(iq)/2d0
            k02= k0**2
            
            ksub = cdsqrt(k02 - 4d0*pi*sld(1))
            kj = cdsqrt(k02 - 4d0*pi*sld(2))
            refl_val = (kj-ksub)/(kj+ksub) * cdexp(-2d0*ksub*kj*sigma(1)**2)
            do ir=2, NLayers-1
                knext = cdsqrt(k02 - 4d0*pi*sld(ir+1))
                rnext = (knext-kj)/(knext+kj)*cdexp(-2d0*knext*kj*sigma(ir)**2)
                p2 = cdexp(2d0*ci*thickness(ir)*kj)
                refl_val = (rnext + refl_val*p2) / ( 1 + rnext*refl_val*p2)
                kj = knext
            end do
            refl(iq) = refl_val
        end do
        !$omp end do
        !$omp end parallel
    end subroutine parrat_amplitude
    
    subroutine parrat(q, sld, sigma, thickness, NLayers, Nq, refl)
        double precision, intent(in), dimension(Nq) :: q
        complex*16, intent(in), dimension(NLayers) :: sld
        double precision, intent(in), dimension(NLayers) :: sigma, thickness
        integer, intent(in):: NLayers, Nq
        double precision, intent(out), dimension(Nq) :: refl

        complex*16, dimension(Nq) :: refl_ampl
        
        call parrat_amplitude(q, sld, sigma, thickness, NLayers, Nq, refl_ampl)
        refl = zabs(refl_ampl)**2
    end subroutine parrat

    subroutine add_layer_to_sld(x, x0, sld_pre, sld_after, rough,&
                                Nx, roughsld_in, roughsld_out)
        double precision, dimension(Nx), intent(in) :: x
        double precision, intent(in) :: x0
        complex*16, intent(in) :: sld_pre, sld_after
        double precision, intent(in) :: rough
        integer, intent(in) :: Nx
        complex*16, dimension(Nx), intent(in) :: roughsld_in
        complex*16, dimension(Nx), intent(out) :: roughsld_out
        
        complex*16 :: diffsld
        integer :: j
        diffsld = sld_after - sld_pre
        roughsld_out = roughsld_in
        if (rough == 0d0) then 
            do j=1, Nx
                if (x(j) < x0) then
                    cycle
                else
                    roughsld_out(j) = roughsld_in(j) + diffsld
                end if
            end do
        else
            do j=1, Nx
                roughsld_out(j) = roughsld_in(j) +&
                        diffsld * 0.5d0*(1d0+erf((x(j)-x0)/(sq2*rough)))
            end do
        end if

    end subroutine add_layer_to_sld
    
    subroutine roughsld_thick_layers(x, complex_sld, sigma, thickness, &
            Nx, Nthickness, rough_sld)
        double precision, intent(in), dimension(Nx) :: x
        complex*16, intent(in), dimension(Nthickness) :: complex_sld
        double precision, intent(in), dimension(Nthickness) :: thickness,&
                                                               sigma
        integer, intent(in) :: Nx, Nthickness
        complex*16, intent(out), dimension(Nx) :: rough_sld
        
        integer :: i,j
        complex*16 :: diffsld
        double precision :: rough
        double precision :: x0, xlast
        
        !substrate
        rough_sld = complex_sld(1)
        x0 = x(1)+thickness(1)
        !sequentially go layer by layer
        do i=1, Nthickness-1
            call add_layer_to_sld(x, x0, complex_sld(i), complex_sld(i+1),&
                            sigma(i),&
                            Nx, rough_sld, rough_sld)
            x0 = x0 + thickness(i+1)
        end do
    end subroutine roughsld_thick_layers

    subroutine sld_and_rough_from_thicknesses(thickness, sld, sigma, x,&
                        Nx, Nthickness, sld_on_x, sigma_on_x)
        double precision, dimension(Nthickness), intent(in) :: thickness, sld, &
                                                            sigma
        double precision, dimension(Nx), intent(in) :: x
        integer, intent(in) :: Nx, Nthickness
        double precision, dimension(Nx), intent(out) :: sld_on_x, sigma_on_x
        
        integer :: i,j
        double precision :: z0, zlast
        
        !substrate
        do j=1, Nx
            if (x(j) < 0d0) then
                sld_on_x(j) = sld(1)
                sigma_on_x(j) = sigma(1)
            else
                exit
            end if
        end do
        
        !layer by layer
        z0 = thickness(2)
        zlast = 0d0
        do i=2, Nthickness
            do j=1, Nx
                if (x(j) < zlast) then
                    cycle
                else if (x(j) < z0) then
                    sld_on_x(j) = sld(i)
                    sigma_on_x(j) = sigma(i)
                else
                    exit
                end if
            end do
            zlast = z0
            z0 = z0 + thickness(i+1)
        end do
    end subroutine sld_and_rough_from_thicknesses
    
    subroutine get_k(sld, k0, k)
        complex*16, intent(in) :: sld
        double precision, intent(in) :: k0
        complex*16, intent(out) :: k
        
        k = cdsqrt(k0**2 - 4d0*pi*sld)
    end subroutine get_k

    subroutine scattering_matrix_fresnel(k1, k2, sigma, S)
        complex*16, intent(in) :: k1, k2
        double precision, intent(in) :: sigma
        !f2py double precision, intent(in), optional:: sigma=0
        complex*16, dimension(2,2), intent(out) :: S
        
        complex*16 :: r, t
        
        r = (k1-k2) / (k1+k2) * cdexp(-2*sigma**2*k1*k2)
        t = 2*cdsqrt(k1*k2) / (k1+k2) * cdexp(5d-1*sigma**2*(k1-k2)**2)
        S(1,1) = r
        S(1,2) = t
        S(2,1) = t
        S(2,2) = -r
    end subroutine scattering_matrix_fresnel

    subroutine scattering_matrix_translation(k1, thickness, S)
        complex*16, intent(in) :: k1
        double precision, intent(in) :: thickness
        complex*16, dimension(2,2), intent(out) :: S
        
        S(1,1) = 0d0
        S(1,2) = cdexp(-ci*thickness*k1)
        S(2,1) = cdexp(-ci*thickness*k1)
        S(2,2) = 0d0
    end subroutine scattering_matrix_translation
    
    subroutine connect_scattering_matrices(S1, S2, S)
        complex*16, dimension(2,2), intent(in) :: S1, S2
        complex*16, dimension(2,2), intent(out) :: S
        
        complex*16 :: nominator
        nominator = 1d0/(1d0 - S1(2,2)*S2(1,1))
        S(1,1) = S1(1,1) + S1(1,2)*S2(1,1)*S1(2,1)*nominator
        S(2,2) = S2(2,2) + S2(2,1)*S1(2,2)*S2(1,2)*nominator
        S(2,1) = S1(2,1)*S2(2,1)*nominator
        S(1,2) = S1(1,2)*S2(1,2)*nominator
    end subroutine connect_scattering_matrices
    
    subroutine scattering_matrix_algorithm(q, sld, thickness, roughness, &
                                           NLayers, S)
        double precision, intent(in) :: q
        complex*16, dimension(NLayers), intent(in) :: sld
        double precision, dimension(NLayers), intent(in) :: thickness,roughness
        integer, intent(in):: NLayers
        complex*16, dimension(2, 2), intent(out) :: S
        
        integer :: ir
        double precision :: k0, k02
        complex*16 :: ksub, kj, knext
        
        complex*16, dimension(2,2) :: Snext, Slast
        

        if (NLayers  == 0) then
            print *, "Error: Pass non empty sld array to scattering matrix algorithm."
            return
        end if
        
        k0 = q/2d0
        k02= k0**2
        ksub = cdsqrt(k02 - 4d0*pi*sld(1))

        if (NLayers  == 1) then
            call scattering_matrix_fresnel(ksub, cdsqrt(k02+0*ci), roughness(1), S)
            return
        end if
        
        !Nlayers > 1:
        kj = cdsqrt(k02 - 4d0*pi*sld(2))
        call scattering_matrix_fresnel(ksub, kj, roughness(1), S)

        do ir=2, NLayers-1
            ! translate and connect
            call scattering_matrix_translation(kj, thickness(ir), Snext)
            call connect_scattering_matrices(S, Snext, S)

            !calc next and connect
            knext = cdsqrt(k02 - 4d0*pi*sld(ir+1))
            call scattering_matrix_fresnel(kj, knext, roughness(ir), Snext)
            call connect_scattering_matrices(S, Snext, S)
            kj = knext
        end do
    end subroutine scattering_matrix_algorithm
    
    recursive subroutine recursive_scatt_adding(S1, N, S)
        complex*16, dimension(2,2), intent(in) :: S1
        integer, intent(in) :: N
        complex*16, dimension(2,2), intent(out) :: S
        
        complex*16, dimension(2,2) :: Shelp


        if (N > 1) then
            call connect_scattering_matrices(S1, S1, Shelp)
            if(mod(N,2) == 0) then
                call recursive_scatt_adding(Shelp, N/2, S)
            else
                call recursive_scatt_adding(Shelp, (N-1)/2, S)
                call connect_scattering_matrices(S1, S, S)
            end if
        else
            S = S1
        end if
    end subroutine
    
    subroutine scattering_matrix_periodic_algorithm(q, &
            sld_sub, sld_per, sld_end,&
            thickness_sub, thickness_per, thickness_end,&
            roughness_sub, roughness_per, roughness_end,&
            Nperiods, &
            NLayers_sub, Nlayers_per, Nlayers_end, S)
            
        double precision, intent(in) :: q
        complex*16, dimension(NLayers_sub), intent(in) :: sld_sub
        complex*16, dimension(NLayers_per), intent(in) :: sld_per
        complex*16, dimension(NLayers_end), intent(in) :: sld_end
        double precision, dimension(NLayers_sub), intent(in) :: thickness_sub,&
                                                                roughness_sub
        double precision, dimension(NLayers_per), intent(in) :: thickness_per,&
                                                                roughness_per
        double precision, dimension(NLayers_end), intent(in) :: thickness_end,&
                                                                roughness_end
        integer, intent(in) :: Nperiods
        integer, intent(in):: NLayers_sub, Nlayers_per, Nlayers_end
        complex*16, dimension(2, 2), intent(out) :: S
        
        integer :: iperiod, ir
        double precision :: k0, k02
        complex*16 :: ksub1, ksub0, kper1, kper0, kend1, kj, knext
        
        complex*16, dimension(2,2) :: S_sub, S_per, S_con_per, S_end, S_help,&
                                      S_help2
        

        k0 = q/2d0
        k02= k0**2
        ksub1 = cdsqrt(k02 - 4d0*pi*sld_sub(1))
        ksub0 = cdsqrt(k02 - 4d0*pi*sld_sub(NLayers_sub))
        kper1 = cdsqrt(k02 - 4d0*pi*sld_per(1))
        kper0 = cdsqrt(k02 - 4d0*pi*sld_per(NLayers_per))
        kend1 = cdsqrt(k02 - 4d0*pi*sld_end(1))
        if (Nlayers_sub == 0 .OR. NLayers_per == 0 .OR. Nlayers_end == 0) then
            print *, "Error: Substrate, Periodic or End SLD is empty. Check it!"
            return
        end if
        
        ! calc Scattering matrix for substrate
        if (NLayers_sub  == 1) then
            ! just connect substrate to periodic layers
            call scattering_matrix_fresnel(ksub1, kper1,&
                                           roughness_sub(1), S_sub)
        else
            ! first calc scattering matrix for substrate
            call scattering_matrix_algorithm(q, sld_sub, thickness_sub,&
                                             roughness_sub, &
                                             NLayers_sub, S_sub)
            ! do translation of last layer in substrate
            call scattering_matrix_translation(ksub0,&
                                               thickness_sub(Nlayers_sub),&
                                               S_help)
            call connect_scattering_matrices(S_sub, S_help, S_sub)
            
            ! connect with first layer of periodic layer
            call scattering_matrix_fresnel(ksub0, kper1,&
                                            roughness_sub(NLayers_sub), S_help)
            call connect_scattering_matrices(S_sub, S_help, S_sub)
        end if
        
        ! calc scattering matrix for periodic cell
        !first translation of first element
        call scattering_matrix_translation(kper1, thickness_per(1), S_per)
        if (Nlayers_per > 1) then
            call scattering_matrix_algorithm(q, sld_per, thickness_per,&
                                             roughness_per, &
                                             NLayers_per, S_help)
            call connect_scattering_matrices(S_per, S_help, S_per)

            ! do translation of last layer in periodic layer
            call scattering_matrix_translation(kper0, &
                                               thickness_per(Nlayers_per),&
                                               S_help)
            call connect_scattering_matrices(S_per, S_help, S_per)
            
            !calculate matrix to connect two periodic layers
            call scattering_matrix_fresnel(kper0, kper1, &
                                         roughness_per(Nlayers_per), S_con_per)
        else
            ! for Nlayers == 1, connecting matrix aquivalent to 0 translation.
            call scattering_matrix_translation(kper1, 0d0, S_con_per)
        end if
        
        ! calculate scattering matrix for end
        ! first connect last element of periodic layer to first layer of end
        call scattering_matrix_fresnel(kper0, kend1,&
                                        roughness_per(Nlayers_per), S_end)
        if (Nlayers_end > 1) then
            ! do translation of first layer in last layer
            call scattering_matrix_translation(kend1, &
                                               thickness_end(1),&
                                               S_help)
            call connect_scattering_matrices(S_end, S_help, S_end)
                                               
            call scattering_matrix_algorithm(q, sld_end, thickness_end,&
                                             roughness_end, &
                                             NLayers_end, S_help)
            call connect_scattering_matrices(S_end, S_help, S_end)
        end if
        
        ! now connect everything
        call connect_scattering_matrices(S_sub, S_per, S)
        
        call connect_scattering_matrices(S_con_per, S_per, S_help)
        if (Nperiods > 1) then
            call recursive_scatt_adding(S_help, Nperiods-1, S_help2)
            call connect_scattering_matrices(S, S_help2, S)
        end if
        call connect_scattering_matrices(S, S_end, S)
    end subroutine scattering_matrix_periodic_algorithm

    subroutine periodic_model(q, sld_sub, sld_per, sld_end,&
                                thickness_sub, thickness_per, thickness_end,&
                                roughness_sub, roughness_per, roughness_end,&
                                Nperiods,&
                                Nq, NLayers_sub, Nlayers_per, Nlayers_end,&
                                refl)
        double precision, intent(in), dimension(Nq) :: q
        complex*16, dimension(NLayers_sub), intent(in) :: sld_sub
        complex*16, dimension(NLayers_per), intent(in) :: sld_per
        complex*16, dimension(NLayers_end), intent(in) :: sld_end
        double precision, dimension(NLayers_sub), intent(in) :: thickness_sub,&
                                                                roughness_sub
        double precision, dimension(NLayers_per), intent(in) :: thickness_per,&
                                                                roughness_per
        double precision, dimension(NLayers_end), intent(in) :: thickness_end,&
                                                                roughness_end
        integer, intent(in) :: Nperiods
        integer, intent(in):: Nq, NLayers_sub, Nlayers_per, Nlayers_end
        double precision, intent(out), dimension(Nq) :: refl

        complex*16, dimension(2,2) :: S
        integer :: iq
        !$omp parallel private(S)
        !$omp do
        do iq=1,Nq
            call scattering_matrix_periodic_algorithm(q(iq), &
                                sld_sub, sld_per, sld_end,&
                                thickness_sub, thickness_per, thickness_end,&
                                roughness_sub, roughness_per, roughness_end,&
                                Nperiods, &
                                NLayers_sub, Nlayers_per, Nlayers_end, S)
            refl(iq) = abs(S(2, 2))**2
        end do
        !$omp end do
        !$omp end parallel
    end subroutine periodic_model
    
    subroutine add_layers_shifted(x, sld_profile1, sld_profile2, &
                                  x0_profile2, Nx, sld_res)
            double precision, dimension(Nx), intent(in) :: x
            complex*16, dimension(Nx), intent(in) :: sld_profile1, sld_profile2
            double precision, intent(in) :: x0_profile2
            integer, intent(in) :: Nx
            complex*16, dimension(Nx), intent(out) :: sld_res
            
            integer :: ix, ix0
            ix0 = 1
            do ix=1, Nx
                if(x(ix) < x0_profile2) then
                    sld_res(ix) = sld_profile1(ix)
                    ix0 = ix0 + 1
                else
                    sld_res(ix) = sld_profile2(ix-ix0+1) !- sld_profile2(1)
                end if
            end do
    end subroutine add_layers_shifted
    
    subroutine periodic_sld(x,&
                        sld_sub, sld_per, sld_end,&
                        thickness_sub, thickness_per, thickness_end,&
                        roughness_sub, roughness_per, roughness_end,&
                        Nperiods,&
                        Nx, NLayers_sub, Nlayers_per, Nlayers_end,&
                        sld)
        double precision, intent(in), dimension(Nx) :: x
        complex*16, dimension(NLayers_sub), intent(in) :: sld_sub
        complex*16, dimension(NLayers_per), intent(in) :: sld_per
        complex*16, dimension(NLayers_end), intent(in) :: sld_end
        double precision, dimension(NLayers_sub), intent(in) :: thickness_sub,&
                                                                roughness_sub
        double precision, dimension(NLayers_per), intent(in) :: thickness_per,&
                                                                roughness_per
        double precision, dimension(NLayers_end), intent(in) :: thickness_end,&
                                                                roughness_end
        integer, intent(in) :: Nperiods
        integer, intent(in):: Nx, NLayers_sub, Nlayers_per, Nlayers_end
        complex*16, intent(out), dimension(Nx) :: sld
        
        integer :: i, iperiod
        complex*16, dimension(NLayers_sub + Nperiods*NLayers_per+ Nlayers_end)&
                                 :: combined_sld
        double precision, dimension(NLayers_sub + Nperiods*NLayers_per+ Nlayers_end)&
                                 :: combined_roughness, combined_thickness
        
        do i=1, Nlayers_sub
            combined_sld(i) = sld_sub(i)
            combined_roughness(i) = roughness_sub(i)
            combined_thickness(i) = thickness_sub(i)
        end do
        
        do iperiod=1, Nperiods
            do i=1, Nlayers_per
                combined_sld(Nlayers_sub + (iperiod-1)*Nlayers_per + i) =&
                             sld_per(i)
                combined_roughness(Nlayers_sub + (iperiod-1)*Nlayers_per + i) =&
                             roughness_per(i)
                combined_thickness(Nlayers_sub + (iperiod-1)*Nlayers_per + i) =&
                             thickness_per(i)
            end do
        end do
        
        do i=1, Nlayers_end
            combined_sld(Nlayers_sub + Nperiods*Nlayers_per + i) =&
                              sld_end(i)
            combined_roughness(Nlayers_sub + Nperiods*Nlayers_per + i) =&
                              roughness_end(i)
            combined_thickness(Nlayers_sub + Nperiods*Nlayers_per + i) =&
                              thickness_end(i)
        end do
        
        call roughsld_thick_layers(x, combined_sld, combined_roughness,&
                    combined_thickness, &
                    Nx, NLayers_sub + Nperiods*NLayers_per+Nlayers_end, sld)
    end subroutine periodic_sld
    
    
    
    subroutine coreshell_sphere_area_fraction(x, x0,&
                                    density, core_radius, shell_thickness,&
                                    Nx,&
                                    area_core, area_shell)
    double precision, intent(in), dimension(Nx) :: x
    double precision, intent(in) :: x0, density, core_radius, shell_thickness
    integer, intent(in) :: Nx
    double precision, intent(out), dimension(Nx) :: area_core, area_shell
    
    double precision :: particle_size, xshift, constshift
    double precision :: area_core_mod, area_shell_mod
    integer :: i
    
    particle_size = core_radius + shell_thickness
    
    constshift = density*(1d0 - (core_radius/particle_size)**2)
    area_core_mod = density*(core_radius/particle_size)**2
    area_shell_mod = density
    do i=1, Nx
        xshift = x(i) - particle_size - x0
        if (xshift < -particle_size) then
            area_core(i) = 0d0
            area_shell(i) = 0d0
        else if (xshift < -core_radius) then
            area_core(i) = 0d0
            area_shell(i) = (1d0- (xshift/particle_size)**2 )*area_shell_mod
        else if (xshift  < core_radius) then
            area_core(i) = (1d0- (xshift/core_radius)**2 )*area_core_mod
            area_shell(i) = constshift
        else if (xshift < particle_size) then 
            area_core(i) = 0d0
            area_shell(i) = (1d0- (xshift/particle_size)**2 )*area_shell_mod
        else
            area_core(i) = 0d0
            area_shell(i) = 0d0
        end if
    end do
    end subroutine coreshell_sphere_area_fraction
    
    subroutine coreshell_sphere_sld(density, core_radius, shell_thickness,&
                                    sld_core, sld_shell, Nx,&
                                    sldmodel)
    double precision, intent(in) :: density, core_radius, shell_thickness
    complex*16, intent(in) :: sld_core, sld_shell
    integer, intent(in) :: Nx
    complex*16, intent(out), dimension(Nx) :: sldmodel
    
    double precision :: particle_size
    double precision :: xshift, sld_resolution
    double precision :: constshift
    complex*16 :: sld_core_mod, sld_shell_mod
    integer :: i
    
    particle_size = core_radius + shell_thickness
    sld_resolution = 2*particle_size/Nx
    
    constshift = density*(1d0 - (core_radius/particle_size)**2)*sld_shell
    sld_core_mod = density*(core_radius/particle_size)**2*sld_core
    sld_shell_mod = density*sld_shell
    do i=1, Nx
        xshift = (i-1)*sld_resolution - particle_size
        if (xshift < -particle_size) then
            sldmodel(i) = 0d0
        else if (xshift < -core_radius) then
            sldmodel(i) = (1d0- (xshift/particle_size)**2 )*sld_shell_mod
        else if (xshift  < core_radius) then
            sldmodel(i) = (1d0- (xshift/core_radius)**2 )*sld_core_mod + constshift
        else if (xshift < particle_size) then 
            sldmodel(i) = (1d0- (xshift/particle_size)**2 )*sld_shell_mod
        else
            sldmodel(i) = 0d0
        end if
    end do
    end subroutine coreshell_sphere_sld

    subroutine N_layers_coated_spheres_part_dist_2(z, z0_shifts, packing_densities,&
                     Rcore, sigRcore, dshell, SLDc, SLDs, SLDbg, result_sld,&
                     Nz, Nlayers, Nsigs, Npts)
        implicit none
        integer, intent(in) :: Nz, Nlayers
        double precision, intent(in), dimension(Nz) :: z
        double precision, intent(in), dimension(Nlayers) :: z0_shifts, packing_densities 
        double precision, intent(in) :: Rcore, sigRcore, dshell, SLDc, SLDs, SLDbg
        integer,intent(in) :: Nsigs, Npts
    !f2py   integer, optional, intent(in):: Nsigs=5
    !f2py   integer, optional, intent(in):: Npts=200
        double precision, intent(out), dimension(Nz) :: result_sld

        double precision, parameter :: pi = 3.141592653589793
        double precision, dimension(Nz) :: core_area, shell_area
        
        integer :: i, iz, iLayer
        double precision :: dSLDc, dSLDs
        double precision :: Rd, Rparticle, area_prefactor, tessel_z, last_z0
        double precision :: zshift
        double precision, dimension(Npts) :: hcore, hshell
        
        double precision :: Rstd, Rcore_min, Rcore_max, Rcore_steps
        double precision :: Rcore_i, Rd_i
        double precision, dimension(Npts) :: p_log
        
        ! Initialize core and shell areas:
        do i=1, Nz
            core_area(i) = 0d0
            shell_area(i) = 0d0
        end do
        
        dSLDc = SLDc - SLDbg
        dSLDs = SLDs - SLDbg
        
        ! Initialize integration parameters
        Rstd = dsqrt(dexp(2*dlog(Rcore) + sigRcore**2)*(dexp(sigRcore**2)-1))
        Rcore_min = Rcore - Nsigs*Rstd
        Rcore_max = Rcore + Nsigs*Rstd
        Rcore_steps = (Rcore_max-Rcore_min)/(Npts-1)
        
        do i=1, Npts
            Rcore_i = Rcore_min + (i-1)*Rcore_steps
            call lognormal(Rcore_i, Rcore, sigRcore, p_log(i))
        end do
        
        Rd = Rcore + 0.5d0*dshell
        Rparticle = Rcore + dshell
        area_prefactor = pi/(2*dsqrt(3d0))
        tessel_z = 2d0*dsqrt(6d0)*Rd/3d0
        
        last_z0 = Rparticle
        do iLayer=1, Nlayers
            last_z0 = last_z0 + z0_shifts(iLayer)
            if (iLayer > 1) last_z0 = last_z0 + tessel_z
            do iz=1, Nz
                zshift = z(iz) - last_z0
                do i=1, Npts
                    Rcore_i = Rcore_min + (i-1)*Rcore_steps
                    Rparticle = Rcore_i + dshell
                    if (zshift < -Rparticle) then
                        hshell(i) = 0d0
                        hcore(i) = 0d0
                    else if (zshift < -Rcore_i) then
                        hshell(i) = packing_densities(iLayer)*&
                                (Rparticle**2 - zshift**2)/Rcore_i**2 * p_log(i)
                        hcore(i) = 0d0
                    else if (zshift  < Rcore_i) then
                        hshell(i) = packing_densities(iLayer)*&
                                    (2*Rcore_i*dshell + dshell**2)/Rcore_i**2 * p_log(i)
                        hcore(i) = packing_densities(iLayer)*&
                                    (Rcore_i**2 - zshift**2)/Rcore_i**2 * p_log(i)
                    else if (zshift < Rparticle) then 
                        hshell(i) = packing_densities(iLayer)*&
                                (Rparticle**2 - zshift**2)/Rcore_i**2 * p_log(i)
                        hcore(i) = 0d0
                    else
                        hshell(i) = 0d0
                        hcore(i) = 0d0
                    end if
                end do
                core_area(iz) = core_area(iz) + trapz_uniform(hcore, Rcore_steps, Npts)
                shell_area(iz) = shell_area(iz) + trapz_uniform(hshell, Rcore_steps, Npts)
            end do
        end do
        
        do iz=1, Nz
            result_sld(iz) = SLDbg + core_area(iz)*dSLDc + shell_area(iz)*dSLDs
        end do
    end subroutine N_layers_coated_spheres_part_dist_2

end module models
