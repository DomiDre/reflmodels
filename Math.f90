module math
implicit none
    double precision, parameter :: pi = acos(-1d0)
    complex*16, parameter :: ci = dcmplx(0.d0,1.d0)
    double precision, parameter :: sq2 = sqrt(2d0)
    double precision, parameter :: sq2pi = sqrt(2d0*pi)
    double precision, parameter :: two_sq3 = 2d0*sqrt(3d0)
    double precision, parameter :: deg = pi/180d0
    
    integer :: n_integration_cuts = 50
contains
    ! Help-Functions
    subroutine lognormal(x, mu, sig, lognormalval)
        double precision, intent(in) :: x, mu, sig
        double precision, intent(out) :: lognormalval
        lognormalval =&
            1d0 / (sq2pi*sig*x) * dexp(-0.5d0*(dlog(x/mu) / sig)**2)
    end subroutine lognormal

    subroutine gaussian(x, mu, sig, gaussianval)
        double precision, intent(in) :: x, mu, sig
        double precision, intent(out) :: gaussianval
        gaussianval =&
            1d0 / (sq2pi*sig) * dexp(-0.5d0*((x-mu) / sig )**2)
    end subroutine gaussian

    subroutine equipartition(x, mu, sig, equipartition_val)
        double precision, intent(in) :: x, mu, sig
        double precision, intent(out) :: equipartition_val
        equipartition_val =&
                    1d0/(sig*two_sq3)
    end subroutine equipartition

    subroutine one(x, mu, sig, oneval)
        double precision, intent(in) :: x, mu, sig
        double precision, intent(out) :: oneval
        oneval = 1d0
    end subroutine one
    
    subroutine integrate_prob_distribution(pmu, psig, pmin, pmax, prob_dist,&
                                           probability)
        double precision, intent(in) :: pmu, psig, pmin, pmax
        external prob_dist
        double precision, intent(out) :: probability

        double precision :: abserr
        integer :: neval, ier, last
        integer, dimension(50) :: iwork
        double precision, dimension(200) :: work
        
        interface
            subroutine prob_dist(x, mu, sig, y)
                double precision, intent(in) :: x, mu, sig
                double precision, intent(out) :: y
            end subroutine
        end interface
        
        if (psig > 0d0) then
            call dqag (probability_function, pmin, pmax, 0d0, 1d-10, 6,&
                    probability, abserr, neval, ier, n_integration_cuts,&
                    4*n_integration_cuts, last, iwork, work)
        else
            probability = 1d0
        end if
        
        contains 
            double precision function probability_function(x)
                double precision :: x
                call prob_dist(x, pmu, psig, probability_function)
            end function probability_function

    end subroutine integrate_prob_distribution
    
    subroutine integrate_size_distribution(qval, p, Np, &
                pint, pmin, pmax, psig, &
                ff_function, prob_dist, ff_intensity)
        double precision, intent(in) :: qval
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np, pint
        double precision, intent(in) :: pmin, pmax, psig
        external ff_function
        external prob_dist
        double precision, intent(out) :: ff_intensity
        
        double precision :: abserr
        integer :: neval, ier, last
        integer, dimension(50) :: iwork
        double precision, dimension(200) :: work
        
        interface
            double precision function ff_function(q, p, Np)
                double precision, intent(in) :: q
                double precision, dimension(Np), intent(in) :: p
                integer, intent(in) :: Np
            end function
            
            subroutine prob_dist(x, mu, sig, y)
                double precision, intent(in) :: x, mu, sig
                double precision, intent(out) :: y
            end subroutine
        end interface
        
        if (psig > 0d0) then
            call dqag (size_dist_function, pmin, pmax, 0d0, 1d-10, 6,&
                    ff_intensity, abserr, neval, ier, n_integration_cuts,&
                    4*n_integration_cuts, last, iwork, work)
        else
            ff_intensity = ff_function(qval, p, Np)
        end if
        
        contains
            double precision function size_dist_function(x)
                double precision :: x
                double precision, dimension(Np) :: hp
                double precision :: prob, formfactor
                hp = p
                hp(pint) = x
                formfactor = ff_function(qval, hp, Np)
                call prob_dist(x, p(pint), psig, prob)
                
                size_dist_function = formfactor * prob
            end function size_dist_function
    end subroutine integrate_size_distribution


    subroutine integrate_two_size_distributions(qval, p, Np, &
                pint1, pmin1, pmax1, psig1, &
                pint2, pmin2, pmax2, psig2, &
                ff_function, prob_dist1, prob_dist2, ff_intensity)
        double precision, intent(in) :: qval
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np, pint1, pint2
        double precision, intent(in) :: pmin1, pmax1, psig1
        double precision, intent(in) :: pmin2, pmax2, psig2
        external ff_function
        external prob_dist1
        external prob_dist2
        double precision, intent(out) :: ff_intensity
        
        double precision :: abserr1, abserr2
        integer :: neval1, ier1, last1, neval2, ier2, last2
        integer, dimension(50) :: iwork1, iwork2
        double precision, dimension(200) :: work1, work2
        double precision, dimension(Np) :: help_params
        
        interface
            double precision function ff_function(q, p, Np)
                double precision, intent(in) :: q
                double precision, dimension(Np), intent(in) :: p
                integer, intent(in) :: Np
            end function
            
            subroutine prob_dist1(x, mu, sig, y)
                double precision, intent(in) :: x, mu, sig
                double precision, intent(out) :: y
            end subroutine
            
            subroutine prob_dist2(x, mu, sig, y)
                double precision, intent(in) :: x, mu, sig
                double precision, intent(out) :: y
            end subroutine
        end interface
        
        help_params = p
        if (psig1 > 0 .AND. psig2 > 0) then
            call dqag(integrate12, pmin1, pmax1, 0d0, 1d-10, 6,&
                    ff_intensity, abserr1, neval1, ier1, n_integration_cuts,&
                    4*n_integration_cuts, last1, iwork1, work1)

        else if (psig1 > 0) then ! means psig2 <= 0
            call dqag(integrate01, pmin1, pmax1, 0d0, 1d-10, 6,&
                    ff_intensity, abserr1, neval1, ier1, n_integration_cuts,&
                    4*n_integration_cuts, last1, iwork1, work1)

        else if (psig2 > 0) then ! means psig1 <= 0
            call dqag(integrate02, pmin2, pmax2, 0d0, 1d-10, 6,&
                    ff_intensity, abserr2, neval2, ier2, n_integration_cuts,&
                    4*n_integration_cuts, last2, iwork2, work2)

        else ! means psig1 and psig2 <= 0
            ff_intensity = ff_function(qval, p, Np)
        end if
        
        contains
            double precision function integrate12(x)
                double precision :: x
                double precision :: function_value, prob
                help_params(pint1) = x

                call dqag (integrate02, pmin2, pmax2, 0d0, 1d-10, 6,&
                    function_value, abserr2, neval2, ier2, n_integration_cuts,&
                    4*n_integration_cuts, last2, iwork2, work2)
                call prob_dist1(x, p(pint1), psig1, prob)

                integrate12 = function_value * prob
            end function integrate12
            
            double precision function integrate01(x)
                double precision :: x
                double precision :: function_value, prob
                help_params(pint1) = x

                function_value = ff_function(qval, help_params, Np)
                call prob_dist1(x, p(pint1), psig1, prob)

                integrate01 = function_value * prob
            end function integrate01
            
            double precision function integrate02(x)
                double precision :: x
                double precision :: function_value, prob
                help_params(pint2) = x

                function_value = ff_function(qval, help_params, Np)
                call prob_dist2(x, p(pint2), psig2, prob)

                integrate02 = function_value * prob
            end function integrate02
    end subroutine integrate_two_size_distributions

    subroutine integrate_three_size_distributions(qval, p, Np, &
                pint1, pmin1, pmax1, psig1, &
                pint2, pmin2, pmax2, psig2, &
                pint3, pmin3, pmax3, psig3, &
                ff_function, prob_dist1, prob_dist2, prob_dist3, ff_intensity)
        double precision, intent(in) :: qval
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np, pint1, pint2, pint3
        double precision, intent(in) :: pmin1, pmax1, psig1
        double precision, intent(in) :: pmin2, pmax2, psig2
        double precision, intent(in) :: pmin3, pmax3, psig3
        external ff_function
        external prob_dist1
        external prob_dist2
        external prob_dist3
        double precision, intent(out) :: ff_intensity
        
        double precision :: abserr1, abserr2, abserr3
        integer :: neval1, ier1, last1, neval2, ier2, last2, neval3, ier3, last3
        integer, dimension(50) :: iwork1, iwork2, iwork3
        double precision, dimension(200) :: work1, work2, work3
        double precision, dimension(Np) :: help_params
        
        interface
            double precision function ff_function(q, p, Np)
                double precision, intent(in) :: q
                double precision, dimension(Np), intent(in) :: p
                integer, intent(in) :: Np
            end function
            
            subroutine prob_dist1(x, mu, sig, y)
                double precision, intent(in) :: x, mu, sig
                double precision, intent(out) :: y
            end subroutine
            
            subroutine prob_dist2(x, mu, sig, y)
                double precision, intent(in) :: x, mu, sig
                double precision, intent(out) :: y
            end subroutine
            
            subroutine prob_dist3(x, mu, sig, y)
                double precision, intent(in) :: x, mu, sig
                double precision, intent(out) :: y
            end subroutine
        end interface
        
        help_params = p
        if (psig1 > 0 .AND. psig2 > 0 .AND. psig3 > 0) then
            call dqag(integrate123, pmin1, pmax1, 0d0, 1d-10, 6,&
                    ff_intensity, abserr1, neval1, ier1, n_integration_cuts,&
                    4*n_integration_cuts, last1, iwork1, work1)
        else if (psig1 > 0 .AND. psig2 > 0) then ! means psig3 <= 0
            call dqag(integrate012, pmin1, pmax1, 0d0, 1d-10, 6,&
                    ff_intensity, abserr1, neval1, ier1, n_integration_cuts,&
                    4*n_integration_cuts, last1, iwork1, work1)
        else if (psig2 > 0 .AND. psig3 > 0) then ! means psig1 <= 0
            call dqag(integrate023, pmin2, pmax2, 0d0, 1d-10, 6,&
                    ff_intensity, abserr2, neval2, ier2, n_integration_cuts,&
                    4*n_integration_cuts, last2, iwork2, work2)
        else if (psig1 > 0 .AND. psig3 > 0) then ! means psig2 <= 0
            call dqag(integrate013, pmin1, pmax1, 0d0, 1d-10, 6,&
                    ff_intensity, abserr1, neval1, ier1, n_integration_cuts,&
                    4*n_integration_cuts, last1, iwork1, work1)
        else if (psig1 > 0) then ! means psig2 <= 0 && psig3 <= 0
            call dqag(integrate001, pmin1, pmax1, 0d0, 1d-10, 6,&
                    ff_intensity, abserr1, neval1, ier1, n_integration_cuts,&
                    4*n_integration_cuts, last1, iwork1, work1)
        else if (psig2 > 0) then ! means psig1 <= 0 && psig3 <= 0
            call dqag(integrate002, pmin2, pmax2, 0d0, 1d-10, 6,&
                    ff_intensity, abserr2, neval2, ier2, n_integration_cuts,&
                    4*n_integration_cuts, last2, iwork2, work2)
        else if (psig3 > 0) then ! means psig1 <= 0 && psig2 <= 0
            call dqag(integrate003, pmin3, pmax3, 0d0, 1d-10, 6,&
                    ff_intensity, abserr3, neval3, ier3, n_integration_cuts,&
                    4*n_integration_cuts, last3, iwork3, work3)
        else ! means psig1 && psig2 && psig3 <= 0
            ff_intensity = ff_function(qval, p, Np)
        end if
        
        contains
            double precision function integrate123(x)
                double precision :: x
                double precision :: function_value, prob
                help_params(pint1) = x

                call dqag (integrate023, pmin2, pmax2, 0d0, 1d-10, 6,&
                    function_value, abserr2, neval2, ier2, n_integration_cuts,&
                    4*n_integration_cuts, last2, iwork2, work2)
                call prob_dist1(x, p(pint1), psig1, prob)

                integrate123 = function_value * prob
            end function integrate123
            
            double precision function integrate023(x)
                double precision :: x
                double precision :: function_value, prob
                help_params(pint2) = x

                call dqag (integrate003, pmin3, pmax3, 0d0, 1d-10, 6,&
                    function_value, abserr3, neval3, ier3, n_integration_cuts,&
                    4*n_integration_cuts, last3, iwork3, work3)
                call prob_dist2(x, p(pint2), psig2, prob)

                integrate023 = function_value * prob
            end function integrate023
            
            double precision function integrate013(x)
                double precision :: x
                double precision :: function_value, prob
                help_params(pint1) = x

                call dqag (integrate003, pmin3, pmax3, 0d0, 1d-10, 6,&
                    function_value, abserr3, neval3, ier3, n_integration_cuts,&
                    4*n_integration_cuts, last3, iwork3, work3)
                call prob_dist1(x, p(pint1), psig1, prob)

                integrate013 = function_value * prob
            end function integrate013
            
            double precision function integrate012(x)
                double precision :: x
                double precision :: function_value, prob
                help_params(pint1) = x

                call dqag (integrate002, pmin2, pmax2, 0d0, 1d-10, 6,&
                    function_value, abserr2, neval2, ier2, n_integration_cuts,&
                    4*n_integration_cuts, last2, iwork2, work2)
                call prob_dist1(x, p(pint1), psig1, prob)

                integrate012 = function_value * prob
            end function integrate012
            
            double precision function integrate001(x)
                double precision :: x
                double precision :: function_value, prob
                help_params(pint1) = x

                function_value = ff_function(qval, help_params, Np)
                call prob_dist1(x, p(pint1), psig1, prob)

                integrate001 = function_value * prob
            end function integrate001
            
            double precision function integrate002(x)
                double precision :: x
                double precision :: function_value, prob
                help_params(pint2) = x

                function_value = ff_function(qval, help_params, Np)
                call prob_dist2(x, p(pint2), psig2, prob)

                integrate002 = function_value * prob
            end function integrate002
            
            double precision function integrate003(x)
                double precision :: x
                double precision :: function_value, prob
                help_params(pint3) = x

                function_value = ff_function(qval, help_params, Np)
                call prob_dist3(x, p(pint3), psig3, prob)

                integrate003 = function_value * prob
            end function integrate003
    end subroutine integrate_three_size_distributions

    subroutine twodim_integral_variable_bounds(xmin, xmax, &
                                               ymin_func, ymax_func, &
                                               p, Np, &
                                               integrand_func, integrand_res)
        double precision, intent(in) :: xmin, xmax
        external ymin_func
        external ymax_func
        double precision, dimension(Np), intent(in) :: p
        integer, intent(in) :: Np
        external integrand_func
        double precision, intent(out) :: integrand_res
        
        double precision :: abserr1, abserr2
        integer :: neval1, ier1, last1, neval2, ier2, last2
        integer, dimension(50) :: iwork1, iwork2
        double precision, dimension(200) :: work1, work2
        double precision, dimension(Np+1) :: help_params
        
        interface
            double precision function ymin_func(x, p, Np)
                double precision, intent(in) :: x
                double precision, dimension(Np), intent(in) :: p
                integer, intent(in) :: Np
            end function
            
            double precision function ymax_func(x, p, Np)
                double precision, intent(in) :: x
                double precision, dimension(Np), intent(in) :: p
                integer, intent(in) :: Np
            end function
            
            double precision function integrand_func(x, y, p, Np)
                double precision, intent(in) :: x, y
                double precision, dimension(Np), intent(in) :: p
                integer, intent(in) :: Np
            end function
        end interface
        
    
        help_params(1:Np) = p
        help_params(Np+1) = 0d0
        call dqag(integrate12, xmin, xmax, 0d0, 1d-10, 6,&
                integrand_res, abserr1, neval1, ier1, n_integration_cuts,&
                4*n_integration_cuts, last1,&
                iwork1, work1)
        
        contains
            double precision function integrate12(x)
                double precision :: x
                double precision :: ymin, ymax
                ymin = ymin_func(x, p, Np)
                ymax = ymax_func(x, p, Np)
                help_params(Np+1) = x
                call dqag (integrate02, ymin, ymax, 0d0, 1d-10, 6,&
                    integrate12, abserr2, neval2, ier2, n_integration_cuts,&
                    4*n_integration_cuts, last2,&
                    iwork2, work2)
            end function integrate12
            
            double precision function integrate02(y)
                double precision :: y
                double precision :: x
                
                x = help_params(Np+1)
                integrate02 = integrand_func(x, y, p, Np)
            end function integrate02
    end subroutine twodim_integral_variable_bounds
    
    subroutine get_cutoff_gaussian(mu, sig, left_cutoff, right_cutoff)
        double precision, intent(in) :: mu, sig
        double precision, intent(out) :: left_cutoff, right_cutoff
        
        if (sig > 0d0) then
            left_cutoff = mu - 5d0*sig
            right_cutoff = mu + 5d0*sig
        else
            left_cutoff = 0d0
            right_cutoff = 0d0
        end if
    end subroutine get_cutoff_gaussian

    subroutine get_cutoff_lognormal(mu, sig, left_cutoff, right_cutoff)
        double precision, intent(in) :: mu, sig
        double precision, intent(out) :: left_cutoff, right_cutoff
        
        double precision :: prob_integrated
        
        double precision :: abserr
        integer :: neval, ier, last
        integer, dimension(50) :: iwork
        double precision, dimension(200) :: work
        
        if (sig > 0d0) then
            left_cutoff = mu*exp(-5*sig)
            right_cutoff = mu*exp(-sig**2) + sig*mu
            prob_integrated = 0d0
            ! Shift right cutoff of integral, until lognormal distribution
            ! integrates to 1 with 1d-6 tolerance.
            do while (abs(prob_integrated-1d0) > 1d-6)
                call dqag(lognormal_function, left_cutoff, right_cutoff, &
                          0d0, 1d-10, 6, prob_integrated, abserr, neval, &
                          ier, n_integration_cuts, 4*n_integration_cuts, last, iwork, work)
                right_cutoff = right_cutoff + sig*mu
            end do
        else
            left_cutoff = 0d0
            right_cutoff = 0d0
        end if
        contains
            double precision function lognormal_function(x)
                double precision :: x
                call lognormal(x, mu, sig, lognormal_function)
            end function lognormal_function
    end subroutine get_cutoff_lognormal
    
    subroutine mean(x, Nx, meanx)
        double precision, dimension(Nx), intent(in) :: x
        integer, intent(in) :: Nx
        double precision, intent(out) :: meanx
        
        meanx = sum(x) / Nx
    end subroutine mean
    
    subroutine simple_linear_fit(x, y, Nx, slope, intercept)
        double precision, dimension(Nx), intent(in) :: x, y
        integer, intent(in) :: Nx
        double precision, intent(out) :: slope, intercept
        
        double precision :: mean_xy, mean_x, mean_y, mean_x2
        
        call mean(x, Nx, mean_x)
        call mean(y, Nx, mean_y)
        call mean(x*y, Nx, mean_xy)
        call mean(x*x, Nx, mean_x2)
        
        slope = (mean_xy-mean_x*mean_y) / (mean_x2-mean_x**2)
        intercept = mean_y - slope*mean_x
    end subroutine simple_linear_fit
    
    subroutine get_interpolated_value(x, y, xval, Nx, yval)
        double precision, dimension(Nx), intent(in) :: x, y
        double precision, intent(in) :: xval
        integer, intent(in) :: Nx
        double precision, intent(out) :: yval
        
        integer :: ix
        double precision :: x_left, x_right, y_left, y_right
        
        if (xval <= x(1)) then
            yval = (y(2)-y(1))*(xval-x(1))/(x(2)-x(1)) + y(1)
            return
        end if
        
        if (xval >= x(Nx)) then
            yval = (y(Nx)-y(Nx-1))*(xval-x(Nx))/(x(Nx)-x(Nx-1)) + y(Nx)
            return
        end if
        
        do ix=1, Nx
            if (x(ix) > xval) then
                x_left = x(ix-1)
                x_right = x(ix)
                y_left = y(ix-1)
                y_right = y(ix)
                exit
            end if
        end do
        yval = (y_right-y_left)*(xval-x_right)/(x_right-x_left) + y_right
    end subroutine get_interpolated_value
        
    subroutine resolution_smear(q, I, sigQ, Nq, Ismear)
        double precision, dimension(Nq), intent(in) :: q, I, sigQ
        integer, intent(in) :: Nq
        double precision, dimension(Nq), intent(out) :: Ismear
        
        integer :: iq, jq, iq_tilde, iq_low, iq_high
        double precision :: cq, csigq, cqtilde, prob_qtilde, Itilde
        double precision :: cqtilde_before, cqtilde_next
        double precision :: sum_prob_qtilde, sum_Itilde
        
        double precision :: q_spacing_low, slope_low, intercept_low
        double precision :: q_spacing_high, slope_high, intercept_high
        double precision :: cqhelp
        if (Nq .eq. 1) then
            Ismear = I
            return
        end if
        
        q_spacing_low = q(2) - q(1)
        slope_low = (I(2) - I(1)) / q_spacing_low
        intercept_low = I(1)-slope_low*q(1)

        q_spacing_high = q(Nq) - q(Nq-1)
        slope_high = (I(Nq) - I(Nq-1)) / q_spacing_high
        intercept_high = I(Nq) - slope_high*q(Nq)
        
        do iq=1, Nq
            csigq = sigQ(iq)
            if (csigq <= 0d0) then
                Ismear(iq) = I(iq)
                cycle
            end if
            
            cq = q(iq)
            iq_low = iq
            cqhelp = cq
            do while((cqhelp > cq-3*csigq))
                iq_low = iq_low - 1
                if (iq_low >= 1) then
                    cqhelp = q(iq_low)
                else
                    cqhelp = q(1) + (iq_low-1)*q_spacing_low
                end if
            end do

            iq_high = iq
            cqhelp = cq
            do while((cqhelp < cq+3*csigq))
                iq_high = iq_high + 1
                if (iq_high <= Nq) then
                    cqhelp = q(iq_high)
                else
                    cqhelp = q(Nq) + (iq_high-Nq)*q_spacing_high
                end if
            end do
            
            ! Do integration from I(q_low) to I(q_high) with gaussian probability
            ! Extrapolate function linearly outside of accessible region
            if(iq_low < 1) then
                cqtilde = q(1) + (iq_low-1) * q_spacing_low
                cqtilde_next = q(1) + iq_low * q_spacing_low
                Itilde = slope_low*cqtilde + intercept_low
            else
                cqtilde = q(iq_low)
                cqtilde_next = q(iq_low+1)
                Itilde = I(iq_low)
            end if

            call gaussian(cq, cqtilde, csigq, prob_qtilde)
            sum_prob_qtilde = prob_qtilde*(cqtilde_next-cqtilde)
            sum_Itilde = prob_qtilde*Itilde*(cqtilde_next-cqtilde)
            
            do jq=iq_low+1, iq_high-1
                cqtilde_before = cqtilde
                cqtilde = cqtilde_next
                
                if(jq <= 0) then
                    cqtilde_next = q(1) + jq * q_spacing_low
                    Itilde = slope_low*cqtilde + intercept_low
                else if (jq >= Nq) then
                    cqtilde_next = q(Nq) + (jq-Nq+1) * q_spacing_high
                    Itilde = slope_high*cqtilde + intercept_high
                else
                    cqtilde_next = q(jq+1)
                    Itilde = I(jq)
                end if
                
                if (Itilde < 0d0) cycle
                call gaussian(cq, cqtilde, csigq, prob_qtilde)
                sum_prob_qtilde = sum_prob_qtilde + prob_qtilde*&
                                                (cqtilde_next-cqtilde_before)
                sum_Itilde = sum_Itilde + prob_qtilde*Itilde*&
                                                (cqtilde_next-cqtilde_before)
            end do

            cqtilde_before = cqtilde
            cqtilde = cqtilde_next
            
            if(iq_high >= Nq) then
                Itilde = slope_high*cqtilde+intercept_high
            else
                Itilde = I(iq_high)
            end if
            call gaussian(cq, cqtilde, csigq, prob_qtilde)
            sum_prob_qtilde = sum_prob_qtilde + prob_qtilde*&
                                            (cqtilde_next-cqtilde)
            sum_Itilde = sum_Itilde + prob_qtilde*Itilde*&
                                            (cqtilde_next-cqtilde)

            Ismear(iq) = sum_Itilde / sum_prob_qtilde
        end do
    end subroutine resolution_smear
    
    subroutine resolution_smear_interpolating(q, I, sigQ, Nq, Ismear)
        double precision, dimension(Nq), intent(in) :: q, I, sigQ
        integer, intent(in) :: Nq
        double precision, dimension(Nq), intent(out) :: Ismear
        
        integer :: iq, jq, iq_tilde, iq_low, iq_high
        double precision :: cq, csigq, cqtilde, prob_qtilde, Itilde
        double precision :: cqtilde_before, cqtilde_next
        double precision :: sum_prob_qtilde, sum_Itilde
        
        double precision :: q_spacing
        double precision :: cq_running, q_low_value, q_high_value
        double precision :: qlow_steps, qhigh_steps
        if (Nq .eq. 1) then
            Ismear = I
            return
        end if
        
        call mean(q(2:Nq) - q(1:Nq-1), Nq-1, q_spacing)
        
        do iq=1, Nq
            csigq = sigQ(iq)
            if (csigq <= 0d0) then
                Ismear(iq) = I(iq)
                cycle
            end if
            cq = q(iq)
            cq_running = cq
            iq_low = iq
            
            do while((cq_running > cq-3*csigq))
                iq_low = iq_low - 1
                cq_running = cq - (iq-iq_low)*q_spacing
            end do

            iq_high = iq
            cq_running = cq
            do while((cq_running < cq+3*csigq))
                iq_high = iq_high + 1
                cq_running = cq + (iq_high-iq)*q_spacing
            end do
            
            ! Do integration from I(q_low) to I(q_high) with gaussian probability
            ! Interpolate values from given arrays
            ! Extrapolate function linearly outside of accessible region
            cqtilde = cq - (iq-iq_low) * q_spacing
            cqtilde_next = cqtilde + q_spacing
            call get_interpolated_value(q, I, cqtilde, Nq, Itilde)
            call gaussian(cq, cqtilde, csigq, prob_qtilde)
            sum_prob_qtilde = prob_qtilde*(cqtilde_next-cqtilde)
            sum_Itilde = prob_qtilde*Itilde*(cqtilde_next-cqtilde)
            
            do jq=iq_low+1, iq_high-1
                cqtilde_before = cqtilde
                cqtilde = cqtilde_next
                cqtilde_next = cqtilde + q_spacing
                call get_interpolated_value(q, I, cqtilde, Nq, Itilde)
                call gaussian(cq, cqtilde, csigq, prob_qtilde)
                sum_prob_qtilde = sum_prob_qtilde + prob_qtilde*&
                                                (cqtilde_next-cqtilde_before)
                sum_Itilde = sum_Itilde + prob_qtilde*Itilde*&
                                                (cqtilde_next-cqtilde_before)
            end do

            cqtilde_before = cqtilde
            cqtilde = cqtilde_next
            call get_interpolated_value(q, I, cqtilde, Nq, Itilde)
            call gaussian(cq, cqtilde, csigq, prob_qtilde)
            sum_prob_qtilde = sum_prob_qtilde + prob_qtilde*&
                                            (cqtilde_next-cqtilde)
            sum_Itilde = sum_Itilde + prob_qtilde*Itilde*&
                                            (cqtilde_next-cqtilde)
            Ismear(iq) = sum_Itilde / sum_prob_qtilde
        end do
    end subroutine resolution_smear_interpolating
    
    subroutine angle_correct(q, thoffset, wavelength, Nq, qcorr)
        integer, intent(in):: Nq
        double precision, intent(in), dimension(Nq) :: q
        double precision, intent(in) :: thoffset, wavelength
        double precision, intent(out), dimension(Nq) :: qcorr

        integer :: iq
        do iq=1, Nq
            qcorr(iq) = 4d0*pi/wavelength * sin(asin(wavelength*q(iq)/4d0/pi) - &
                    thoffset*deg)
        end do
    end subroutine angle_correct    

    subroutine footprint_correct(q, reflin, samplelen, beamwidth, &
            wavelength, Nq, reflcorr)
        double precision, intent(in), dimension(Nq) :: q, reflin
        double precision, intent(in) :: samplelen, beamwidth, wavelength
        integer, intent(in) :: Nq
        double precision, intent(out), dimension(Nq) :: reflcorr
        
        double precision :: sinth, sinthb
        integer :: iq
        sinthb = beamwidth/samplelen ! assume uniform beam profile
        do iq=1, Nq
            sinth = q(iq)*wavelength/(4d0*pi)
            if (sinth < sinthb) then
                reflcorr(iq) = reflin(iq)*sinth/sinthb
            else
                reflcorr(iq) = reflin(iq) !*erf(samplelen/beamwidth*sinth/2d0/sq2)
            end if
        end do
    end subroutine footprint_correct

    double precision function trapz_uniform(y, dx, Ny)
        implicit none
        integer, intent(in) :: Ny
        double precision, intent(in), dimension(Ny) :: y
        double precision, intent(in) :: dx
        integer :: i
        trapz_uniform = 0.5d0*(y(1) + y(Ny))
        do i=2, Ny-1
            trapz_uniform = trapz_uniform + y(i)
        end do
        trapz_uniform = dx * trapz_uniform
    end function trapz_uniform
    
end module math
