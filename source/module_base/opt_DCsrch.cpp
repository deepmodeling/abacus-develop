#include "opt_DCsrch.h"

#include <math.h>
#include <string.h>

// This file is translated from fortran codes dcstep.f of scipy.
// The structure and all annotation of the original file have been retained.
// See original source at https://github.com/scipy/scipy/blob/main/scipy/optimize/minpack2/dcstep.f.
// sunliang 2022-05-30

namespace ModuleBase
{
int dcsrch(double& stp,
           double& f,
           double& g,
           double& ftol,
           double& gtol,
           double& xtol,
           char* task,
           double& stpmin,
           double& stpmax,
           int* isave,
           double* dsave)
{
    // c     **********
    // c
    // c     Subroutine dcsrch
    // c
    // c     This subroutine finds a step that satisfies a sufficient
    // c     decrease condition and a curvature condition.
    // c
    // c     Each call of the subroutine updates an interval with
    // c     endpoints stx and sty. The interval is initially chosen
    // c     so that it contains a minimizer of the modified function
    // c
    // c           psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).
    // c
    // c     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
    // c     interval is chosen so that it contains a minimizer of f.
    // c
    // c     The algorithm is designed to find a step that satisfies
    // c     the sufficient decrease condition
    // c
    // c           f(stp) <= f(0) + ftol*stp*f'(0),
    // c
    // c     and the curvature condition
    // c
    // c           abs(f'(stp)) <= gtol*abs(f'(0)).
    // c
    // c     If ftol is less than gtol and if, for example, the function
    // c     is bounded below, then there is always a step which satisfies
    // c     both conditions.
    // c
    // c     If no step can be found that satisfies both conditions, then
    // c     the algorithm stops with a warning. In this case stp only
    // c     satisfies the sufficient decrease condition.
    // c
    // c     A typical invocation of dcsrch has the following outline:
    // c
    // c     Evaluate the function at stp = 0.0d0; store in f.
    // c     Evaluate the gradient at stp = 0.0d0; store in g.
    // c     Choose a starting step stp.
    // c
    // c     task = 'START'
    // c  10 continue
    // c        call dcsrch(stp,f,g,ftol,gtol,xtol,task,stpmin,stpmax,
    // c    +               isave,dsave)
    // c        if (task .eq. 'FG') then
    // c           Evaluate the function and the gradient at stp
    // c           go to 10
    // c           end if
    // c
    // c     NOTE: The user must not alter work arrays between calls.
    // c
    // c     The subroutine statement is
    // c
    // c       subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
    // c                         task,isave,dsave)
    // c     where
    // c
    // c       stp is a double precision variable.
    // c         On entry stp is the current estimate of a satisfactory
    // c            step. On initial entry, a positive initial estimate
    // c            must be provided.
    // c         On exit stp is the current estimate of a satisfactory step
    // c            if task = 'FG'. If task = 'CONV' then stp satisfies
    // c            the sufficient decrease and curvature condition.
    // c
    // c       f is a double precision variable.
    // c         On initial entry f is the value of the function at 0.
    // c            On subsequent entries f is the value of the
    // c            function at stp.
    // c         On exit f is the value of the function at stp.
    // c
    // c       g is a double precision variable.
    // c         On initial entry g is the derivative of the function at 0.
    // c            On subsequent entries g is the derivative of the
    // c            function at stp.
    // c         On exit g is the derivative of the function at stp.
    // c
    // c       ftol is a double precision variable.
    // c         On entry ftol specifies a nonnegative tolerance for the
    // c            sufficient decrease condition.
    // c         On exit ftol is unchanged.
    // c
    // c       gtol is a double precision variable.
    // c         On entry gtol specifies a nonnegative tolerance for the
    // c            curvature condition.
    // c         On exit gtol is unchanged.
    // c
    // c       xtol is a double precision variable.
    // c         On entry xtol specifies a nonnegative relative tolerance
    // c            for an acceptable step. The subroutine exits with a
    // c            warning if the relative difference between sty and stx
    // c            is less than xtol.
    // c         On exit xtol is unchanged.
    // c
    // c       task is a character variable of length at least 60.
    // c         On initial entry task must be set to 'START'.
    // c         On exit task indicates the required action:
    // c
    // c            If task(1:2) = 'FG' then evaluate the function and
    // c            derivative at stp and call dcsrch again.
    // c
    // c            If task(1:4) = 'CONV' then the search is successful.
    // c
    // c            If task(1:4) = 'WARN' then the subroutine is not able
    // c            to satisfy the convergence conditions. The exit value of
    // c            stp contains the best point found during the search.
    // c
    // c            If task(1:5) = 'ERROR' then there is an error in the
    // c            input arguments.
    // c
    // c         On exit with convergence, a warning or an error, the
    // c            variable task contains additional information.
    // c
    // c       stpmin is a double precision variable.
    // c         On entry stpmin is a nonnegative lower bound for the step.
    // c         On exit stpmin is unchanged.
    // c
    // c       stpmax is a double precision variable.
    // c         On entry stpmax is a nonnegative upper bound for the step.
    // c         On exit stpmax is unchanged.
    // c
    // c       isave is an integer work array of dimension 2.
    // c
    // c       dsave is a double precision work array of dimension 13.
    // c
    // c     Subprograms called
    // c
    // c       MINPACK-2 ... dcstep
    // c
    // c     MINPACK-1 Project. June 1983.
    // c     Argonne National Laboratory.
    // c     Jorge J. More' and David J. Thuente.
    // c
    // c     MINPACK-2 Project. November 1993.
    // c     Argonne National Laboratory and University of Minnesota.
    // c     Brett M. Averick, Richard G. Carter, and Jorge J. More'.
    // c
    // c     **********
    double zero = 0.0;
    double p5 = 0.5;
    double p66 = 0.66;
    double xtrapl = 1.1;
    double xtrapu = 4.0;

    bool brackt;
    int stage;
    double finit, ftest, fm, fx, fxm, fy, fym, ginit, gtest, gm, gx, gxm, gy, gym, stx, sty, stmin, stmax, width,
        width1;

    extern /* Subroutine */ void dcstep(double&,
                                        double&,
                                        double&,
                                        double&,
                                        double&,
                                        double&,
                                        double&,
                                        double&,
                                        double&,
                                        bool&,
                                        double&,
                                        double&);
    // c     Initialization block.
    if (strncmp(task, "START", 5) == 0)
    {
        // c        Check the input arguments for errors.
        if (stp < stpmin)
        {
            strcpy(task, "ERROR: STP .LT. STPMIN");
        }
        if (stp > stpmax)
        {
            strcpy(task, "ERROR: STP .GT. STPMAX");
        }
        if (g >= 0.)
        {
            strcpy(task, "ERROR: INITIAL G .GE. ZERO");
        }
        if (ftol < 0.)
        {
            strcpy(task, "ERROR: FTOL .LT. ZERO");
        }
        if (gtol < 0.)
        {
            strcpy(task, "ERROR: GTOL .LT. ZERO");
        }
        if (xtol < 0.)
        {
            strcpy(task, "ERROR: XTOL .LT. ZERO");
        }
        if (stpmin < 0.)
        {
            strcpy(task, "ERROR: STPMIN .LT. ZERO");
        }
        if (stpmax < stpmin)
        {
            strcpy(task, "ERROR: STPMAX .LT. STPMIN");
        }

        // c        Exit if there are errors on input.

        if (strncmp(task, "ERROR", 5) == 0)
        {
            return 0;
        }
        // c        Initialize local variables.
        brackt = false;
        stage = 1;
        finit = f;
        ginit = g;
        gtest = ftol * ginit;
        width = stpmax - stpmin;
        width1 = width / p5;

        // c        The variables stx, fx, gx contain the values of the step,
        // c        function, and derivative at the best step.
        // c        The variables sty, fy, gy contain the value of the step,
        // c        function, and derivative at sty.
        // c        The variables stp, f, g contain the values of the step,
        // c        function, and derivative at stp.

        stx = zero;
        fx = finit;
        gx = ginit;
        sty = zero;
        fy = finit;
        gy = ginit;
        stmin = zero;
        stmax = stp + stp * xtrapu;
        strcpy(task, "FG");
        goto L10;
    }
    else
    {

        // c        Restore local variables.

        if (isave[1] == 1)
        {
            brackt = true;
        }
        else
        {
            brackt = false;
        }
        stage = isave[2];
        ginit = dsave[1];
        gtest = dsave[2];
        gx = dsave[3];
        gy = dsave[4];
        finit = dsave[5];
        fx = dsave[6];
        fy = dsave[7];
        stx = dsave[8];
        sty = dsave[9];
        stmin = dsave[10];
        stmax = dsave[11];
        width = dsave[12];
        width1 = dsave[13];
    }

    // c     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
    // c     algorithm enters the second stage.

    ftest = finit + stp * gtest;
    if (stage == 1 && f <= ftest && g >= zero)
        stage = 2;

    // c     Test for warnings.

    if (brackt && (stp <= stmin || stp >= stmax))
    {
        strcpy(task, "WARNING: ROUNDING ERRORS PREVENT PROGRESS");
    }
    if (brackt && stmax - stmin <= xtol * stmax)
    {
        strcpy(task, "WARNING: XTOL TEST SATISFIED");
    }
    if (stp == stpmax && f <= ftest && g <= gtest)
    {
        strcpy(task, "WARNING: STP = STPMAX");
    }
    if (stp == stpmin && (f > ftest || g >= gtest))
    {
        strcpy(task, "WARNING: STP = STPMIN");
    }

    // c     Test for convergence.

    if (f <= ftest && std::abs(g) <= gtol * (-ginit))
    {
        strcpy(task, "CONVERGENCE");
        // strcpy(task, "CONVERGENCE", 11);
    }

    // c     Test for termination.

    if (strncmp(task, "WARN", 4) == 0 || strncmp(task, "CONV", 4) == 0)
    {
        goto L10;
    }

    // c     A modified function is used to predict the step during the
    // c     first stage if a lower function value has been obtained but
    // c     the decrease is not sufficient.

    if (stage == 1 && f <= fx && f > ftest)
    {

        // c        Define the modified function and derivative values.

        fm = f - stp * gtest;
        fxm = fx - stx * gtest;
        fym = fy - sty * gtest;
        gm = g - gtest;
        gxm = gx - gtest;
        gym = gy - gtest;

        // c        Call dcstep to update stx, sty, and to compute the new step.

        dcstep(stx, fxm, gxm, sty, fym, gym, stp, fm, gm, brackt, stmin, stmax);

        // c        Reset the function and derivative values for f.

        fx = fxm + stx * gtest;
        fy = fym + sty * gtest;
        gx = gxm + gtest;
        gy = gym + gtest;
    }
    else
    {

        // c       Call dcstep to update stx, sty, and to compute the new step.

        dcstep(stx, fx, gx, sty, fy, gy, stp, f, g, brackt, stmin, stmax);
    }
    // c     Decide if a bisection step is needed.
    if (brackt)
    {
        if (std::abs(sty - stx) >= p66 * width1)
            stp = stx + p5 * (sty - stx);
        width1 = width;
        width = std::abs(sty - stx);
    }
    // c     Set the minimum and maximum steps allowed for stp.

    if (brackt)
    {
        stmin = std::min(stx, sty);
        stmax = std::max(stx, sty);
    }
    else
    {
        stmin = stp + xtrapl * (stp - stx);
        stmax = stp + xtrapu * (stp - stx);
    }

    // c     Force the step to be within the bounds stpmax and stpmin.
    stp = std::max(stp, stpmin);
    stp = std::min(stp, stpmax);
    // c     If further progress is not possible, let stp be the best
    // c     point obtained during the search.
    if ((brackt && (stp <= stmin || stp >= stmax)) || (brackt && stmax - stmin <= xtol * stmax))
    {
        stp = stx;
    }
    // c     Obtain another function and derivative.

    strcpy(task, "FG");
L10:
    // c     Save local variables.
    if (brackt)
    {
        isave[1] = 1;
    }
    else
    {
        isave[1] = 0;
    }
    isave[2] = stage;
    dsave[1] = ginit;
    dsave[2] = gtest;
    dsave[3] = gx;
    dsave[4] = gy;
    dsave[5] = finit;
    dsave[6] = fx;
    dsave[7] = fy;
    dsave[8] = stx;
    dsave[9] = sty;
    dsave[10] = stmin;
    dsave[11] = stmax;
    dsave[12] = width;
    dsave[13] = width1;
    return 0;
}

/* Subroutine */ void dcstep(double& stx,
                             double& fx,
                             double& dx,
                             double& sty,
                             double& fy,
                             double& dy,
                             double& stp,
                             double& fp,
                             double& dp,
                             bool& brackt,
                             double& stpmin,
                             double& stpmax)
{
    // c     **********
    // c
    // c     Subroutine dcstep
    // c
    // c     This subroutine computes a safeguarded step for a search
    // c     procedure and updates an interval that contains a step that
    // c     satisfies a sufficient decrease and a curvature condition.
    // c
    // c     The parameter stx contains the step with the least function
    // c     value. If brackt is set to .true. then a minimizer has
    // c     been bracketed in an interval with endpoints stx and sty.
    // c     The parameter stp contains the current step.
    // c     The subroutine assumes that if brackt is set to .true. then
    // c
    // c           min(stx,sty) < stp < max(stx,sty),
    // c
    // c     and that the derivative at stx is negative in the direction
    // c     of the step.
    // c
    // c     The subroutine statement is
    // c
    // c       subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
    // c                         stpmin,stpmax)
    // c
    // c     where
    // c
    // c       stx is a double precision variable.
    // c         On entry stx is the best step obtained so far and is an
    // c            endpoint of the interval that contains the minimizer.
    // c         On exit stx is the updated best step.
    // c
    // c       fx is a double precision variable.
    // c         On entry fx is the function at stx.
    // c         On exit fx is the function at stx.
    // c
    // c       dx is a double precision variable.
    // c         On entry dx is the derivative of the function at
    // c            stx. The derivative must be negative in the direction of
    // c            the step, that is, dx and stp - stx must have opposite
    // c            signs.
    // c         On exit dx is the derivative of the function at stx.
    // c
    // c       sty is a double precision variable.
    // c         On entry sty is the second endpoint of the interval that
    // c            contains the minimizer.
    // c         On exit sty is the updated endpoint of the interval that
    // c            contains the minimizer.
    // c
    // c       fy is a double precision variable.
    // c         On entry fy is the function at sty.
    // c         On exit fy is the function at sty.
    // c
    // c       dy is a double precision variable.
    // c         On entry dy is the derivative of the function at sty.
    // c         On exit dy is the derivative of the function at the exit sty.
    // c
    // c       stp is a double precision variable.
    // c         On entry stp is the current step. If brackt is set to .true.
    // c            then on input stp must be between stx and sty.
    // c         On exit stp is a new trial step.
    // c
    // c       fp is a double precision variable.
    // c         On entry fp is the function at stp
    // c         On exit fp is unchanged.
    // c
    // c       dp is a double precision variable.
    // c         On entry dp is the the derivative of the function at stp.
    // c         On exit dp is unchanged.
    // c
    // c       brackt is an logical variable.
    // c         On entry brackt specifies if a minimizer has been bracketed.
    // c            Initially brackt must be set to .false.
    // c         On exit brackt specifies if a minimizer has been bracketed.
    // c            When a minimizer is bracketed brackt is set to .true.
    // c
    // c       stpmin is a double precision variable.
    // c         On entry stpmin is a lower bound for the step.
    // c         On exit stpmin is unchanged.
    // c
    // c       stpmax is a double precision variable.
    // c         On entry stpmax is an upper bound for the step.
    // c         On exit stpmax is unchanged.
    // c
    // c     MINPACK-1 Project. June 1983
    // c     Argonne National Laboratory.
    // c     Jorge J. More' and David J. Thuente.
    // c
    // c     MINPACK-2 Project. November 1993.
    // c     Argonne National Laboratory and University of Minnesota.
    // c     Brett M. Averick and Jorge J. More'.
    // c
    // c     **********
    double zero = 0.;
    double p66 = 0.66;
    double two = 2.;
    double three = 3.;

    double gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta;

    sgnd = dp * (dx / std::abs(dx));

    // c     First case: A higher function value. The minimum is bracketed.
    // c     If the cubic step is closer to stx than the quadratic step, the
    // c     cubic step is taken, otherwise the average of the cubic and
    // c     quadratic steps is taken.

    if (fp > fx)
    {
        theta = three * (fx - fp) / (stp - stx) + dx + dp;
        double temps = std::max(std::abs(theta), std::abs(dx)); // get max(std::abs(theta),std::abs(dx),std::abs(dp))
        s = std::max(temps, std::abs(dp));
        gamma = s * sqrt(pow(theta / s, 2) - (dx / s) * (dp / s));
        if (stp < stx)
            gamma = -gamma;
        p = (gamma - dx) + theta;
        q = ((gamma - dx) + gamma) + dp;
        r = p / q;
        stpc = stx + r * (stp - stx);
        stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / two) * (stp - stx);
        if (std::abs(stpc - stx) < std::abs(stpq - stx))
        {
            stpf = stpc;
        }
        else
        {
            stpf = stpc + (stpq - stpc) / two;
        }
        brackt = true;
    }

    // c     Second case: A lower function value and derivatives of opposite
    // c     sign. The minimum is bracketed. If the cubic step is farther from
    // c     stp than the secant step, the cubic step is taken, otherwise the
    // c     secant step is taken.

    else if (sgnd < zero)
    {
        theta = three * (fx - fp) / (stp - stx) + dx + dp;
        double temps = std::max(std::abs(theta), std::abs(dx)); // get max(std::abs(theta),std::abs(dx),std::abs(dp))
        s = std::max(temps, std::abs(dp));
        gamma = s * sqrt(pow(theta / s, 2) - (dx / s) * (dp / s));
        if (stp > stx)
            gamma = -gamma;
        p = (gamma - dp) + theta;
        q = ((gamma - dp) + gamma) + dx;
        r = p / q;
        stpc = stp + r * (stx - stp);
        stpq = stp + (dp / (dp - dx)) * (stx - stp);
        if (std::abs(stpc - stp) > std::abs(stpq - stp))
        {
            stpf = stpc;
        }
        else
        {
            stpf = stpq;
        }
        brackt = true;
    }

    // c     Third case: A lower function value, derivatives of the same sign,
    // c     and the magnitude of the derivative decreases.

    else if (std::abs(dp) < std::abs(dx))
    {
        // c        The cubic step is computed only if the cubic tends to infinity
        // c        in the direction of the step or if the minimum of the cubic
        // c        is beyond stp. Otherwise the cubic step is defined to be the
        // c        secant step.
        theta = three * (fx - fp) / (stp - stx) + dx + dp;
        double temps = std::max(std::abs(theta), std::abs(dx)); // get max(std::abs(theta),std::abs(dx),std::abs(dp))
        s = std::max(temps, std::abs(dp));
        // c        The case gamma = 0 only arises if the cubic does not tend
        // c        to infinity in the direction of the step.
        gamma = s * sqrt(std::max(zero, pow(theta / s, 2) - (dx / s) * (dp / s)));
        if (stp > stx)
            gamma = -gamma;
        p = (gamma - dp) + theta;
        q = (gamma + (dx - dp)) + gamma;
        r = p / q;
        if (r < zero && gamma != zero)
        {
            stpc = stp + r * (stx - stp);
        }
        else if (stp > stx)
        {
            stpc = stpmax;
        }
        else
        {
            stpc = stpmin;
        }
        stpq = stp + (dp / (dp - dx)) * (stx - stp);

        if (brackt)
        {
            // c           A minimizer has been bracketed. If the cubic step is
            // c           closer to stp than the secant step, the cubic step is
            // c           taken, otherwise the secant step is taken.
            if (std::abs(stpc - stp) < std::abs(stpq - stp))
            {
                stpf = stpc;
            }
            else
            {
                stpf = stpq;
            }
            if (stp > stx)
            {
                stpf = std::min(stp + p66 * (sty - stp), stpf);
            }
            else
            {
                stpf = std::max(stp + p66 * (sty - stp), stpf);
            }
        }
        else
        {
            // c           A minimizer has not been bracketed. If the cubic step is
            // c           farther from stp than the secant step, the cubic step is
            // c           taken, otherwise the secant step is taken.
            if (std::abs(stpc - stp) > std::abs(stpq - stp))
            {
                stpf = stpc;
            }
            else
            {
                stpf = stpq;
            }
            stpf = std::min(stpmax, stpf);
            stpf = std::max(stpmin, stpf);
        }
    }
    // c     Fourth case: A lower function value, derivatives of the same sign,
    // c     and the magnitude of the derivative does not decrease. If the
    // c     minimum is not bracketed, the step is either stpmin or stpmax,
    // c     otherwise the cubic step is taken.
    else
    {
        if (brackt)
        {
            theta = three * (fp - fy) / (sty - stp) + dy + dp;
            double temps
                = std::max(std::abs(theta), std::abs(dy)); // get max(std::abs(theta),std::abs(dy),std::abs(dp))
            s = std::max(temps, std::abs(dp));
            gamma = s * sqrt(pow(theta / s, 2) - (dy / s) * (dp / s));
            if (stp > sty)
                gamma = -gamma;
            p = (gamma - dp) + theta;
            q = ((gamma - dp) + gamma) + dy;
            r = p / q;
            stpc = stp + r * (sty - stp);
            stpf = stpc;
        }
        else if (stp > stx)
        {
            stpf = stpmax;
        }
        else
        {
            stpf = stpmin;
        }
    }
    // c     Update the interval which contains a minimizer.
    if (fp > fx)
    {
        sty = stp;
        fy = fp;
        dy = dp;
    }
    else
    {
        if (sgnd < zero)
        {
            sty = stx;
            fy = fx;
            dy = dx;
        }
        stx = stp;
        fx = fp;
        dx = dp;
    }
    // c     Compute the new step.
    stp = stpf;
}

void Opt_DCsrch::dcSrch(double& f, double& g, double& rstp, char* rtask)
{
    dcsrch(rstp,
           f,
           g,
           this->ftol_,
           this->gtol_,
           this->xtol_,
           rtask,
           this->stpmin_,
           this->stpmax_,
           this->isave_,
           this->dsave_);
}
} // namespace ModuleBase