/*************************************************************************
Copyright (c) 2005-2009, Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************/

#ifndef _autogk_h
#define _autogk_h

#include "ap.h"
#include "amp.h"
#include "tsort.h"
#include "hblas.h"
#include "reflections.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "hsschur.h"
#include "evd.h"
#include "gammafunc.h"
#include "gq.h"
#include "gkq.h"
namespace autogk
{
    /*************************************************************************
    Integration report:
    * TerminationType = completetion code:
        * -5    non-convergence of Gauss-Kronrod nodes
                calculation subroutine.
        * -1    incorrect parameters were specified
        *  1    OK
    * Rep.NFEV countains number of function calculations
    * Rep.NIntervals contains number of intervals [a,b]
      was partitioned into.
    *************************************************************************/
    template<unsigned int Precision>
    class autogkreport
    {
    public:
        int terminationtype;
        int nfev;
        int nintervals;
    };


    template<unsigned int Precision>
    class autogkinternalstate
    {
    public:
        amp::ampf<Precision> a;
        amp::ampf<Precision> b;
        amp::ampf<Precision> eps;
        amp::ampf<Precision> xwidth;
        amp::ampf<Precision> x;
        amp::ampf<Precision> f;
        int info;
        amp::ampf<Precision> r;
        ap::template_2d_array< amp::ampf<Precision> > heap;
        int heapsize;
        int heapwidth;
        int heapused;
        amp::ampf<Precision> sumerr;
        amp::ampf<Precision> sumabs;
        ap::template_1d_array< amp::ampf<Precision> > qn;
        ap::template_1d_array< amp::ampf<Precision> > wg;
        ap::template_1d_array< amp::ampf<Precision> > wk;
        ap::template_1d_array< amp::ampf<Precision> > wr;
        int n;
        amp::rcommstate<Precision> rstate;
    };


    /*************************************************************************
    This structure stores internal state of the integration algorithm  between
    subsequent calls of the AutoGKIteration() subroutine.
    *************************************************************************/
    template<unsigned int Precision>
    class autogkstate
    {
    public:
        amp::ampf<Precision> a;
        amp::ampf<Precision> b;
        amp::ampf<Precision> alpha;
        amp::ampf<Precision> beta;
        amp::ampf<Precision> xwidth;
        amp::ampf<Precision> x;
        amp::ampf<Precision> xminusa;
        amp::ampf<Precision> bminusx;
        amp::ampf<Precision> f;
        int wrappermode;
        autogkinternalstate<Precision> internalstate;
        amp::rcommstate<Precision> rstate;
        amp::ampf<Precision> v;
        int terminationtype;
        int nfev;
        int nintervals;
    };




    template<unsigned int Precision>
    void autogksmooth(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        autogkstate<Precision>& state);
    template<unsigned int Precision>
    void autogksmoothw(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        amp::ampf<Precision> xwidth,
        autogkstate<Precision>& state);
    template<unsigned int Precision>
    void autogksingular(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        amp::ampf<Precision> alpha,
        amp::ampf<Precision> beta,
        autogkstate<Precision>& state);
    template<unsigned int Precision>
    bool autogkiteration(autogkstate<Precision>& state);
    template<unsigned int Precision>
    void autogkresults(const autogkstate<Precision>& state,
        amp::ampf<Precision>& v,
        autogkreport<Precision>& rep);
    template<unsigned int Precision>
    void autogkinternalprepare(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        amp::ampf<Precision> eps,
        amp::ampf<Precision> xwidth,
        autogkinternalstate<Precision>& state);
    template<unsigned int Precision>
    bool autogkinternaliteration(autogkinternalstate<Precision>& state);
    template<unsigned int Precision>
    void mheappop(ap::template_2d_array< amp::ampf<Precision> >& heap,
        int heapsize,
        int heapwidth);
    template<unsigned int Precision>
    void mheappush(ap::template_2d_array< amp::ampf<Precision> >& heap,
        int heapsize,
        int heapwidth);
    template<unsigned int Precision>
    void mheapresize(ap::template_2d_array< amp::ampf<Precision> >& heap,
        int& heapsize,
        int newheapsize,
        int heapwidth);


    /*************************************************************************
    Integration of a smooth function F(x) on a finite interval [a,b].

    Fast-convergent algorithm based on a Gauss-Kronrod formula is used. Result
    is calculated with accuracy close to the machine precision.

    Algorithm works well only with smooth integrands.  It  may  be  used  with
    continuous non-smooth integrands, but with  less  performance.

    It should never be used with integrands which have integrable singularities
    at lower or upper limits - algorithm may crash. Use AutoGKSingular in such
    cases.

    INPUT PARAMETERS:
        A, B    -   interval boundaries (A<B, A=B or A>B)
        
    OUTPUT PARAMETERS
        State   -   structure which stores algorithm state between  subsequent
                    calls of AutoGKIteration.  Used for reverse communication.
                    This structure should be  passed  to  the  AutoGKIteration
                    subroutine.

    SEE ALSO
        AutoGKSmoothW, AutoGKSingular, AutoGKIteration, AutoGKResults.
        

      -- ALGLIB --
         Copyright 06.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void autogksmooth(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        autogkstate<Precision>& state)
    {
        autogksmoothw<Precision>(a, b, amp::ampf<Precision>("0.0"), state);
    }


    /*************************************************************************
    Integration of a smooth function F(x) on a finite interval [a,b].

    This subroutine is same as AutoGKSmooth(), but it guarantees that interval
    [a,b] is partitioned into subintervals which have width at most XWidth.

    Subroutine  can  be  used  when  integrating nearly-constant function with
    narrow "bumps" (about XWidth wide). If "bumps" are too narrow, AutoGKSmooth
    subroutine can overlook them.

    INPUT PARAMETERS:
        A, B    -   interval boundaries (A<B, A=B or A>B)

    OUTPUT PARAMETERS
        State   -   structure which stores algorithm state between  subsequent
                    calls of AutoGKIteration.  Used for reverse communication.
                    This structure should be  passed  to  the  AutoGKIteration
                    subroutine.

    SEE ALSO
        AutoGKSmooth, AutoGKSingular, AutoGKIteration, AutoGKResults.


      -- ALGLIB --
         Copyright 06.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void autogksmoothw(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        amp::ampf<Precision> xwidth,
        autogkstate<Precision>& state)
    {
        state.wrappermode = 0;
        state.a = a;
        state.b = b;
        state.xwidth = xwidth;
        state.rstate.ra.setbounds(0, 10);
        state.rstate.stage = -1;
    }


    /*************************************************************************
    Integration on a finite interval [A,B].
    Integrand have integrable singularities at A/B.

    F(X) must diverge as "(x-A)^alpha" at A, as "(B-x)^beta" at B,  with known
    alpha/beta (alpha>-1, beta>-1).  If alpha/beta  are  not known,  estimates
    from below can be used (but these estimates should be greater than -1 too).

    One  of  alpha/beta variables (or even both alpha/beta) may be equal to 0,
    which means than function F(x) is non-singular at A/B. Anyway (singular at
    bounds or not), function F(x) is supposed to be continuous on (A,B).

    Fast-convergent algorithm based on a Gauss-Kronrod formula is used. Result
    is calculated with accuracy close to the machine precision.

    INPUT PARAMETERS:
        A, B    -   interval boundaries (A<B, A=B or A>B)
        Alpha   -   power-law coefficient of the F(x) at A,
                    Alpha>-1
        Beta    -   power-law coefficient of the F(x) at B,
                    Beta>-1

    OUTPUT PARAMETERS
        State   -   structure which stores algorithm state between  subsequent
                    calls of AutoGKIteration.  Used for reverse communication.
                    This structure should be  passed  to  the  AutoGKIteration
                    subroutine.

    SEE ALSO
        AutoGKSmooth, AutoGKSmoothW, AutoGKIteration, AutoGKResults.


      -- ALGLIB --
         Copyright 06.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void autogksingular(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        amp::ampf<Precision> alpha,
        amp::ampf<Precision> beta,
        autogkstate<Precision>& state)
    {
        state.wrappermode = 1;
        state.a = a;
        state.b = b;
        state.alpha = alpha;
        state.beta = beta;
        state.xwidth = amp::ampf<Precision>("0.0");
        state.rstate.ra.setbounds(0, 10);
        state.rstate.stage = -1;
    }


    /*************************************************************************
    One step of adaptive integration process.

    Called after initialization with one of AutoGKXXX subroutines.
    See HTML documentation for examples.

    Input parameters:
        State   -   structure which stores algorithm state between  calls  and
                    which  is  used  for   reverse   communication.   Must  be
                    initialized with one of AutoGKXXX subroutines.

    If suborutine returned False, iterative proces has converged. If subroutine
    returned True, caller should calculate function value State.F  at  State.X
    and call AutoGKIteration again.

    NOTE:

    When integrating "difficult" functions with integrable singularities like

        F(x) = (x-A)^alpha * (B-x)^beta

    subroutine may require the value of F at points which are too close to A/B.
    Sometimes to calculate integral with high enough precision we  may need to
    calculate F(A+delta) when delta is less than machine  epsilon.  In  finite
    precision arithmetics A+delta will be effectively equal to A,  so  we  may
    find us in situation when  we  are  trying  to  calculate  something  like
    1/sqrt(1-1).

    To avoid  such  situations,  AutoGKIteration  subroutine  fills  not  only
    State.X  field,  but  also   State.XMinusA   (which  equals  to  X-A)  and
    State.BMinusX  (which  equals to B-X) fields.  If X is too close to A or B
    (X-A<0.001*A, or B-X<0.001*B, for example) use  these  fields  instead  of
    State.X


      -- ALGLIB --
         Copyright 07.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool autogkiteration(autogkstate<Precision>& state)
    {
        bool result;
        amp::ampf<Precision> s;
        amp::ampf<Precision> tmp;
        amp::ampf<Precision> eps;
        amp::ampf<Precision> a;
        amp::ampf<Precision> b;
        amp::ampf<Precision> x;
        amp::ampf<Precision> t;
        amp::ampf<Precision> alpha;
        amp::ampf<Precision> beta;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v2;


        
        //
        // Reverse communication preparations
        // I know it looks ugly, but it works the same way
        // anywhere from C++ to Python.
        //
        // This code initializes locals by:
        // * random values determined during code
        //   generation - on first subroutine call
        // * values from previous call - on subsequent calls
        //
        if( state.rstate.stage>=0 )
        {
            s = state.rstate.ra(0);
            tmp = state.rstate.ra(1);
            eps = state.rstate.ra(2);
            a = state.rstate.ra(3);
            b = state.rstate.ra(4);
            x = state.rstate.ra(5);
            t = state.rstate.ra(6);
            alpha = state.rstate.ra(7);
            beta = state.rstate.ra(8);
            v1 = state.rstate.ra(9);
            v2 = state.rstate.ra(10);
        }
        else
        {
            s = -983;
            tmp = -989;
            eps = -834;
            a = 900;
            b = -287;
            x = 364;
            t = 214;
            alpha = -338;
            beta = -686;
            v1 = 912;
            v2 = 585;
        }
        if( state.rstate.stage==0 )
        {
            goto lbl_0;
        }
        if( state.rstate.stage==1 )
        {
            goto lbl_1;
        }
        if( state.rstate.stage==2 )
        {
            goto lbl_2;
        }
        
        //
        // Routine body
        //
        eps = 0;
        a = state.a;
        b = state.b;
        alpha = state.alpha;
        beta = state.beta;
        state.terminationtype = -1;
        state.nfev = 0;
        state.nintervals = 0;
        
        //
        // smooth function  at a finite interval
        //
        if( state.wrappermode!=0 )
        {
            goto lbl_3;
        }
        
        //
        // special case
        //
        if( a==b )
        {
            state.terminationtype = 1;
            state.v = 0;
            result = false;
            return result;
        }
        
        //
        // general case
        //
        autogkinternalprepare<Precision>(a, b, eps, state.xwidth, state.internalstate);
    lbl_5:
        if( ! autogkinternaliteration<Precision>(state.internalstate) )
        {
            goto lbl_6;
        }
        x = state.internalstate.x;
        state.x = x;
        state.xminusa = x-a;
        state.bminusx = b-x;
        state.rstate.stage = 0;
        goto lbl_rcomm;
    lbl_0:
        state.nfev = state.nfev+1;
        state.internalstate.f = state.f;
        goto lbl_5;
    lbl_6:
        state.v = state.internalstate.r;
        state.terminationtype = state.internalstate.info;
        state.nintervals = state.internalstate.heapused;
        result = false;
        return result;
    lbl_3:
        
        //
        // function with power-law singularities at the ends of a finite interval
        //
        if( state.wrappermode!=1 )
        {
            goto lbl_7;
        }
        
        //
        // test coefficients
        //
        if( alpha<=-1 || beta<=-1 )
        {
            state.terminationtype = -1;
            state.v = 0;
            result = false;
            return result;
        }
        
        //
        // special cases
        //
        if( a==b )
        {
            state.terminationtype = 1;
            state.v = 0;
            result = false;
            return result;
        }
        
        //
        // reduction to general form
        //
        if( a<b )
        {
            s = +1;
        }
        else
        {
            s = -1;
            tmp = a;
            a = b;
            b = tmp;
            tmp = alpha;
            alpha = beta;
            beta = tmp;
        }
        alpha = amp::minimum<Precision>(alpha, amp::ampf<Precision>(0));
        beta = amp::minimum<Precision>(beta, amp::ampf<Precision>(0));
        
        //
        // first, integrate left half of [a,b]:
        //     integral(f(x)dx, a, (b+a)/2) =
        //     = 1/(1+alpha) * integral(t^(-alpha/(1+alpha))*f(a+t^(1/(1+alpha)))dt, 0, (0.5*(b-a))^(1+alpha))
        //
        autogkinternalprepare<Precision>(amp::ampf<Precision>(0), amp::pow<Precision>(amp::ampf<Precision>("0.5")*(b-a), 1+alpha), eps, state.xwidth, state.internalstate);
    lbl_9:
        if( ! autogkinternaliteration<Precision>(state.internalstate) )
        {
            goto lbl_10;
        }
        
        //
        // Fill State.X, State.XMinusA, State.BMinusX.
        // Latter two are filled correctly even if B<A.
        //
        x = state.internalstate.x;
        t = amp::pow<Precision>(x, 1/(1+alpha));
        state.x = a+t;
        if( s>0 )
        {
            state.xminusa = t;
            state.bminusx = b-(a+t);
        }
        else
        {
            state.xminusa = a+t-b;
            state.bminusx = -t;
        }
        state.rstate.stage = 1;
        goto lbl_rcomm;
    lbl_1:
        if( alpha!=0 )
        {
            state.internalstate.f = state.f*amp::pow<Precision>(x, -alpha/(1+alpha))/(1+alpha);
        }
        else
        {
            state.internalstate.f = state.f;
        }
        state.nfev = state.nfev+1;
        goto lbl_9;
    lbl_10:
        v1 = state.internalstate.r;
        state.nintervals = state.nintervals+state.internalstate.heapused;
        
        //
        // then, integrate right half of [a,b]:
        //     integral(f(x)dx, (b+a)/2, b) =
        //     = 1/(1+beta) * integral(t^(-beta/(1+beta))*f(b-t^(1/(1+beta)))dt, 0, (0.5*(b-a))^(1+beta))
        //
        autogkinternalprepare<Precision>(amp::ampf<Precision>(0), amp::pow<Precision>(amp::ampf<Precision>("0.5")*(b-a), 1+beta), eps, state.xwidth, state.internalstate);
    lbl_11:
        if( ! autogkinternaliteration<Precision>(state.internalstate) )
        {
            goto lbl_12;
        }
        
        //
        // Fill State.X, State.XMinusA, State.BMinusX.
        // Latter two are filled correctly (X-A, B-X) even if B<A.
        //
        x = state.internalstate.x;
        t = amp::pow<Precision>(x, 1/(1+beta));
        state.x = b-t;
        if( s>0 )
        {
            state.xminusa = b-t-a;
            state.bminusx = t;
        }
        else
        {
            state.xminusa = -t;
            state.bminusx = a-(b-t);
        }
        state.rstate.stage = 2;
        goto lbl_rcomm;
    lbl_2:
        if( beta!=0 )
        {
            state.internalstate.f = state.f*amp::pow<Precision>(x, -beta/(1+beta))/(1+beta);
        }
        else
        {
            state.internalstate.f = state.f;
        }
        state.nfev = state.nfev+1;
        goto lbl_11;
    lbl_12:
        v2 = state.internalstate.r;
        state.nintervals = state.nintervals+state.internalstate.heapused;
        
        //
        // final result
        //
        state.v = s*(v1+v2);
        state.terminationtype = 1;
        result = false;
        return result;
    lbl_7:
        result = false;
        return result;
        
        //
        // Saving state
        //
    lbl_rcomm:
        result = true;
        state.rstate.ra(0) = s;
        state.rstate.ra(1) = tmp;
        state.rstate.ra(2) = eps;
        state.rstate.ra(3) = a;
        state.rstate.ra(4) = b;
        state.rstate.ra(5) = x;
        state.rstate.ra(6) = t;
        state.rstate.ra(7) = alpha;
        state.rstate.ra(8) = beta;
        state.rstate.ra(9) = v1;
        state.rstate.ra(10) = v2;
        return result;
    }


    /*************************************************************************
    Adaptive integration results

    Called after AutoGKIteration returned False.

    Input parameters:
        State   -   algorithm state (used by AutoGKIteration).

    Output parameters:
        V       -   integral(f(x)dx,a,b)
        Rep     -   optimization report (see AutoGKReport description)

      -- ALGLIB --
         Copyright 14.11.2007 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void autogkresults(const autogkstate<Precision>& state,
        amp::ampf<Precision>& v,
        autogkreport<Precision>& rep)
    {
        v = state.v;
        rep.terminationtype = state.terminationtype;
        rep.nfev = state.nfev;
        rep.nintervals = state.nintervals;
    }


    /*************************************************************************
    Internal AutoGK subroutine
    eps<0   - error
    eps=0   - automatic eps selection

    width<0 -   error
    width=0 -   no width requirements
    *************************************************************************/
    template<unsigned int Precision>
    void autogkinternalprepare(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        amp::ampf<Precision> eps,
        amp::ampf<Precision> xwidth,
        autogkinternalstate<Precision>& state)
    {
        
        //
        // Save settings
        //
        state.a = a;
        state.b = b;
        state.eps = eps;
        state.xwidth = xwidth;
        
        //
        // Prepare RComm structure
        //
        state.rstate.ia.setbounds(0, 3);
        state.rstate.ra.setbounds(0, 8);
        state.rstate.stage = -1;
    }


    /*************************************************************************
    Internal AutoGK subroutine
    *************************************************************************/
    template<unsigned int Precision>
    bool autogkinternaliteration(autogkinternalstate<Precision>& state)
    {
        bool result;
        amp::ampf<Precision> c1;
        amp::ampf<Precision> c2;
        int i;
        int j;
        amp::ampf<Precision> intg;
        amp::ampf<Precision> intk;
        amp::ampf<Precision> inta;
        amp::ampf<Precision> v;
        amp::ampf<Precision> ta;
        amp::ampf<Precision> tb;
        int ns;
        amp::ampf<Precision> qeps;
        int info;


        
        //
        // Reverse communication preparations
        // I know it looks ugly, but it works the same way
        // anywhere from C++ to Python.
        //
        // This code initializes locals by:
        // * random values determined during code
        //   generation - on first subroutine call
        // * values from previous call - on subsequent calls
        //
        if( state.rstate.stage>=0 )
        {
            i = state.rstate.ia(0);
            j = state.rstate.ia(1);
            ns = state.rstate.ia(2);
            info = state.rstate.ia(3);
            c1 = state.rstate.ra(0);
            c2 = state.rstate.ra(1);
            intg = state.rstate.ra(2);
            intk = state.rstate.ra(3);
            inta = state.rstate.ra(4);
            v = state.rstate.ra(5);
            ta = state.rstate.ra(6);
            tb = state.rstate.ra(7);
            qeps = state.rstate.ra(8);
        }
        else
        {
            i = 497;
            j = -271;
            ns = -581;
            info = 745;
            c1 = -533;
            c2 = -77;
            intg = 678;
            intk = -293;
            inta = 316;
            v = 647;
            ta = -756;
            tb = 830;
            qeps = -871;
        }
        if( state.rstate.stage==0 )
        {
            goto lbl_0;
        }
        if( state.rstate.stage==1 )
        {
            goto lbl_1;
        }
        if( state.rstate.stage==2 )
        {
            goto lbl_2;
        }
        
        //
        // Routine body
        //
        
        //
        // initialize quadratures.
        // use 15-point Gauss-Kronrod formula.
        //
        state.n = 15;
        gkq::gkqgenerategausslegendre<Precision>(state.n, info, state.qn, state.wk, state.wg);
        if( info<0 )
        {
            state.info = -5;
            state.r = 0;
            result = false;
            return result;
        }
        state.wr.setlength(state.n);
        for(i=0; i<=state.n-1; i++)
        {
            if( i==0 )
            {
                state.wr(i) = amp::ampf<Precision>("0.5")*amp::abs<Precision>(state.qn(1)-state.qn(0));
                continue;
            }
            if( i==state.n-1 )
            {
                state.wr(state.n-1) = amp::ampf<Precision>("0.5")*amp::abs<Precision>(state.qn(state.n-1)-state.qn(state.n-2));
                continue;
            }
            state.wr(i) = amp::ampf<Precision>("0.5")*amp::abs<Precision>(state.qn(i-1)-state.qn(i+1));
        }
        
        //
        // special case
        //
        if( state.a==state.b )
        {
            state.info = 1;
            state.r = 0;
            result = false;
            return result;
        }
        
        //
        // test parameters
        //
        if( state.eps<0 || state.xwidth<0 )
        {
            state.info = -1;
            state.r = 0;
            result = false;
            return result;
        }
        state.info = 1;
        if( state.eps==0 )
        {
            state.eps = 1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        }
        
        //
        // First, prepare heap
        // * column 0   -   absolute error
        // * column 1   -   integral of a F(x) (calculated using Kronrod extension nodes)
        // * column 2   -   integral of a |F(x)| (calculated using modified rect. method)
        // * column 3   -   left boundary of a subinterval
        // * column 4   -   right boundary of a subinterval
        //
        if( state.xwidth!=0 )
        {
            goto lbl_3;
        }
        
        //
        // no maximum width requirements
        // start from one big subinterval
        //
        state.heapwidth = 5;
        state.heapsize = 1;
        state.heapused = 1;
        state.heap.setlength(state.heapsize, state.heapwidth);
        c1 = amp::ampf<Precision>("0.5")*(state.b-state.a);
        c2 = amp::ampf<Precision>("0.5")*(state.b+state.a);
        intg = 0;
        intk = 0;
        inta = 0;
        i = 0;
    lbl_5:
        if( i>state.n-1 )
        {
            goto lbl_7;
        }
        
        //
        // obtain F
        //
        state.x = c1*state.qn(i)+c2;
        state.rstate.stage = 0;
        goto lbl_rcomm;
    lbl_0:
        v = state.f;
        
        //
        // Gauss-Kronrod formula
        //
        intk = intk+v*state.wk(i);
        if( i%2==1 )
        {
            intg = intg+v*state.wg(i);
        }
        
        //
        // Integral |F(x)|
        // Use rectangles method
        //
        inta = inta+amp::abs<Precision>(v)*state.wr(i);
        i = i+1;
        goto lbl_5;
    lbl_7:
        intk = intk*(state.b-state.a)*amp::ampf<Precision>("0.5");
        intg = intg*(state.b-state.a)*amp::ampf<Precision>("0.5");
        inta = inta*(state.b-state.a)*amp::ampf<Precision>("0.5");
        state.heap(0,0) = amp::abs<Precision>(intg-intk);
        state.heap(0,1) = intk;
        state.heap(0,2) = inta;
        state.heap(0,3) = state.a;
        state.heap(0,4) = state.b;
        state.sumerr = state.heap(0,0);
        state.sumabs = amp::abs<Precision>(inta);
        goto lbl_4;
    lbl_3:
        
        //
        // maximum subinterval should be no more than XWidth.
        // so we create Ceil((B-A)/XWidth)+1 small subintervals
        //
        ns = amp::ceil<Precision>(amp::abs<Precision>(state.b-state.a)/state.xwidth)+1;
        state.heapsize = ns;
        state.heapused = ns;
        state.heapwidth = 5;
        state.heap.setlength(state.heapsize, state.heapwidth);
        state.sumerr = 0;
        state.sumabs = 0;
        j = 0;
    lbl_8:
        if( j>ns-1 )
        {
            goto lbl_10;
        }
        ta = state.a+j*(state.b-state.a)/ns;
        tb = state.a+(j+1)*(state.b-state.a)/ns;
        c1 = amp::ampf<Precision>("0.5")*(tb-ta);
        c2 = amp::ampf<Precision>("0.5")*(tb+ta);
        intg = 0;
        intk = 0;
        inta = 0;
        i = 0;
    lbl_11:
        if( i>state.n-1 )
        {
            goto lbl_13;
        }
        
        //
        // obtain F
        //
        state.x = c1*state.qn(i)+c2;
        state.rstate.stage = 1;
        goto lbl_rcomm;
    lbl_1:
        v = state.f;
        
        //
        // Gauss-Kronrod formula
        //
        intk = intk+v*state.wk(i);
        if( i%2==1 )
        {
            intg = intg+v*state.wg(i);
        }
        
        //
        // Integral |F(x)|
        // Use rectangles method
        //
        inta = inta+amp::abs<Precision>(v)*state.wr(i);
        i = i+1;
        goto lbl_11;
    lbl_13:
        intk = intk*(tb-ta)*amp::ampf<Precision>("0.5");
        intg = intg*(tb-ta)*amp::ampf<Precision>("0.5");
        inta = inta*(tb-ta)*amp::ampf<Precision>("0.5");
        state.heap(j,0) = amp::abs<Precision>(intg-intk);
        state.heap(j,1) = intk;
        state.heap(j,2) = inta;
        state.heap(j,3) = ta;
        state.heap(j,4) = tb;
        state.sumerr = state.sumerr+state.heap(j,0);
        state.sumabs = state.sumabs+amp::abs<Precision>(inta);
        j = j+1;
        goto lbl_8;
    lbl_10:
    lbl_4:
        
        //
        // method iterations
        //
    lbl_14:
        if( false )
        {
            goto lbl_15;
        }
        
        //
        // additional memory if needed
        //
        if( state.heapused==state.heapsize )
        {
            mheapresize<Precision>(state.heap, state.heapsize, 4*state.heapsize, state.heapwidth);
        }
        
        //
        // TODO: every 20 iterations recalculate errors/sums
        // TODO: one more criterion to prevent infinite loops with too strict Eps
        //
        if( state.sumerr<=state.eps*state.sumabs )
        {
            state.r = 0;
            for(j=0; j<=state.heapused-1; j++)
            {
                state.r = state.r+state.heap(j,1);
            }
            result = false;
            return result;
        }
        
        //
        // Exclude interval with maximum absolute error
        //
        mheappop<Precision>(state.heap, state.heapused, state.heapwidth);
        state.sumerr = state.sumerr-state.heap(state.heapused-1,0);
        state.sumabs = state.sumabs-state.heap(state.heapused-1,2);
        
        //
        // Divide interval, create subintervals
        //
        ta = state.heap(state.heapused-1,3);
        tb = state.heap(state.heapused-1,4);
        state.heap(state.heapused-1,3) = ta;
        state.heap(state.heapused-1,4) = amp::ampf<Precision>("0.5")*(ta+tb);
        state.heap(state.heapused,3) = amp::ampf<Precision>("0.5")*(ta+tb);
        state.heap(state.heapused,4) = tb;
        j = state.heapused-1;
    lbl_16:
        if( j>state.heapused )
        {
            goto lbl_18;
        }
        c1 = amp::ampf<Precision>("0.5")*(state.heap(j,4)-state.heap(j,3));
        c2 = amp::ampf<Precision>("0.5")*(state.heap(j,4)+state.heap(j,3));
        intg = 0;
        intk = 0;
        inta = 0;
        i = 0;
    lbl_19:
        if( i>state.n-1 )
        {
            goto lbl_21;
        }
        
        //
        // F(x)
        //
        state.x = c1*state.qn(i)+c2;
        state.rstate.stage = 2;
        goto lbl_rcomm;
    lbl_2:
        v = state.f;
        
        //
        // Gauss-Kronrod formula
        //
        intk = intk+v*state.wk(i);
        if( i%2==1 )
        {
            intg = intg+v*state.wg(i);
        }
        
        //
        // Integral |F(x)|
        // Use rectangles method
        //
        inta = inta+amp::abs<Precision>(v)*state.wr(i);
        i = i+1;
        goto lbl_19;
    lbl_21:
        intk = intk*(state.heap(j,4)-state.heap(j,3))*amp::ampf<Precision>("0.5");
        intg = intg*(state.heap(j,4)-state.heap(j,3))*amp::ampf<Precision>("0.5");
        inta = inta*(state.heap(j,4)-state.heap(j,3))*amp::ampf<Precision>("0.5");
        state.heap(j,0) = amp::abs<Precision>(intg-intk);
        state.heap(j,1) = intk;
        state.heap(j,2) = inta;
        state.sumerr = state.sumerr+state.heap(j,0);
        state.sumabs = state.sumabs+state.heap(j,2);
        j = j+1;
        goto lbl_16;
    lbl_18:
        mheappush<Precision>(state.heap, state.heapused-1, state.heapwidth);
        mheappush<Precision>(state.heap, state.heapused, state.heapwidth);
        state.heapused = state.heapused+1;
        goto lbl_14;
    lbl_15:
        result = false;
        return result;
        
        //
        // Saving state
        //
    lbl_rcomm:
        result = true;
        state.rstate.ia(0) = i;
        state.rstate.ia(1) = j;
        state.rstate.ia(2) = ns;
        state.rstate.ia(3) = info;
        state.rstate.ra(0) = c1;
        state.rstate.ra(1) = c2;
        state.rstate.ra(2) = intg;
        state.rstate.ra(3) = intk;
        state.rstate.ra(4) = inta;
        state.rstate.ra(5) = v;
        state.rstate.ra(6) = ta;
        state.rstate.ra(7) = tb;
        state.rstate.ra(8) = qeps;
        return result;
    }


    template<unsigned int Precision>
    void mheappop(ap::template_2d_array< amp::ampf<Precision> >& heap,
        int heapsize,
        int heapwidth)
    {
        int i;
        int p;
        amp::ampf<Precision> t;
        int maxcp;


        if( heapsize==1 )
        {
            return;
        }
        for(i=0; i<=heapwidth-1; i++)
        {
            t = heap(heapsize-1,i);
            heap(heapsize-1,i) = heap(0,i);
            heap(0,i) = t;
        }
        p = 0;
        while( 2*p+1<heapsize-1 )
        {
            maxcp = 2*p+1;
            if( 2*p+2<heapsize-1 )
            {
                if( heap(2*p+2,0)>heap(2*p+1,0) )
                {
                    maxcp = 2*p+2;
                }
            }
            if( heap(p,0)<heap(maxcp,0) )
            {
                for(i=0; i<=heapwidth-1; i++)
                {
                    t = heap(p,i);
                    heap(p,i) = heap(maxcp,i);
                    heap(maxcp,i) = t;
                }
                p = maxcp;
            }
            else
            {
                break;
            }
        }
    }


    template<unsigned int Precision>
    void mheappush(ap::template_2d_array< amp::ampf<Precision> >& heap,
        int heapsize,
        int heapwidth)
    {
        int i;
        int p;
        amp::ampf<Precision> t;
        int parent;


        if( heapsize==0 )
        {
            return;
        }
        p = heapsize;
        while( p!=0 )
        {
            parent = (p-1)/2;
            if( heap(p,0)>heap(parent,0) )
            {
                for(i=0; i<=heapwidth-1; i++)
                {
                    t = heap(p,i);
                    heap(p,i) = heap(parent,i);
                    heap(parent,i) = t;
                }
                p = parent;
            }
            else
            {
                break;
            }
        }
    }


    template<unsigned int Precision>
    void mheapresize(ap::template_2d_array< amp::ampf<Precision> >& heap,
        int& heapsize,
        int newheapsize,
        int heapwidth)
    {
        ap::template_2d_array< amp::ampf<Precision> > tmp;
        int i;


        tmp.setlength(heapsize, heapwidth);
        for(i=0; i<=heapsize-1; i++)
        {
            amp::vmove(tmp.getrow(i, 0, heapwidth-1), heap.getrow(i, 0, heapwidth-1));
        }
        heap.setlength(newheapsize, heapwidth);
        for(i=0; i<=heapsize-1; i++)
        {
            amp::vmove(heap.getrow(i, 0, heapwidth-1), tmp.getrow(i, 0, heapwidth-1));
        }
        heapsize = newheapsize;
    }
} // namespace

#endif
