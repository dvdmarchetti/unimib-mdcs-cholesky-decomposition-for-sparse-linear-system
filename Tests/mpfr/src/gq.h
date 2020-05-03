/*************************************************************************
Copyright (c) 2005-2007, Sergey Bochkanov (ALGLIB project).

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

#ifndef _gq_h
#define _gq_h

#include "ap.h"
#include "amp.h"
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
namespace gq
{
    template<unsigned int Precision>
    void gqgeneraterec(const ap::template_1d_array< amp::ampf<Precision> >& alpha,
        const ap::template_1d_array< amp::ampf<Precision> >& beta,
        amp::ampf<Precision> mu0,
        int n,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w);
    template<unsigned int Precision>
    void gqgenerategausslobattorec(ap::template_1d_array< amp::ampf<Precision> > alpha,
        ap::template_1d_array< amp::ampf<Precision> > beta,
        amp::ampf<Precision> mu0,
        amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        int n,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w);
    template<unsigned int Precision>
    void gqgenerategaussradaurec(ap::template_1d_array< amp::ampf<Precision> > alpha,
        ap::template_1d_array< amp::ampf<Precision> > beta,
        amp::ampf<Precision> mu0,
        amp::ampf<Precision> a,
        int n,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w);
    template<unsigned int Precision>
    void gqgenerategausslegendre(int n,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w);
    template<unsigned int Precision>
    void gqgenerategaussjacobi(int n,
        amp::ampf<Precision> alpha,
        amp::ampf<Precision> beta,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w);
    template<unsigned int Precision>
    void gqgenerategausslaguerre(int n,
        amp::ampf<Precision> alpha,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w);
    template<unsigned int Precision>
    void gqgenerategausshermite(int n,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w);


    /*************************************************************************
    Computation of nodes and weights for a Gauss quadrature formula

    The algorithm generates the N-point Gauss quadrature formula  with  weight
    function given by coefficients alpha and beta  of  a  recurrence  relation
    which generates a system of orthogonal polynomials:

    P-1(x)   =  0
    P0(x)    =  1
    Pn+1(x)  =  (x-alpha(n))*Pn(x)  -  beta(n)*Pn-1(x)

    and zeroth moment Mu0

    Mu0 = integral(W(x)dx,a,b)

    INPUT PARAMETERS:
        Alpha   �   array[0..N-1], alpha coefficients
        Beta    �   array[0..N-1], beta coefficients
                    Zero-indexed element is not used and may be arbitrary.
                    Beta[I]>0.
        Mu0     �   zeroth moment of the weight function.
        N       �   number of nodes of the quadrature formula, N>=1

    OUTPUT PARAMETERS:
        Info    -   error code:
                    * -3    internal eigenproblem solver hasn't converged
                    * -2    Beta[i]<=0
                    * -1    incorrect N was passed
                    *  1    OK
        X       -   array[0..N-1] - array of quadrature nodes,
                    in ascending order.
        W       -   array[0..N-1] - array of quadrature weights.

      -- ALGLIB --
         Copyright 2005-2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void gqgeneraterec(const ap::template_1d_array< amp::ampf<Precision> >& alpha,
        const ap::template_1d_array< amp::ampf<Precision> >& beta,
        amp::ampf<Precision> mu0,
        int n,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w)
    {
        int i;
        ap::template_1d_array< amp::ampf<Precision> > d;
        ap::template_1d_array< amp::ampf<Precision> > e;
        ap::template_2d_array< amp::ampf<Precision> > z;


        if( n<1 )
        {
            info = -1;
            return;
        }
        info = 1;
        
        //
        // Initialize
        //
        d.setlength(n);
        e.setlength(n);
        for(i=1; i<=n-1; i++)
        {
            d(i-1) = alpha(i-1);
            if( beta(i)<=0 )
            {
                info = -2;
                return;
            }
            e(i-1) = amp::sqrt<Precision>(beta(i));
        }
        d(n-1) = alpha(n-1);
        
        //
        // EVD
        //
        if( !evd::smatrixtdevd<Precision>(d, e, n, 3, z) )
        {
            info = -3;
            return;
        }
        
        //
        // Generate
        //
        x.setlength(n);
        w.setlength(n);
        for(i=1; i<=n; i++)
        {
            x(i-1) = d(i-1);
            w(i-1) = mu0*amp::sqr<Precision>(z(0,i-1));
        }
    }


    /*************************************************************************
    Computation of nodes and weights for a Gauss-Lobatto quadrature formula

    The algorithm generates the N-point Gauss-Lobatto quadrature formula  with
    weight function given by coefficients alpha and beta of a recurrence which
    generates a system of orthogonal polynomials.

    P-1(x)   =  0
    P0(x)    =  1
    Pn+1(x)  =  (x-alpha(n))*Pn(x)  -  beta(n)*Pn-1(x)

    and zeroth moment Mu0

    Mu0 = integral(W(x)dx,a,b)

    INPUT PARAMETERS:
        Alpha   �   array[0..N-2], alpha coefficients
        Beta    �   array[0..N-2], beta coefficients.
                    Zero-indexed element is not used, may be arbitrary.
                    Beta[I]>0
        Mu0     �   zeroth moment of the weighting function.
        A       �   left boundary of the integration interval.
        B       �   right boundary of the integration interval.
        N       �   number of nodes of the quadrature formula, N>=3
                    (including the left and right boundary nodes).

    OUTPUT PARAMETERS:
        Info    -   error code:
                    * -3    internal eigenproblem solver hasn't converged
                    * -2    Beta[i]<=0
                    * -1    incorrect N was passed
                    *  1    OK
        X       -   array[0..N-1] - array of quadrature nodes,
                    in ascending order.
        W       -   array[0..N-1] - array of quadrature weights.

      -- ALGLIB --
         Copyright 2005-2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void gqgenerategausslobattorec(ap::template_1d_array< amp::ampf<Precision> > alpha,
        ap::template_1d_array< amp::ampf<Precision> > beta,
        amp::ampf<Precision> mu0,
        amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        int n,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w)
    {
        int i;
        ap::template_1d_array< amp::ampf<Precision> > d;
        ap::template_1d_array< amp::ampf<Precision> > e;
        ap::template_2d_array< amp::ampf<Precision> > z;
        amp::ampf<Precision> pim1a;
        amp::ampf<Precision> pia;
        amp::ampf<Precision> pim1b;
        amp::ampf<Precision> pib;
        amp::ampf<Precision> t;
        amp::ampf<Precision> a11;
        amp::ampf<Precision> a12;
        amp::ampf<Precision> a21;
        amp::ampf<Precision> a22;
        amp::ampf<Precision> b1;
        amp::ampf<Precision> b2;
        amp::ampf<Precision> alph;
        amp::ampf<Precision> bet;


        if( n<=2 )
        {
            info = -1;
            return;
        }
        info = 1;
        
        //
        // Initialize, D[1:N+1], E[1:N]
        //
        n = n-2;
        d.setlength(n+2);
        e.setlength(n+1);
        for(i=1; i<=n+1; i++)
        {
            d(i-1) = alpha(i-1);
        }
        for(i=1; i<=n; i++)
        {
            if( beta(i)<=0 )
            {
                info = -2;
                return;
            }
            e(i-1) = amp::sqrt<Precision>(beta(i));
        }
        
        //
        // Caclulate Pn(a), Pn+1(a), Pn(b), Pn+1(b)
        //
        beta(0) = 0;
        pim1a = 0;
        pia = 1;
        pim1b = 0;
        pib = 1;
        for(i=1; i<=n+1; i++)
        {
            
            //
            // Pi(a)
            //
            t = (a-alpha(i-1))*pia-beta(i-1)*pim1a;
            pim1a = pia;
            pia = t;
            
            //
            // Pi(b)
            //
            t = (b-alpha(i-1))*pib-beta(i-1)*pim1b;
            pim1b = pib;
            pib = t;
        }
        
        //
        // Calculate alpha'(n+1), beta'(n+1)
        //
        a11 = pia;
        a12 = pim1a;
        a21 = pib;
        a22 = pim1b;
        b1 = a*pia;
        b2 = b*pib;
        if( amp::abs<Precision>(a11)>amp::abs<Precision>(a21) )
        {
            a22 = a22-a12*a21/a11;
            b2 = b2-b1*a21/a11;
            bet = b2/a22;
            alph = (b1-bet*a12)/a11;
        }
        else
        {
            a12 = a12-a22*a11/a21;
            b1 = b1-b2*a11/a21;
            bet = b1/a12;
            alph = (b2-bet*a22)/a21;
        }
        if( bet<0 )
        {
            info = -3;
            return;
        }
        d(n+1) = alph;
        e(n) = amp::sqrt<Precision>(bet);
        
        //
        // EVD
        //
        if( !evd::smatrixtdevd<Precision>(d, e, n+2, 3, z) )
        {
            info = -3;
            return;
        }
        
        //
        // Generate
        //
        x.setlength(n+2);
        w.setlength(n+2);
        for(i=1; i<=n+2; i++)
        {
            x(i-1) = d(i-1);
            w(i-1) = mu0*amp::sqr<Precision>(z(0,i-1));
        }
    }


    /*************************************************************************
    Computation of nodes and weights for a Gauss-Radau quadrature formula

    The algorithm generates the N-point Gauss-Radau  quadrature  formula  with
    weight function given by the coefficients alpha and  beta  of a recurrence
    which generates a system of orthogonal polynomials.

    P-1(x)   =  0
    P0(x)    =  1
    Pn+1(x)  =  (x-alpha(n))*Pn(x)  -  beta(n)*Pn-1(x)

    and zeroth moment Mu0

    Mu0 = integral(W(x)dx,a,b)

    INPUT PARAMETERS:
        Alpha   �   array[0..N-2], alpha coefficients.
        Beta    �   array[0..N-1], beta coefficients
                    Zero-indexed element is not used.
                    Beta[I]>0
        Mu0     �   zeroth moment of the weighting function.
        A       �   left boundary of the integration interval.
        N       �   number of nodes of the quadrature formula, N>=2
                    (including the left boundary node).

    OUTPUT PARAMETERS:
        Info    -   error code:
                    * -3    internal eigenproblem solver hasn't converged
                    * -2    Beta[i]<=0
                    * -1    incorrect N was passed
                    *  1    OK
        X       -   array[0..N-1] - array of quadrature nodes,
                    in ascending order.
        W       -   array[0..N-1] - array of quadrature weights.


      -- ALGLIB --
         Copyright 2005-2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void gqgenerategaussradaurec(ap::template_1d_array< amp::ampf<Precision> > alpha,
        ap::template_1d_array< amp::ampf<Precision> > beta,
        amp::ampf<Precision> mu0,
        amp::ampf<Precision> a,
        int n,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w)
    {
        int i;
        ap::template_1d_array< amp::ampf<Precision> > d;
        ap::template_1d_array< amp::ampf<Precision> > e;
        ap::template_2d_array< amp::ampf<Precision> > z;
        amp::ampf<Precision> polim1;
        amp::ampf<Precision> poli;
        amp::ampf<Precision> t;


        if( n<2 )
        {
            info = -1;
            return;
        }
        info = 1;
        
        //
        // Initialize, D[1:N], E[1:N]
        //
        n = n-1;
        d.setlength(n+1);
        e.setlength(n);
        for(i=1; i<=n; i++)
        {
            d(i-1) = alpha(i-1);
            if( beta(i)<=0 )
            {
                info = -2;
                return;
            }
            e(i-1) = amp::sqrt<Precision>(beta(i));
        }
        
        //
        // Caclulate Pn(a), Pn-1(a), and D[N+1]
        //
        beta(0) = 0;
        polim1 = 0;
        poli = 1;
        for(i=1; i<=n; i++)
        {
            t = (a-alpha(i-1))*poli-beta(i-1)*polim1;
            polim1 = poli;
            poli = t;
        }
        d(n) = a-beta(n)*polim1/poli;
        
        //
        // EVD
        //
        if( !evd::smatrixtdevd<Precision>(d, e, n+1, 3, z) )
        {
            info = -3;
            return;
        }
        
        //
        // Generate
        //
        x.setbounds(0, n);
        w.setbounds(0, n);
        for(i=1; i<=n+1; i++)
        {
            x(i-1) = d(i-1);
            w(i-1) = mu0*amp::sqr<Precision>(z(0,i-1));
        }
    }


    /*************************************************************************
    Returns nodes/weights for Gauss-Legendre quadrature on [-1,1] with N
    nodes.

    INPUT PARAMETERS:
        N           -   number of nodes, >=1

    OUTPUT PARAMETERS:
        Info        -   error code:
                        * -4    an  error   was   detected   when  calculating
                                weights/nodes.  N  is  too  large   to  obtain
                                weights/nodes  with  high   enough   accuracy.
                                Try  to   use   multiple   precision  version.
                        * -3    internal eigenproblem solver hasn't  converged
                        * -1    incorrect N was passed
                        * +1    OK
        X           -   array[0..N-1] - array of quadrature nodes,
                        in ascending order.
        W           -   array[0..N-1] - array of quadrature weights.


      -- ALGLIB --
         Copyright 12.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void gqgenerategausslegendre(int n,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w)
    {
        ap::template_1d_array< amp::ampf<Precision> > alpha;
        ap::template_1d_array< amp::ampf<Precision> > beta;
        int i;


        if( n<1 )
        {
            info = -1;
            return;
        }
        alpha.setlength(n);
        beta.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            alpha(i) = 0;
        }
        beta(0) = 2;
        for(i=1; i<=n-1; i++)
        {
            beta(i) = 1/(4-1/amp::sqr<Precision>(amp::ampf<Precision>(i)));
        }
        gqgeneraterec<Precision>(alpha, beta, beta(0), n, info, x, w);
        
        //
        // test basic properties to detect errors
        //
        if( info>0 )
        {
            if( x(0)<-1 || x(n-1)>+1 )
            {
                info = -4;
            }
            for(i=0; i<=n-2; i++)
            {
                if( x(i)>=x(i+1) )
                {
                    info = -4;
                }
            }
        }
    }


    /*************************************************************************
    Returns  nodes/weights  for  Gauss-Jacobi quadrature on [-1,1] with weight
    function W(x)=Power(1-x,Alpha)*Power(1+x,Beta).

    INPUT PARAMETERS:
        N           -   number of nodes, >=1
        Alpha       -   power-law coefficient, Alpha>-1
        Beta        -   power-law coefficient, Beta>-1

    OUTPUT PARAMETERS:
        Info        -   error code:
                        * -4    an  error  was   detected   when   calculating
                                weights/nodes. Alpha or  Beta  are  too  close
                                to -1 to obtain weights/nodes with high enough
                                accuracy, or, may be, N is too large.  Try  to
                                use multiple precision version.
                        * -3    internal eigenproblem solver hasn't converged
                        * -1    incorrect N/Alpha/Beta was passed
                        * +1    OK
        X           -   array[0..N-1] - array of quadrature nodes,
                        in ascending order.
        W           -   array[0..N-1] - array of quadrature weights.


      -- ALGLIB --
         Copyright 12.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void gqgenerategaussjacobi(int n,
        amp::ampf<Precision> alpha,
        amp::ampf<Precision> beta,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w)
    {
        ap::template_1d_array< amp::ampf<Precision> > a;
        ap::template_1d_array< amp::ampf<Precision> > b;
        amp::ampf<Precision> alpha2;
        amp::ampf<Precision> beta2;
        amp::ampf<Precision> apb;
        amp::ampf<Precision> t;
        int i;
        amp::ampf<Precision> s;


        if( n<1 || alpha<=-1 || beta<=-1 )
        {
            info = -1;
            return;
        }
        a.setlength(n);
        b.setlength(n);
        apb = alpha+beta;
        a(0) = (beta-alpha)/(apb+2);
        t = (apb+1)*amp::log<Precision>(amp::ampf<Precision>(2))+gammafunc::lngamma<Precision>(alpha+1, s)+gammafunc::lngamma<Precision>(beta+1, s)-gammafunc::lngamma<Precision>(apb+2, s);
        if( t>amp::log<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber()) )
        {
            info = -4;
            return;
        }
        b(0) = amp::exp<Precision>(t);
        if( n>1 )
        {
            alpha2 = amp::sqr<Precision>(alpha);
            beta2 = amp::sqr<Precision>(beta);
            a(1) = (beta2-alpha2)/((apb+2)*(apb+4));
            b(1) = 4*(alpha+1)*(beta+1)/((apb+3)*amp::sqr<Precision>(apb+2));
            for(i=2; i<=n-1; i++)
            {
                a(i) = amp::ampf<Precision>("0.25")*(beta2-alpha2)/(i*i*(1+amp::ampf<Precision>("0.5")*apb/i)*(1+amp::ampf<Precision>("0.5")*(apb+2)/i));
                b(i) = amp::ampf<Precision>("0.25")*(1+alpha/i)*(1+beta/i)*(1+apb/i)/((1+amp::ampf<Precision>("0.5")*(apb+1)/i)*(1+amp::ampf<Precision>("0.5")*(apb-1)/i)*amp::sqr<Precision>(1+amp::ampf<Precision>("0.5")*apb/i));
            }
        }
        gqgeneraterec<Precision>(a, b, b(0), n, info, x, w);
        
        //
        // test basic properties to detect errors
        //
        if( info>0 )
        {
            if( x(0)<-1 || x(n-1)>+1 )
            {
                info = -4;
            }
            for(i=0; i<=n-2; i++)
            {
                if( x(i)>=x(i+1) )
                {
                    info = -4;
                }
            }
        }
    }


    /*************************************************************************
    Returns  nodes/weights  for  Gauss-Laguerre  quadrature  on  [0,+inf) with
    weight function W(x)=Power(x,Alpha)*Exp(-x)

    INPUT PARAMETERS:
        N           -   number of nodes, >=1
        Alpha       -   power-law coefficient, Alpha>-1

    OUTPUT PARAMETERS:
        Info        -   error code:
                        * -4    an  error  was   detected   when   calculating
                                weights/nodes. Alpha is too  close  to  -1  to
                                obtain weights/nodes with high enough accuracy
                                or, may  be,  N  is  too  large.  Try  to  use
                                multiple precision version.
                        * -3    internal eigenproblem solver hasn't converged
                        * -1    incorrect N/Alpha was passed
                        * +1    OK
        X           -   array[0..N-1] - array of quadrature nodes,
                        in ascending order.
        W           -   array[0..N-1] - array of quadrature weights.


      -- ALGLIB --
         Copyright 12.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void gqgenerategausslaguerre(int n,
        amp::ampf<Precision> alpha,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w)
    {
        ap::template_1d_array< amp::ampf<Precision> > a;
        ap::template_1d_array< amp::ampf<Precision> > b;
        amp::ampf<Precision> t;
        int i;
        amp::ampf<Precision> s;


        if( n<1 || alpha<=-1 )
        {
            info = -1;
            return;
        }
        a.setlength(n);
        b.setlength(n);
        a(0) = alpha+1;
        t = gammafunc::lngamma<Precision>(alpha+1, s);
        if( t>=amp::log<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber()) )
        {
            info = -4;
            return;
        }
        b(0) = amp::exp<Precision>(t);
        if( n>1 )
        {
            for(i=1; i<=n-1; i++)
            {
                a(i) = 2*i+alpha+1;
                b(i) = i*(i+alpha);
            }
        }
        gqgeneraterec<Precision>(a, b, b(0), n, info, x, w);
        
        //
        // test basic properties to detect errors
        //
        if( info>0 )
        {
            if( x(0)<0 )
            {
                info = -4;
            }
            for(i=0; i<=n-2; i++)
            {
                if( x(i)>=x(i+1) )
                {
                    info = -4;
                }
            }
        }
    }


    /*************************************************************************
    Returns  nodes/weights  for  Gauss-Hermite  quadrature on (-inf,+inf) with
    weight function W(x)=Exp(-x*x)

    INPUT PARAMETERS:
        N           -   number of nodes, >=1

    OUTPUT PARAMETERS:
        Info        -   error code:
                        * -4    an  error  was   detected   when   calculating
                                weights/nodes.  May be, N is too large. Try to
                                use multiple precision version.
                        * -3    internal eigenproblem solver hasn't converged
                        * -1    incorrect N/Alpha was passed
                        * +1    OK
        X           -   array[0..N-1] - array of quadrature nodes,
                        in ascending order.
        W           -   array[0..N-1] - array of quadrature weights.


      -- ALGLIB --
         Copyright 12.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void gqgenerategausshermite(int n,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& w)
    {
        ap::template_1d_array< amp::ampf<Precision> > a;
        ap::template_1d_array< amp::ampf<Precision> > b;
        int i;


        if( n<1 )
        {
            info = -1;
            return;
        }
        a.setlength(n);
        b.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            a(i) = 0;
        }
        b(0) = amp::sqrt<Precision>(4*amp::atan<Precision>(amp::ampf<Precision>(1)));
        if( n>1 )
        {
            for(i=1; i<=n-1; i++)
            {
                b(i) = amp::ampf<Precision>("0.5")*i;
            }
        }
        gqgeneraterec<Precision>(a, b, b(0), n, info, x, w);
        
        //
        // test basic properties to detect errors
        //
        if( info>0 )
        {
            for(i=0; i<=n-2; i++)
            {
                if( x(i)>=x(i+1) )
                {
                    info = -4;
                }
            }
        }
    }
} // namespace

#endif
