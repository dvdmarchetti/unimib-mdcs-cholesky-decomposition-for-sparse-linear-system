
#ifndef _testmincgunit_h
#define _testmincgunit_h

#include <stdio.h>
#include "ap.h"
#include "amp.h"
#include "linmin.h"
#include "mincg.h"
namespace testmincgunit
{
    template<unsigned int Precision>
    bool testmincg(bool silent);
    template<unsigned int Precision>
    void testfunc1(mincg::mincgstate<Precision>& state);
    template<unsigned int Precision>
    void testfunc2(mincg::mincgstate<Precision>& state);
    template<unsigned int Precision>
    void testfunc3(mincg::mincgstate<Precision>& state);
    template<unsigned int Precision>
    bool testmincgunit_test_silent();
    template<unsigned int Precision>
    bool testmincgunit_test();


    template<unsigned int Precision>
    bool testmincg(bool silent)
    {
        bool result;
        bool waserrors;
        bool referror;
        bool eqerror;
        bool linerror1;
        bool linerror2;
        bool converror;
        bool othererrors;
        int n;
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > xe;
        ap::template_1d_array< amp::ampf<Precision> > b;
        ap::template_1d_array< amp::ampf<Precision> > xlast;
        amp::ampf<Precision> fprev;
        amp::ampf<Precision> xprev;
        amp::ampf<Precision> stpmax;
        int i;
        int j;
        amp::ampf<Precision> v;
        ap::template_2d_array< amp::ampf<Precision> > a;
        mincg::mincgstate<Precision> state;
        mincg::mincgreport<Precision> rep;
        int cgtype;


        waserrors = false;
        referror = false;
        linerror1 = false;
        linerror2 = false;
        eqerror = false;
        converror = false;
        othererrors = false;
        for(cgtype=0; cgtype<=1; cgtype++)
        {
            
            //
            // Reference problem
            //
            x.setbounds(0, 2);
            n = 3;
            x(0) = 100*amp::ampf<Precision>::getRandom()-50;
            x(1) = 100*amp::ampf<Precision>::getRandom()-50;
            x(2) = 100*amp::ampf<Precision>::getRandom()-50;
            mincg::mincgcreate<Precision>(n, x, state);
            mincg::mincgsetcgtype<Precision>(state, cgtype);
            while( mincg::mincgiteration<Precision>(state) )
            {
                state.f = amp::sqr<Precision>(state.x(0)-2)+amp::sqr<Precision>(state.x(1))+amp::sqr<Precision>(state.x(2)-state.x(0));
                state.g(0) = 2*(state.x(0)-2)+2*(state.x(0)-state.x(2));
                state.g(1) = 2*state.x(1);
                state.g(2) = 2*(state.x(2)-state.x(0));
            }
            mincg::mincgresults<Precision>(state, x, rep);
            referror = referror || rep.terminationtype<=0 || amp::abs<Precision>(x(0)-2)>amp::ampf<Precision>("0.001") || amp::abs<Precision>(x(1))>amp::ampf<Precision>("0.001") || amp::abs<Precision>(x(2)-2)>amp::ampf<Precision>("0.001");
            
            //
            // 1D problem #1
            //
            x.setbounds(0, 0);
            n = 1;
            x(0) = 100*amp::ampf<Precision>::getRandom()-50;
            mincg::mincgcreate<Precision>(n, x, state);
            mincg::mincgsetcgtype<Precision>(state, cgtype);
            while( mincg::mincgiteration<Precision>(state) )
            {
                state.f = -amp::cos<Precision>(state.x(0));
                state.g(0) = amp::sin<Precision>(state.x(0));
            }
            mincg::mincgresults<Precision>(state, x, rep);
            linerror1 = linerror1 || rep.terminationtype<=0 || amp::abs<Precision>(x(0)/amp::pi<Precision>()-amp::round<Precision>(x(0)/amp::pi<Precision>()))>amp::ampf<Precision>("0.001");
            
            //
            // 1D problem #2
            //
            x.setbounds(0, 0);
            n = 1;
            x(0) = 100*amp::ampf<Precision>::getRandom()-50;
            mincg::mincgcreate<Precision>(n, x, state);
            mincg::mincgsetcgtype<Precision>(state, cgtype);
            while( mincg::mincgiteration<Precision>(state) )
            {
                state.f = amp::sqr<Precision>(state.x(0))/(1+amp::sqr<Precision>(state.x(0)));
                state.g(0) = (2*state.x(0)*(1+amp::sqr<Precision>(state.x(0)))-amp::sqr<Precision>(state.x(0))*2*state.x(0))/amp::sqr<Precision>(1+amp::sqr<Precision>(state.x(0)));
            }
            mincg::mincgresults<Precision>(state, x, rep);
            linerror2 = linerror2 || rep.terminationtype<=0 || amp::abs<Precision>(x(0))>amp::ampf<Precision>("0.001");
            
            //
            // Linear equations
            //
            for(n=1; n<=10; n++)
            {
                
                //
                // Prepare task
                //
                a.setbounds(0, n-1, 0, n-1);
                x.setbounds(0, n-1);
                xe.setbounds(0, n-1);
                b.setbounds(0, n-1);
                for(i=0; i<=n-1; i++)
                {
                    xe(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        a(i,j) = 2*amp::ampf<Precision>::getRandom()-1;
                    }
                    a(i,i) = a(i,i)+3*amp::sign<Precision>(a(i,i));
                }
                for(i=0; i<=n-1; i++)
                {
                    v = amp::vdotproduct(a.getrow(i, 0, n-1), xe.getvector(0, n-1));
                    b(i) = v;
                }
                
                //
                // Solve task
                //
                for(i=0; i<=n-1; i++)
                {
                    x(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                mincg::mincgcreate<Precision>(n, x, state);
                mincg::mincgsetcgtype<Precision>(state, cgtype);
                while( mincg::mincgiteration<Precision>(state) )
                {
                    state.f = 0;
                    for(i=0; i<=n-1; i++)
                    {
                        state.g(i) = 0;
                    }
                    for(i=0; i<=n-1; i++)
                    {
                        v = amp::vdotproduct(a.getrow(i, 0, n-1), state.x.getvector(0, n-1));
                        state.f = state.f+amp::sqr<Precision>(v-b(i));
                        for(j=0; j<=n-1; j++)
                        {
                            state.g(j) = state.g(j)+2*(v-b(i))*a(i,j);
                        }
                    }
                }
                mincg::mincgresults<Precision>(state, x, rep);
                eqerror = eqerror || rep.terminationtype<=0;
                for(i=0; i<=n-1; i++)
                {
                    eqerror = eqerror || amp::abs<Precision>(x(i)-xe(i))>amp::ampf<Precision>("0.001");
                }
            }
            
            //
            // Testing convergence properties
            //
            x.setbounds(0, 2);
            n = 3;
            for(i=0; i<=2; i++)
            {
                x(i) = 6*amp::ampf<Precision>::getRandom()-3;
            }
            mincg::mincgcreate<Precision>(n, x, state);
            mincg::mincgsetcond<Precision>(state, amp::ampf<Precision>("0.001"), amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), 0);
            mincg::mincgsetcgtype<Precision>(state, cgtype);
            while( mincg::mincgiteration<Precision>(state) )
            {
                testfunc3<Precision>(state);
            }
            mincg::mincgresults<Precision>(state, x, rep);
            converror = converror || rep.terminationtype!=4;
            for(i=0; i<=2; i++)
            {
                x(i) = 6*amp::ampf<Precision>::getRandom()-3;
            }
            mincg::mincgcreate<Precision>(n, x, state);
            mincg::mincgsetcond<Precision>(state, amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.001"), amp::ampf<Precision>("0.0"), 0);
            mincg::mincgsetcgtype<Precision>(state, cgtype);
            while( mincg::mincgiteration<Precision>(state) )
            {
                testfunc3<Precision>(state);
            }
            mincg::mincgresults<Precision>(state, x, rep);
            converror = converror || rep.terminationtype!=1;
            for(i=0; i<=2; i++)
            {
                x(i) = 6*amp::ampf<Precision>::getRandom()-3;
            }
            mincg::mincgcreate<Precision>(n, x, state);
            mincg::mincgsetcond<Precision>(state, amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.001"), 0);
            mincg::mincgsetcgtype<Precision>(state, cgtype);
            while( mincg::mincgiteration<Precision>(state) )
            {
                testfunc3<Precision>(state);
            }
            mincg::mincgresults<Precision>(state, x, rep);
            converror = converror || rep.terminationtype!=2;
            for(i=0; i<=2; i++)
            {
                x(i) = 2*amp::ampf<Precision>::getRandom()-1;
            }
            mincg::mincgcreate<Precision>(n, x, state);
            mincg::mincgsetcond<Precision>(state, amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), 10);
            mincg::mincgsetcgtype<Precision>(state, cgtype);
            while( mincg::mincgiteration<Precision>(state) )
            {
                testfunc3<Precision>(state);
            }
            mincg::mincgresults<Precision>(state, x, rep);
            converror = converror || !(rep.terminationtype==5 && rep.iterationscount==10 || rep.terminationtype==7);
            
            //
            // Other properties:
            // 1. test reports (F should form monotone sequence)
            // 2. test maximum step
            //
            n = 50;
            x.setlength(n);
            xlast.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                x(i) = 1;
            }
            mincg::mincgcreate<Precision>(n, x, state);
            mincg::mincgsetcond<Precision>(state, amp::ampf<Precision>(0), amp::ampf<Precision>(0), amp::ampf<Precision>(0), 100);
            mincg::mincgsetxrep<Precision>(state, true);
            fprev = amp::ampf<Precision>::getAlgoPascalMaxNumber();
            while( mincg::mincgiteration<Precision>(state) )
            {
                if( state.needfg )
                {
                    state.f = 0;
                    for(i=0; i<=n-1; i++)
                    {
                        state.f = state.f+amp::sqr<Precision>((1+i)*state.x(i));
                        state.g(i) = 2*(1+i)*state.x(i);
                    }
                }
                if( state.xupdated )
                {
                    othererrors = othererrors || state.f>fprev;
                    if( fprev==amp::ampf<Precision>::getAlgoPascalMaxNumber() )
                    {
                        for(i=0; i<=n-1; i++)
                        {
                            othererrors = othererrors || state.x(i)!=x(i);
                        }
                    }
                    fprev = state.f;
                    amp::vmove(xlast.getvector(0, n-1), state.x.getvector(0, n-1));
                }
            }
            mincg::mincgresults<Precision>(state, x, rep);
            for(i=0; i<=n-1; i++)
            {
                othererrors = othererrors || x(i)!=xlast(i);
            }
            n = 1;
            x.setlength(n);
            x(0) = 100;
            stpmax = amp::ampf<Precision>("0.05")+amp::ampf<Precision>("0.05")*amp::ampf<Precision>::getRandom();
            mincg::mincgcreate<Precision>(n, x, state);
            mincg::mincgsetcond<Precision>(state, amp::ampf<Precision>("1.0E-9"), amp::ampf<Precision>(0), amp::ampf<Precision>(0), 0);
            mincg::mincgsetstpmax<Precision>(state, stpmax);
            mincg::mincgsetxrep<Precision>(state, true);
            xprev = x(0);
            while( mincg::mincgiteration<Precision>(state) )
            {
                if( state.needfg )
                {
                    state.f = amp::exp<Precision>(state.x(0))+amp::exp<Precision>(-state.x(0));
                    state.g(0) = amp::exp<Precision>(state.x(0))-amp::exp<Precision>(-state.x(0));
                    othererrors = othererrors || amp::abs<Precision>(state.x(0)-xprev)>(1+amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalEpsilon()))*stpmax;
                }
                if( state.xupdated )
                {
                    othererrors = othererrors || amp::abs<Precision>(state.x(0)-xprev)>(1+amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalEpsilon()))*stpmax;
                    xprev = state.x(0);
                }
            }
        }
        
        //
        // end
        //
        waserrors = referror || eqerror || linerror1 || linerror2 || converror || othererrors;
        if( !silent )
        {
            printf("TESTING CG OPTIMIZATION\n");
            printf("REFERENCE PROBLEM:                        ");
            if( referror )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("LIN-1 PROBLEM:                            ");
            if( linerror1 )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("LIN-2 PROBLEM:                            ");
            if( linerror2 )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("LINEAR EQUATIONS:                         ");
            if( eqerror )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("CONVERGENCE PROPERTIES:                   ");
            if( converror )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            printf("OTHER PROPERTIES:                         ");
            if( othererrors )
            {
                printf("FAILED\n");
            }
            else
            {
                printf("OK\n");
            }
            if( waserrors )
            {
                printf("TEST FAILED\n");
            }
            else
            {
                printf("TEST PASSED\n");
            }
            printf("\n\n");
        }
        result = !waserrors;
        return result;
    }


    /*************************************************************************
    Calculate test function #1
    *************************************************************************/
    template<unsigned int Precision>
    void testfunc1(mincg::mincgstate<Precision>& state)
    {
        if( state.x(0)<100 )
        {
            state.f = amp::sqr<Precision>(amp::exp<Precision>(state.x(0))-2)+amp::sqr<Precision>(state.x(1))+amp::sqr<Precision>(state.x(2)-state.x(0));
            state.g(0) = 2*(amp::exp<Precision>(state.x(0))-2)*amp::exp<Precision>(state.x(0))+2*(state.x(0)-state.x(2));
            state.g(1) = 2*state.x(1);
            state.g(2) = 2*(state.x(2)-state.x(0));
        }
        else
        {
            state.f = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber());
            state.g(0) = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber());
            state.g(1) = 0;
            state.g(2) = 0;
        }
    }


    /*************************************************************************
    Calculate test function #2

    Simple variation of #1, much more nonlinear, which makes unlikely premature
    convergence of algorithm .
    *************************************************************************/
    template<unsigned int Precision>
    void testfunc2(mincg::mincgstate<Precision>& state)
    {
        if( state.x(0)<100 )
        {
            state.f = amp::sqr<Precision>(amp::exp<Precision>(state.x(0))-2)+amp::sqr<Precision>(amp::sqr<Precision>(state.x(1)))+amp::sqr<Precision>(state.x(2)-state.x(0));
            state.g(0) = 2*(amp::exp<Precision>(state.x(0))-2)*amp::exp<Precision>(state.x(0))+2*(state.x(0)-state.x(2));
            state.g(1) = 4*state.x(1)*amp::sqr<Precision>(state.x(1));
            state.g(2) = 2*(state.x(2)-state.x(0));
        }
        else
        {
            state.f = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber());
            state.g(0) = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber());
            state.g(1) = 0;
            state.g(2) = 0;
        }
    }


    /*************************************************************************
    Calculate test function #3

    Simple variation of #1, much more nonlinear, with non-zero value at minimum.
    It achieve two goals:
    * makes unlikely premature convergence of algorithm .
    * solves some issues with EpsF stopping condition which arise when
      F(minimum) is zero

    *************************************************************************/
    template<unsigned int Precision>
    void testfunc3(mincg::mincgstate<Precision>& state)
    {
        amp::ampf<Precision> s;


        s = amp::ampf<Precision>("0.001");
        if( state.x(0)<100 )
        {
            state.f = amp::sqr<Precision>(amp::exp<Precision>(state.x(0))-2)+amp::sqr<Precision>(amp::sqr<Precision>(state.x(1))+s)+amp::sqr<Precision>(state.x(2)-state.x(0));
            state.g(0) = 2*(amp::exp<Precision>(state.x(0))-2)*amp::exp<Precision>(state.x(0))+2*(state.x(0)-state.x(2));
            state.g(1) = 2*(amp::sqr<Precision>(state.x(1))+s)*2*state.x(1);
            state.g(2) = 2*(state.x(2)-state.x(0));
        }
        else
        {
            state.f = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber());
            state.g(0) = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber());
            state.g(1) = 0;
            state.g(2) = 0;
        }
    }


    /*************************************************************************
    Silent unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testmincgunit_test_silent()
    {
        bool result;


        result = testmincg<Precision>(true);
        return result;
    }


    /*************************************************************************
    Unit test
    *************************************************************************/
    template<unsigned int Precision>
    bool testmincgunit_test()
    {
        bool result;


        result = testmincg<Precision>(false);
        return result;
    }
} // namespace

#endif
