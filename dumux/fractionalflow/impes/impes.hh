// $Id$
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUNE_IMPES_HH
#define DUNE_IMPES_HH

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include "dumux/fractionalflow/fractionalflow.hh"

/**
 * @file
 * @brief  IMPES scheme
 * @author Bernd Flemisch, Markus Wolff
 */

namespace Dune
{
/**
 * \ingroup fracflow
 * @brief IMplicit Pressure Explicit Saturation (IMPES) scheme for the solution of
 * coupled diffusion/transport problems
 */

template<class GridView, class Diffusion, class Transport, class VC> class IMPES: public FractionalFlow<
        GridView, Diffusion, Transport, VC>
{
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
typedef    typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef Dune::FractionalFlow<GridView, Diffusion, Transport, VC> FractionalFlow;
    typedef typename FractionalFlow::RepresentationType PressType;
    typedef typename FractionalFlow::Scalar Scalar;

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

public:
    typedef typename Transport::RepresentationType RepresentationType;

    virtual void initial()
    {
        Scalar t = 0;
        //initial saturations
        this->transport.initialTransport();
        //call function with true to get a first initialisation of the pressure field
        this->diffusion.pressure(true,t);
        this->diffusion.calculateVelocity(t);

        return;
    }

    //calculate saturation defect
    virtual int update(const Scalar t, Scalar& dt, RepresentationType& updateVec,
            Scalar cFLFactor = 1)
    {
        PressType pressOldIter(this->transport.problem().variables().pressure());
        PressType pressHelp(this->transport.problem().variables().gridSizeDiffusion());
        int satSize = this->transport.problem().variables().gridSizeTransport();
        RepresentationType saturation(this->transport.problem().variables().saturation());
        RepresentationType satOldIter(this->transport.problem().variables().saturation());
        RepresentationType satHelp(satSize);
        RepresentationType satDiff(satSize);
        RepresentationType updateOldIter(satSize);
        RepresentationType updateHelp(satSize);
        RepresentationType updateDiff(satSize);

        //update constitutive functions
        this->transport.updateMaterialLaws();

        bool converg = false;
        int iter = 0;
        int iterTot = 0;
        updateOldIter = 0;
        while (!converg)
        {
            iter++;
            iterTot++;

            // update pressure: give false as the pressure field is already initialised
            this->diffusion.pressure(false , t);
            this->diffusion.calculateVelocity(t);

            //calculate saturation defect
            this->transport.update(t, dt, updateVec,cFLFactor, true);

            if (iterFlag)
            { // only needed if iteration has to be done
                this->transport.problem().variables().pressure() *= omega;
                pressHelp = pressOldIter;
                pressHelp *= (1-omega);
                this->transport.problem().variables().pressure() += pressHelp;
                pressOldIter = this->transport.problem().variables().pressure();

                updateHelp = updateVec;
                saturation = this->transport.problem().variables().saturation();
                saturation += (updateHelp *= (dt*cFLFactor));
                saturation *= omega;
                satHelp = satOldIter;
                satHelp *= (1-omega);
                saturation += satHelp;
                updateDiff = updateVec;
                updateDiff -= updateOldIter;
                satOldIter = saturation;
                updateOldIter = updateVec;
            }
            // break criteria for iteration loop
            if (iterFlag == 2 && dt * updateDiff.two_norm() / saturation.two_norm() <= maxDefect )
            {
                converg = true;
            }
            else if (iterFlag == 2 && (saturation.infinity_norm() > 1 || saturation.two_norm()> 1))
            {
                converg = false;
            }
            else if (iterFlag == 2 && iter> nIter )
            {
                std::cout << "Nonlinear loop in IMPES.update exceeded nIter = "
                << nIter << " iterations."<< std::endl;
                return 1;
            }
            else if (iterFlag == 1 && iter> nIter )
            {
                converg = true;
            }
            else if (iterFlag==0)
            {
                converg = true;
            }
        }
        // outputs
        if (iterFlag==2)
        std::cout << "Iteration steps: "<< iterTot << std::endl;
        std::cout.setf(std::ios::scientific, std::ios::floatfield);

        return 0;
    }

    virtual void vtkout(const char* name, int k) const
    {
        this->transport.problem().variables().vtkout(name, k);
        return;
    }

    //! Construct an IMPES object.
    IMPES(Diffusion& diffusion, Transport& transport, int flag = 0, int nIt = 2,
            Scalar maxDef = 1e-5, Scalar om = 1) :
    FractionalFlow(diffusion, transport),
    iterFlag(flag), nIter(nIt), maxDefect(maxDef), omega(om)
    {
    }

protected:
    const int iterFlag;
    const int nIter;
    const Scalar maxDefect;
    const Scalar omega;
};
}
#endif
