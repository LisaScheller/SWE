/**
 * @file
 *  This file is part of SWE1D
 *
 *  SWE1D is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SWE1D is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with SWE1D.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Diese Datei ist Teil von SWE1D.
 *
 *  SWE1D ist Freie Software: Sie koennen es unter den Bedingungen
 *  der GNU General Public License, wie von der Free Software Foundation,
 *  Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
 *  veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
 *
 *  SWE1D wird in der Hoffnung, dass es nuetzlich sein wird, aber
 *  OHNE JEDE GEWAEHELEISTUNG, bereitgestellt; sogar ohne die implizite
 *  Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FUER EINEN BESTIMMTEN
 *  ZWECK. Siehe die GNU General Public License fuer weitere Details.
 *
 *  Sie sollten eine Kopie der GNU General Public License zusammen mit diesem
 *  Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.
 * 
 * @copyright 2013 Technische Universitaet Muenchen
 * @author Sebastian Rettenberger <rettenbs@in.tum.de>
 */

#include "types.h"
#include "WavePropagation.h"
#include "scenarios/dambreak.h"
//#include "writer/ConsoleWriter.h"
#include "writer/VtkWriter.h"
#include "tools/args.h"
#include "scenarios/gaussian.h"
#include "NodalAdvection.h"
#include "scenarios/AdvectionGalerkin.h"
#include "writer/AdvectionWriter.h"
#include "writer/GalerkinWriter.h"
#include "writer/DGWriter.h"
#include "scenarios/SWEGalerkin.h"
#include "GalerkinSWE.h"
#include "AdvectionFiniteVolumes.h"

#include <cstring>


int main(int argc, char** argv)
{
    //Choose which scenarios or solvers you want to use
    bool dambreak = false;
    bool gauss = false;
    bool Adv = false;
    bool AdvDG = true;
    bool SWE = false;
    bool SWEDG = false;
    bool SWEDGGauss = false;

    //Choose between normal mode and error computing (no vtk output while error computing)
    bool errorComputing = true;

	// Parse command line parameters
	tools::Args args(argc, argv);
    int numberOfIntervals = args.size();

	// Initialize scenarios

    scenarios::DamBreak scDambreak(numberOfIntervals);

    scenarios::Gaussian scGauss(numberOfIntervals);

    scenarios::SWEGalerkin scSWEDG(numberOfIntervals);

    scenarios::AdvectionGalerkin scAdvDG(numberOfIntervals);

	// Water height and momentum

    //SWE Finite Volume
    vect hFV(numberOfIntervals+2, 0.0);
    vect huFV(numberOfIntervals+2, 0.0);
    vect ahFV(numberOfIntervals+2, 0.0);
    vect ahuFV(numberOfIntervals+2, 0.0);

    //Advection finite volumes
    vect u(numberOfIntervals+2, 0.0);
    vect uu(numberOfIntervals+2, 0.0);
    vect u0(numberOfIntervals+2, 0.0);

    //Advection equation DG
    vecu h(numberOfIntervals+2);
    vecu hu(numberOfIntervals+2);
    vecq qAdv(numberOfIntervals+2);
    vecu hAdvAnalytic(numberOfIntervals+2);
    vecq qAdvAnalytic(numberOfIntervals+2);
    //Wavespeed advection DG
    T a = -1.25;

    //SWE DG
    vecq q(numberOfIntervals+2);


	// Initialize water height and momentum
    if(SWE && dambreak) {
        for (unsigned int i = 0; i < numberOfIntervals + 2; i++) {
            hFV.at(i) = scDambreak.getHeight(i);
            ahFV.at(i) = scDambreak.getHeight(i);
            //huFV.at(i) = scDambreak.getU(i);
        }
    }
    else if (SWE && gauss){
        for (unsigned int i = 0; i < numberOfIntervals + 2; i++) {
            hFV.at(i) = scGauss.getHeight(i);
            ahFV.at(i) = scGauss.getHeight(i);
        }
    }
    if (Adv){
        for (unsigned int i = 0; i < numberOfIntervals + 2; i++) {
            u.at(i) = scGauss.getHeight(i);
        }
    }
    if (AdvDG){
        for (unsigned int i = 0; i<numberOfIntervals+2; i++){
            h.at(i) = scAdvDG.getHeight(i);
            hAdvAnalytic.at(i) = scAdvDG.getHeight(i);
            qAdv.at(i).h = h.at(i);
            qAdv.at(i).hu.u0 = 0.0;
            qAdv.at(i).hu.u1 = 0.0;
            qAdvAnalytic.at(i).h = hAdvAnalytic.at(i);
            qAdvAnalytic.at(i).hu.u0 = 0.0;
            qAdvAnalytic.at(i).hu.u1 = 0.0;
        }
    }
    if (SWEDG && SWEDGGauss) {
        for (unsigned int i = 0; i < numberOfIntervals + 2; i++) {
            q.at(i) = scSWEDG.getHeightGauss(i);
        }
    }
    else if (SWEDG){
        for (unsigned int i = 0; i < numberOfIntervals + 2; i++) {
            q.at(i) = scSWEDG.getHeight(i);
        }
    }

    // Create a writer that is responsible for printing out values
    writer::VtkWriter writerG("SWE1DGauss", scGauss.getCellSize());
    writer::VtkWriter writer("swe1d", scDambreak.getCellSize());
    writer::VtkWriter analyticWriter("analytic", scDambreak.getCellSize());
    writer::VtkWriter advWriter("Adv", scGauss.getCellSize());
    writer::DGWriter dgWriter("AdvDG", scAdvDG.getCellSize());
    writer::DGWriter analyticalAdvWriter("AdvDGAnalytic", scAdvDG.getCellSize());
    writer::DGWriter writerSWEG("SWEDGGauss", scSWEDG.getCellSize());
    writer::DGWriter writerSWEDG("SWEDG", scSWEDG.getCellSize());

    // Helper classes for computing the wave propagation and computing the solutions
	WavePropagation wavePropagation(hFV, huFV, ahFV, ahuFV, numberOfIntervals, scDambreak.getCellSize());
    AdvectionFiniteVolumes advection(a, u, numberOfIntervals, scGauss.getCellSize());
    NodalAdvection nodalAdvection(a,h,numberOfIntervals,scAdvDG.getCellSize());
    GalerkinSWE galerkin(q,numberOfIntervals,scSWEDG.getCellSize());

	// Write initial data
	tools::Logger::logger.info("Initial data");

	// Current time of simulation
    T t = 0;

    if (!errorComputing) {
        if (SWE && dambreak) {
            writer.write(t, hFV, huFV, numberOfIntervals);
            analyticWriter.write(t, ahFV, ahuFV, numberOfIntervals);

        } else if (SWE && gauss) {
            writerG.write(t, hFV, huFV, numberOfIntervals);
        }
        if (Adv){
            advWriter.write(t, u, uu, numberOfIntervals);
        }
        if (AdvDG) {
            dgWriter.write(t, qAdv, numberOfIntervals);
            if (errorComputing) {
                analyticalAdvWriter.write(t, qAdvAnalytic, numberOfIntervals);
            }
        }
        if (SWEDG && SWEDGGauss) {
            writerSWEG.write(t, q, numberOfIntervals);
        } else if (SWEDG) {
            writerSWEDG.write(t, q, numberOfIntervals);
        }
    }

    //Error variables needed for error computing
    T errorH = 0.0;
    T errorHu = 0.0;
    T error = 0.0;

    //Variables to compute the amount of needed steps until Advection equation reaches its starting point again
    //Time for one full round
    T tRound = 100/std::abs(a);
    //Save starting conditions
    vecu h0 = h;
    u0 = u;
    //Initalize variables for calculation of number of timesteps
    T timestepsize = 0.0;
    int iRound = 0;
    int maxTimeSteps = args.timeSteps();

    //Variable for computing the number of timesteps until a certain time is reached (for error computing)
    T errorTime = 0.5;
    T timestep = 0.0;
    T neededTime = 0.0;
    int iNeeded = 0;


	for (unsigned int i = 0; i < maxTimeSteps+1; i++) {
		// Do one time step
		tools::Logger::logger << "Computing timestep " << i
				<< " at time " << t << std::endl;

        if (SWE){
            // Update boundaries
            wavePropagation.setOutflowBoundaryConditions();

            // Compute numerical flux on each edge
            T maxTimeStep = wavePropagation.computeLaxFriedrichsFlux(t);
            //T maxTimeStep = wavePropagation.computeUnstableFlux(t);

            if (t <= errorTime && errorComputing){
                maxTimeSteps  = i + 10;
            }
            else if (t > errorTime && errorComputing){
                maxTimeSteps = i;
            }
            // Update unknowns from net updates (Choose between Lax-Friedrichs and Local Lax Friedrichs in WavePropagation.cpp)
            vecu res = wavePropagation.updateUnknownsLaxFriedrichs(maxTimeStep);
            //vecu res = wavePropagation.updateUnknownsUnstable(maxTimeStep);

            for (unsigned int j = 0; j<numberOfIntervals+2; j++){
                hFV.at(j) = res.at(j).u0;
                huFV.at(j) = res.at(j).u1;
            }

            //Save values of analytic solution if error shall be computed
            vecu resA = wavePropagation.updateAnalyticalSolution(maxTimeStep);
            for (unsigned int j = 0; j<numberOfIntervals+2; j++){
                ahFV.at(j) = resA.at(j).u0;
                ahuFV.at(j) = resA.at(j).u1;
                }


            // Update time
            t += maxTimeStep;

            if (errorComputing){
                //Compute difference between exact solution and numerical method at errorTime
                if (t>=errorTime && std::abs(t-errorTime)<maxTimeStep) {
                    ahFV = wavePropagation.getExactSolutionH(t);
                    ahuFV = wavePropagation.getExactSolutionHu(t);
                    T deltaH = 0.0;
                    T deltaHu = 0.0;
                    T errorAtI = 0.0;
                    vect errH(numberOfIntervals+2, 0.0);
                    vect  errHu(numberOfIntervals+2, 0.0);
                    T eH  = 0.0;
                    T eHu = 0.0;
                    for (unsigned int i = 1; i < numberOfIntervals + 1; i++) {
                        deltaH = std::abs(hFV.at(i) - ahFV.at(i));
                        deltaHu = std::abs(huFV.at(i) - ahuFV.at(i));
                        errorAtI = deltaH+deltaHu;
                        errH.at(i) = deltaH;
                        errHu.at(i) = deltaHu;
                    }
                    //Variante 1 Betrag der Vektoren als ganzes
                    for (int j = 1; j< numberOfIntervals+1; j++){
                        eH += errH.at(j);
                        eHu += errHu.at(j);
                    }
                    error = scDambreak.getCellSize()*(eH+eHu);
                }
            }

            // Write new values
            if (!errorComputing) {
                if (gauss) {
                    writerG.write(t, hFV, huFV, numberOfIntervals);
                } else if (dambreak) {
                    writer.write(t, hFV, huFV, numberOfIntervals);
                    analyticWriter.write(t, ahFV, ahuFV, numberOfIntervals);

                }
            }
        }

        if (Adv){
            if (iRound > 0 && i == iRound+1){
                error = advection.computeError(u0);
            }

            //Update boundaries
            advection.setBoundaryConditions();

            //Compute numerical flux on each edge
            T advFVTimeStep = advection.computeLocalLaxFriedrichsFlux(t);

            if (i == 0 && errorComputing){
                timestepsize = advFVTimeStep;
                iRound = ((int)tRound/timestepsize)+1;
                maxTimeSteps = iRound + 10;
            }

            advection.updateUnknownsLocalLaxFriedrichs(advFVTimeStep);

            t+= advFVTimeStep;

            vect tmp = advection.setU();
            for (unsigned int j = 0; j<numberOfIntervals+2; j++){
                u.at(j) = tmp.at(j);
            }

            //Write values to vtk
            if(!errorComputing){
                advWriter.write(t, u, uu, numberOfIntervals);
            }

        }
        if (AdvDG){
            if (iRound > 0 && i == iRound+1){
                error = nodalAdvection.computeError(h0);
            }

            // Update boundaries
            nodalAdvection.setBoundaryConditions();

            // Compute numerical flux on each edge
            T advTimeStep = nodalAdvection.computeLocalLaxFriedrichsFluxes(t);

            if (i == 0 && errorComputing){
                timestepsize = advTimeStep;
                iRound = ((int)tRound/timestepsize)+1;
                maxTimeSteps = iRound + 10;
            }

            //Compute time derivative of u
            nodalAdvection.computeTimeDerivative();

            //Compute the euler integration step
            nodalAdvection.computeEulerStep(advTimeStep);

            // Update time
            t += advTimeStep;

            // Write new values
            h = nodalAdvection.setH();
            for (unsigned int k = 0; k<numberOfIntervals+2; k++){
                qAdv.at(k).h.u0 = h.at(k).u0;
                qAdv.at(k).h.u1 = h.at(k).u1;
                qAdv.at(k).hu.u0 = 0.0;
                qAdv.at(k).hu.u1 = 0.0;
            }
            for (unsigned int p = 0; p<numberOfIntervals+2; p++){
                hAdvAnalytic.at(p).u0 = nodalAdvection.getExactSolution(t, numberOfIntervals).at(p);
                hAdvAnalytic.at(p).u1 = nodalAdvection.getExactSolution(t, numberOfIntervals).at(p);
                qAdvAnalytic.at(p).h.u0 = hAdvAnalytic.at(p).u0;
                qAdvAnalytic.at(p).h.u1 = hAdvAnalytic.at(p).u1;
            }

            if (!errorComputing) {
                dgWriter.write(t, qAdv, numberOfIntervals);
                analyticalAdvWriter.write(t, qAdvAnalytic, numberOfIntervals);
            }
        }

        if (SWEDG){
            // Update boundaries
            galerkin.setBoundaryConditions();

            // Compute numerical flux on each edge
            T galTimeStep = galerkin.computeLocalLaxFriedrichsFluxes(t);

            //Compute time derivative of h and hu
            galerkin.computeTimeDerivative();

            //Compute Euler integration step
            galerkin.computeEulerStep(galTimeStep);

            // Update time
            t += galTimeStep;

            // Write new values
            q = galerkin.setQ();
            if (!errorComputing) {
                if (SWEDGGauss) {
                    writerSWEG.write(t, q, numberOfIntervals);
                } else{
                    writerSWEDG.write(t, q, numberOfIntervals);
                }
            }
        }
    }
    if (errorComputing && SWE && dambreak){
        std::cout << "Error at time approx. 0.5s is " << error <<std::endl;
    }
    if (errorComputing && AdvDG){
        std::cout << "Error after one circle is " << error << std::endl;
    }
    if (errorComputing && Adv){
        std::cout << "Error after one circle is " << error << std::endl;
    }

	return 0;

}
