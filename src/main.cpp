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

#include <cstring>

T meanOfN(int n, T *q, int j){
    T sum = 0.0;
    for (int i = 0; i<n; i++){
        sum += q[(n*j)-i];
    }
    return sum/n;
}



int main(int argc, char** argv)
{
	// Parse command line parameters
	tools::Args args(argc, argv);
    int factor = 4;
    int numberOfIntervals = (factor*args.size());


	// Scenario (Choose between Dambreak and gaussian)
	scenarios::DamBreak scenario(numberOfIntervals);
	//scenarios::Gaussian scenario(args.size());

	// Allocate memory
	// Water height

	T *h = new T[numberOfIntervals+2];
    T *ah = new T[numberOfIntervals+2];
	// Momentum
	T *hu = new T[numberOfIntervals+2];
    T *ahu = new T[numberOfIntervals+2];

    T *h_num = new T[args.size()+2];
    T *hu_num = new T[args.size()+2];
    T *ah_num = new T[args.size()+2];
    T *ahu_num = new T[args.size()+2];

	// Initialize water height and momentum
	for (unsigned int i = 0; i < numberOfIntervals+2; i++) {
        //int j = (int) i / 2;
        h[i] = scenario.getHeight(i);
    }
	memset(hu, 0, sizeof(T)*(numberOfIntervals+2));

    for (unsigned int i = 0; i < numberOfIntervals+2; i++) {
        //int j = (int) i/2;
        ah[i] = scenario.getHeight(i);
    }
    memset(ahu, 0, sizeof(T)*(numberOfIntervals+2));

    for(unsigned int i = 0; i<args.size()+2; i++){
        if (i==0){
            h_num[i]=h[i];
        }
        else if (i==args.size()+1){
            h_num[i] = h[numberOfIntervals+1];
        }
        //else h_num[i]=(h[2*i]+h[(2*i)-1])/2;
        //else h_num[i]=meanOfN(factor, h_num, i);
        else {
            T sum = 0.0;
            for (int j = 0; j<factor; j++){
                sum += h[(factor*i)-j];
            }
            h_num[i]=sum/factor;
        }
    }
    memset(hu_num, 0, sizeof(T)*(args.size()+2));
    for(unsigned int i = 0; i<args.size()+2; i++){
        if (i==0){
            ah_num[i]=h[i];
        }
        else if (i==args.size()+1){
            ah_num[i] = h[numberOfIntervals+1];
        }
        //else ah_num[i]=(h[2*i]+h[(2*i)-1])/2;
        else ah_num[i]=meanOfN(factor, ah, i);
    }
    memset(ahu_num, 0, sizeof(T)*(args.size()+2));

	// Create a writer that is responsible printing out values
	//writer::ConsoleWriter writer;
	writer::VtkWriter writer("swe1d", scenario.getCellSize());
    writer::VtkWriter analyticalWriter("analytical", scenario.getCellSize());
    writer::VtkWriter nodalAdvectionWriter("nodaladvection", 1);

	// Helper class computing the wave propagation
	WavePropagation wavePropagation(h, hu, ah, ahu, numberOfIntervals, scenario.getCellSize());
    //Helper class for compution the solution of advection equation using nodal DG
    NodalAdvection nodalAdvection();

	// Write initial data
	tools::Logger::logger.info("Initial data");

	// Current time of simulation
    T t = 0;




	writer.write(t, h_num, hu_num, args.size());
    analyticalWriter.write(t, ah_num, ahu_num, args.size());

    T error1 = 0.0;
    T error2 = 0.0;
    T error3 = 0.0;

	for (unsigned int i = 0; i < args.timeSteps(); i++) {
		// Do one time step
		tools::Logger::logger << "Computing timestep " << i
				<< " at time " << t << std::endl;

		// Update boundaries
		wavePropagation.setOutflowBoundaryConditions();
        nodalAdvection.setBoundaryConditions();

		// Compute numerical flux on each edge
		//T maxTimeStep = wavePropagation.computeNumericalFluxes();
        T maxTimeStep = wavePropagation.computeLaxFriedrichsFlux(t);
        T advTimeStep = nodalAdvection.compute(t);

		// Update unknowns from net updates (Choose between unstable and Lax Friedrichs)
		//wavePropagation.updateUnknowns(maxTimeStep);
        //wavePropagation.updateUnknownsUnstable(maxTimeStep);
        wavePropagation.updateUnknownsLaxFriedrichs(maxTimeStep);
        nodalAdvection.updateUnknowns(advTimeStep);



		// Update time
		//t += maxTimeStep;
        t += advTimeStep;

        h_num[0]=h[0];
        hu_num[0]=hu[0];
        h_num[args.size()+1]=h[numberOfIntervals+1];
        hu_num[args.size()+1]=hu[numberOfIntervals+1];

        ah_num[0]=ah[0];
        ahu_num[0]=ahu[0];
        ah_num[args.size()+1]=ah[numberOfIntervals+1];
        ahu_num[args.size()+1]=ahu[numberOfIntervals+1];
        for(int j=1; j<args.size()+1; j++) {
            /*h_num[j] = meanOfN(factor,h_num,j);*/
            T sum = 0.0;
            for (int k = 0; k<factor; k++){
                sum += h[(factor*j)-k];
            }
            h_num[j]=sum/factor;
            hu_num[j] = meanOfN(factor, hu, j) ;
            ah_num[j] = meanOfN(factor, ah, j);
            ahu_num[j] = meanOfN(factor, ahu,j);
            /*h_num[j] = (h[2 * j] + h[(2 * j) - 1]) / 2;
            hu_num[j] = (hu[2 * j] + hu[(2 * j) - 1] / 2);
            ah_num[j] = (ah[2 * j] + ah[(2 * j) - 1]) / 2;
            ahu_num[j] = (ahu[2 * j] + ahu[(2 * j) - 1] / 2);*/

        }
        //Compute difference between exact solution and numerical method at t=2,5
        if(i==57){
            error1 = wavePropagation.computeError();
            /*T res = 0.0;
            for (int i = 0; i<args.size()+1; i++){
                res += sqrtf((std::abs(h_num[i]-ah_num[i])*std::abs(h_num[i]-ah_num[i]))+(std::abs(hu_num[i]-ahu_num[i])*std::abs(hu_num[i]-ahu_num[i])));
            }
            error1 = res;*/

        }
        if(i==10){
            error2 = wavePropagation.computeError();
            /*T res = 0.0;
            for (int i = 0; i<args.size()+1; i++){
                res += sqrtf((std::abs(h_num[i]-ah_num[i])*std::abs(h_num[i]-ah_num[i]))+(std::abs(hu_num[i]-ahu_num[i])*std::abs(hu_num[i]-ahu_num[i])));
            }
            error2 = res;*/

        }
        if(i==21){
            error3 = wavePropagation.computeError();
            /*T res = 0.0;
            for (int i = 0; i<args.size()+1; i++){
                res += sqrtf((std::abs(h_num[i]-ah_num[i])*std::abs(h_num[i]-ah_num[i]))+(std::abs(hu_num[i]-ahu_num[i])*std::abs(hu_num[i]-ahu_num[i])));
            }
            error3 = res;*/

        }


		// Write new values
		//writer.write(t, h, hu, numberOfIntervals/2);
        writer.write(t, h_num, hu_num, numberOfIntervals / factor);
        //analyticalWriter.write(t, ah, ahu, numberOfIntervals/2);
        analyticalWriter.write(t, ah_num, ahu_num, numberOfIntervals/factor);
        nodalAdvectionWriter.write();


	}


	// Free allocated memory
	delete [] h;
	delete [] hu;
    delete [] ah;
    delete [] ahu;
    delete [] h_num;
    delete [] hu_num;

    std::cout << "Error at time approx. 0.5s is " << error1
              <<  std::endl;
    std::cout << "Error at time approx. 1.0s is " << error2
              <<  std::endl;
    std::cout << "Error at time approx. 2.5s is " << error3
              <<  std::endl;
	return 0;

}
