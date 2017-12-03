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

#include "WavePropagation.h"
#include "LaxFriedrichsSolver.h"
#include <algorithm>
#include <cmath>
 /*T WavePropagation::computeNumericalFluxes()
 {
 	float maxWaveSpeed = 0.f;

 	// Loop over all edges
 	for (unsigned int i = 1; i < m_size+2; i++) {
 		T maxEdgeSpeed;

 		// Compute net updates
 		m_solver.computeNetUpdates( m_h[i-1], m_h[i],
 				m_hu[i-1], m_hu[i],
				0.f, 0.f,	// Bathymetry
				m_hNetUpdatesLeft[i-1], m_hNetUpdatesRight[i-1],
 				m_huNetUpdatesLeft[i-1], m_huNetUpdatesRight[i-1],
 				maxEdgeSpeed );

 		// Update maxWaveSpeed
 		if (maxEdgeSpeed > maxWaveSpeed)
 			maxWaveSpeed = maxEdgeSpeed;
 	}

 	// Compute CFL condition
	T maxTimeStep = m_cellSize/maxWaveSpeed * .4f;

	return maxTimeStep;
 }*/

T WavePropagation::computeLaxFriedrichsFlux(T t) {
  return solver.computeLocalLaxFriedrichsFlux(t);
     //return solver.computeLaxFriedrichsFlux();
     //return solver.computeLaxFriedrichsFlux2(t);

}

T WavePropagation::computeUnstableFlux(T t) {
    return solver.computeUnstableFlux(t);
}



vecu WavePropagation::updateUnknownsUnstable(T dt){
    return solver.updateUnknownsUnstable(dt);

}

vecu WavePropagation::updateUnknownsLaxFriedrichs(T dt){
   return solver.updateUnknownsLocalLaxFriedrichs(dt);
    //solver.updateUnknownsLaxFriedrichs(dt);
    //solver.updateUnknownsLaxFriedrichsDirect(dt);
    //return solver.updateUnknownsLaxFriedrichs2(dt);
    solver.solveAnalytically(dt);

}

vecu WavePropagation::updateAnalyticalSolution(T dt){
    return solver.getAnalyticalSolution(dt);
}

void WavePropagation::updateUnknowns(T dt)
{
	// Loop over all inner cells
	for (unsigned int i = 1; i < m_size+1; i++) {
        m_h.at(i) -=  dt/m_cellSize * (m_hNetUpdatesRight.at(i-1) + m_hNetUpdatesLeft.at(i));
        m_hu.at(i) -= dt/m_cellSize * (m_huNetUpdatesRight.at(i-1)+ m_huNetUpdatesLeft.at(i));
	}
}

void WavePropagation::setOutflowBoundaryConditions()
{
	m_h.at(0) = m_h.at(1); m_h.at(m_size+1) = m_h.at(m_size);
	m_hu.at(0) = m_hu.at(1); m_hu.at(m_size+1) = m_hu.at(m_size);
}

T WavePropagation::computeError() {
    T res = 0.0;
    for (int i = 0; i<m_size+1; i++){
        //res += sqrtf((std::abs(m_cellSize*(m_h[i]-solver.a_h[i]))*std::abs(m_cellSize*(m_h[i]-solver.a_h[i])))+(std::abs(m_cellSize*(m_hu[i]-solver.a_hu[i]))*std::abs(m_cellSize*(m_hu[i]-solver.a_hu[i]))));
        res += std::abs((m_cellSize*(m_h.at(i)-solver.a_h.at(i)))*std::abs(m_cellSize*(m_h.at(i)-solver.a_h.at(i))))+(std::abs(m_cellSize*(m_hu.at(i)-solver.a_hu.at(i)))*std::abs(m_cellSize*(m_hu.at(i)-solver.a_hu.at(i))));
        //res += std::abs(m_cellSize*(m_h[i]-solver.a_h[i]))+std::abs(m_cellSize*(m_hu[i]-solver.a_hu[i]));
    }
    return sqrt(res);
}

vect WavePropagation::getExactSolutionH(T t){
    vect res(m_size+2);
    vecu sol = solver.getAnalyticalSolution(t);
    for(int i = 0; i<m_size+2; i++){
        res.at(i) = sol.at(i).u0;
    }
    return res;
}

vect WavePropagation::getExactSolutionHu(T t){
    vect res(m_size+2);
    vecu sol = solver.getAnalyticalSolution(t);
    for(int i = 0; i<m_size+2; i++){
        res.at(i) = sol.at(i).u1;
    }
    return res;
}


