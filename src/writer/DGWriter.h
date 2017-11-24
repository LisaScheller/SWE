//
// Created by lisa on 23.11.17.
//

#ifndef SWE1D_DGWRITER_H
#define SWE1D_DGWRITER_H
#include "../types.h"

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>

namespace writer
{

/**
 * A writer class that generates vtk files
 */
    class DGWriter
    {
    private:
        // base name of the vtp collectiond and vtk files
        std::string m_basename;

        // cell size
        T m_cellSize;

        // current time step
        unsigned int m_timeStep;

        // vtp stream
        std::ofstream *m_vtpFile;


    public:
        // constructor
        DGWriter( const std::string& basename = "DG", const T cellSize = 1)
                : m_basename(basename),
                  m_cellSize(cellSize),
                  m_timeStep(0)
        {
            // initialize vtp stream
            std::ostringstream l_vtpFileName;
            l_vtpFileName << m_basename << ".vtp";

            m_vtpFile = new std::ofstream( l_vtpFileName.str().c_str() );

            // write vtp header
            *m_vtpFile
                    << "<?xml version=\"1.0\"?>" << std::endl
                    << "<VTKFile type=\"Collection\" version=\"0.1\">" << std::endl
                    << "<Collection>" << std::endl;

        }

        // destructor (free memory)
        ~DGWriter() {
            // close vtp file
            *m_vtpFile
                    << "</Collection>" << std::endl
                    << "</VTKFile>" << std::endl;
            delete m_vtpFile;
        }

        /**
         * Writes all values to vtk file
         *
         * @param size Number of cells (without boundary values)
         */
        void write(const T time, const vecq q, unsigned int size)
        {
            // generate vtk file name
            std::string l_fileName = generateFileName();

            // add current time to vtp collection
            *m_vtpFile << "<DataSet timestep=\""
                       << time
                       << "0\" group=\"\" part=\"0\" file=\""
                       << l_fileName
                       << "\"/> " << std::endl;

            // write vtk file
            std::ofstream vtkFile(l_fileName.c_str());
            assert(vtkFile.good());

            // vtk xml header
            vtkFile << "<VTKFile type=\"UnstructuredGrid\">" << std::endl;



            vtkFile << "<UnstructuredGrid>" << std::endl;
            vtkFile << "<Piece NumberOfPoints=\""<< size*2<< "\" NumberOfCells=\""<< size <<"\">" << std::endl;
            vtkFile << "<Points>" << std::endl;
            vtkFile << "<DataArray Name=\"xyz\" NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">" << std::endl;

            //  points
            vtkFile << m_cellSize * 0 << " " << 0.0 << " " << 0.0 <<std::endl;
            for (int i=0; i < size; i++) {
                vtkFile << m_cellSize * i << " " << 0.0 << " " << 0.0 <<std::endl;
                vtkFile << m_cellSize * i << " " << 0.0 << " " << 0.0 <<std::endl;
            }
            vtkFile << m_cellSize * size << " " << 0.0 << " " << 0.0 <<std::endl;

            vtkFile << "</DataArray>" << std::endl;
            vtkFile << "</Points>" << std::endl;

            vtkFile << "<Cells>" << std::endl;
            vtkFile << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
            for (unsigned int i = 0; i< size+1; i++){
                vtkFile << 2*i << " " << (2*i)+1 << std::endl;
            }
            vtkFile << "</DataArray>" << std::endl;
            vtkFile << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
            for (unsigned int i = 0; i< size+1; i++){
                vtkFile<< 2+(2*i) << std::endl;
            }
            vtkFile << "</DataArray>" << std::endl;
            vtkFile << "<DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">" << std::endl;
            for (unsigned int i = 0; i< size+1; i++){
                vtkFile << 3 << std::endl;
            }
            vtkFile << "</DataArray>" << std::endl;
            vtkFile << "</Cells>" << std::endl;

            vtkFile << "<PointData>" << std::endl;

            // water surface height
            vtkFile << "<DataArray Name=\"h\" NumberOfComponents=\"1\" type=\"Float64\"  format=\"ascii\">" << std::endl;

            for (int i=0; i < size+1; i++) {
                vtkFile << q.at(i).h.u0 << std::endl;
                vtkFile << q.at(i).h.u1 << std::endl;
            }
            vtkFile << "</DataArray>" << std::endl;

            // momentum
            vtkFile << "<DataArray Name=\"hu\" NumberOfComponents=\"1\" type=\"Float64\" format=\"ascii\">" << std::endl;

            for (int i=0; i < size+1; i++){
                vtkFile << q.at(i).hu.u0 << std::endl;
                vtkFile << q.at(i).hu.u1 << std::endl;
            }

            vtkFile << "</DataArray>" << std::endl;

            vtkFile << "</PointData>" << std::endl;

            vtkFile << "<CellData>" << std::endl;
            vtkFile << "</CellData>" << std::endl;

            // bathymetry
            //vtkFile << "<DataArray Name=\"B\" type=\"Float32\" format=\"ascii\">" << std::endl;
            //for (int i=1; i<size+1; i++)
            //		vtkFile << b[i] << std::endl;
            //vtkFile << "</DataArray>" << std::endl;

            vtkFile << "</Piece>" << std::endl;

            vtkFile << "</UnstructuredGrid>" << std::endl
                    << "</VTKFile>" << std::endl;

            // increment time step
            m_timeStep++;
        }

    private:
        /**
         * @return The generated filename containing the time step and the real name
         */
        std::string generateFileName()
        {
            std::ostringstream name;
            name << m_basename << '_' << m_timeStep << ".vtu";

            return name.str();
        }
    };

}
#endif //SWE1D_DGWRITER_H
