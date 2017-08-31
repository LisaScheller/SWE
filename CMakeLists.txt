cmake_minimum_required(VERSION 3.8)
project(SWE1D)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        src/scenarios/dambreak.h
        src/tools/args.cpp
        src/tools/args.h
        src/tools/logger.cpp
        src/tools/logger.h
        src/writer/ConsoleWriter.h
        src/writer/VtkWriter.h
        src/main.cpp
        src/types.h
        src/WavePropagation.cpp
        src/WavePropagation.h
        submodules/solvers/src/solver/AugRie.hpp
        submodules/solvers/src/solver/AugRieGeoClaw.hpp
        submodules/solvers/src/solver/FWave.hpp
        submodules/solvers/src/solver/FWaveCuda.h
        submodules/solvers/src/solver/FWaveVec.hpp
        submodules/solvers/src/solver/Hybrid.hpp
        submodules/solvers/src/solver/WavePropagation.hpp
        )

include_directories(src/tools src/writer submodules/solvers/src/solver)

add_executable(SWE1D ${SOURCE_FILES})