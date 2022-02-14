#pragma once

#include "ifile_writer.hpp"
#include "file_utils.hpp"

namespace sphexa
{
template<typename Dataset>
struct NohDataFileWriter : IFileWriter<Dataset>
{
    void dumpParticleDataToBinFile(const Dataset& d, const std::string& path) const override
    {
        try
        {
            printf("Dumping particles data to file at path: %s\n", path.c_str());
            fileutils::writeParticleDataToBinFile(path,
                                                  d.x,
                                                  d.y,
                                                  d.z,
                                                  d.vx,
                                                  d.vy,
                                                  d.vz,
                                                  d.ro,
                                                  d.u,
                                                  d.p,
                                                  d.h,
                                                  d.m,
                                                  d.temp,
                                                  d.mue,
                                                  d.mui,
                                                  d.du,
                                                  d.du_m1,
                                                  d.dt,
                                                  d.dt_m1,
                                                  d.x_m1,
                                                  d.y_m1,
                                                  d.z_m1,
                                                  d.c,
                                                  d.grad_P_x,
                                                  d.grad_P_y,
                                                  d.grad_P_z);
        }
        catch (FileNotOpenedException& ex)
        {
            fprintf(stderr, "ERROR: %s. Terminating\n", ex.what());
            exit(EXIT_FAILURE);
        }
    }

    void dumpParticleDataToAsciiFile(const Dataset& d, int firstIndex, int lastIndex,
                                     const std::string& path) const override
    {
        try
        {
            const char separator = ' ';

            printf("Dumping particles data to ASCII file at path: %s\n", path.c_str());
            fileutils::writeParticleDataToAsciiFile(firstIndex,
                                                    lastIndex,
                                                    path,
                                                    separator,
                                                    d.x,
                                                    d.y,
                                                    d.z,
                                                    d.vx,
                                                    d.vy,
                                                    d.vz,
                                                    d.ro,
                                                    d.u,
                                                    d.p,
                                                    d.h,
                                                    d.m,
                                                    d.temp,
                                                    d.mue,
                                                    d.mui,
                                                    d.du,
                                                    d.du_m1,
                                                    d.dt,
                                                    d.dt_m1,
                                                    d.x_m1,
                                                    d.y_m1,
                                                    d.z_m1,
                                                    d.c,
                                                    d.grad_P_x,
                                                    d.grad_P_y,
                                                    d.grad_P_z);
        }
        catch (FileNotOpenedException& ex)
        {
            fprintf(stderr, "ERROR: %s. Terminating\n", ex.what());
            exit(EXIT_FAILURE);
        }
    }

    void dumpParticleDataToAsciiFile(const Dataset& d, const std::vector<int>& clist,
                                     const std::string& path) const override
    {
        try
        {
            const char separator = ' ';

            printf("Dumping particles data to ASCII file at path: %s\n", path.c_str());
            fileutils::writeParticleDataToAsciiFile(clist,
                                                    path,
                                                    separator,
                                                    d.x,
                                                    d.y,
                                                    d.z,
                                                    d.vx,
                                                    d.vy,
                                                    d.vz,
                                                    d.ro,
                                                    d.u,
                                                    d.p,
                                                    d.h,
                                                    d.m,
                                                    d.temp,
                                                    d.mue,
                                                    d.mui,
                                                    d.du,
                                                    d.du_m1,
                                                    d.dt,
                                                    d.dt_m1,
                                                    d.x_m1,
                                                    d.y_m1,
                                                    d.z_m1,
                                                    d.c,
                                                    d.grad_P_x,
                                                    d.grad_P_y,
                                                    d.grad_P_z);
        }
        catch (FileNotOpenedException& ex)
        {
            fprintf(stderr, "ERROR: %s. Terminating\n", ex.what());
            exit(EXIT_FAILURE);
        }
    }

    void dumpCheckpointDataToBinFile(const Dataset& d, const std::string& path) const override
    {
        try
        {
            printf("Saving checkpoint at path: %s\n", path.c_str());
            fileutils::writeParticleCheckpointDataToBinFile(d,
                                                            path,
                                                            d.x,
                                                            d.y,
                                                            d.z,
                                                            d.vx,
                                                            d.vy,
                                                            d.vz,
                                                            d.ro,
                                                            d.u,
                                                            d.p,
                                                            d.h,
                                                            d.m,
                                                            d.temp,
                                                            d.mue,
                                                            d.mui,
                                                            d.du,
                                                            d.du_m1,
                                                            d.dt,
                                                            d.dt_m1,
                                                            d.x_m1,
                                                            d.y_m1,
                                                            d.z_m1,
                                                            d.c,
                                                            d.grad_P_x,
                                                            d.grad_P_y,
                                                            d.grad_P_z);
        }
        catch (FileNotOpenedException& ex)
        {
            fprintf(stderr, "ERROR: %s. Terminating\n", ex.what());
            exit(EXIT_FAILURE);
        }
    }
};

#ifdef USE_MPI

template<typename Dataset>
struct NohDataMPIFileWriter : IFileWriter<Dataset>
{
#ifdef SPH_EXA_HAVE_H5PART
    void dumpParticleDataToH5File(const Dataset& d, int firstIndex, int lastIndex,
                                  const std::string& path) const override
    {
        fileutils::writeParticleDataToBinFileWithH5Part(d, firstIndex, lastIndex, path);
    }

    void dumpParticleDataToH5File(const Dataset& d, const std::vector<int>& clist,
                                  const std::string& path) const override
    {
        try
        {
            fileutils::writeParticleDataToBinFileWithH5Part(d, clist, path, d.x, d.y, d.z, d.h, d.ro);
        }
        catch (MPIFileNotOpenedException& ex)
        {
            if (d.rank == 0) fprintf(stderr, "ERROR: %s. Terminating\n", ex.what());
            MPI_Abort(d.comm, ex.mpierr);
        }
    }
#endif

    void dumpParticleDataToBinFile(const Dataset& d, const std::string& path) const override
    {
        try
        {
            fileutils::writeParticleDataToBinFileWithMPI(d,
                                                         path,
                                                         d.x,
                                                         d.y,
                                                         d.z,
                                                         d.vx,
                                                         d.vy,
                                                         d.vz,
                                                         d.ro,
                                                         d.u,
                                                         d.p,
                                                         d.h,
                                                         d.m,
                                                         d.temp,
                                                         d.mue,
                                                         d.mui,
                                                         d.du,
                                                         d.du_m1,
                                                         d.dt,
                                                         d.dt_m1,
                                                         d.x_m1,
                                                         d.y_m1,
                                                         d.z_m1,
                                                         d.c,
                                                         d.grad_P_x,
                                                         d.grad_P_y,
                                                         d.grad_P_z);
        }
        catch (MPIFileNotOpenedException& ex)
        {
            if (d.rank == 0) fprintf(stderr, "ERROR: %s. Terminating\n", ex.what());
            MPI_Abort(d.comm, ex.mpierr);
        }
    }

    void dumpParticleDataToAsciiFile(const Dataset& d, int firstIndex, int lastIndex,
                                     const std::string& path) const override
    {
        const char separator = ' ';

        for (int turn = 0; turn < d.nrank; turn++)
        {
            if (turn == d.rank)
            {
                try
                {
                    fileutils::writeParticleDataToAsciiFile(firstIndex,
                                                            lastIndex,
                                                            path,
                                                            d.rank != 0,
                                                            separator,
                                                            d.x,
                                                            d.y,
                                                            d.z,
                                                            d.vx,
                                                            d.vy,
                                                            d.vz,
                                                            d.ro,
                                                            d.u,
                                                            d.p,
                                                            d.h,
                                                            d.m,
                                                            d.temp,
                                                            d.mue,
                                                            d.mui,
                                                            d.du,
                                                            d.du_m1,
                                                            d.dt,
                                                            d.dt_m1,
                                                            d.x_m1,
                                                            d.y_m1,
                                                            d.z_m1,
                                                            d.c,
                                                            d.grad_P_x,
                                                            d.grad_P_y,
                                                            d.grad_P_z);
                }
                catch (MPIFileNotOpenedException& ex)
                {
                    if (d.rank == 0) fprintf(stderr, "ERROR: %s. Terminating\n", ex.what());
                    MPI_Abort(d.comm, ex.mpierr);
                }

                MPI_Barrier(MPI_COMM_WORLD);
            }
            else
            {
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }
    }

    void dumpCheckpointDataToBinFile(const Dataset& d, const std::string& path) const override
    {
        if (d.rank == 0) printf("Dumping checkpoint...\n");

        dumpParticleDataToBinFile(d, path);

        if (d.rank == 0) printf("Dumped checkpoint...\n");
    }
};
#endif

} // namespace sphexa