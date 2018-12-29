#include <iostream>
#include <fstream>
#include <string>

#include "sphexa.hpp"
#include "Evrard.hpp"

using namespace std;
using namespace sphexa;

int neighbors_sum(const std::vector<std::vector<int>> &neighbors)
{
    int sum = 0;
    #pragma omp parallel for reduction (+:sum)
    for(unsigned int i=0; i<neighbors.size(); i++)
        sum += neighbors[i].size();

    #ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    #endif

    return sum;
}

int main()
{
    typedef double Real;
    typedef Octree<Real> Tree;
    typedef Evrard<Real> Dataset;

    #ifdef USE_MPI
        MPI_Init(NULL, NULL);
        Dataset d(1e6, "bigfiles/Evrard3D.bin", MPI_COMM_WORLD);
    #else
        Dataset d(1e6, "bigfiles/Evrard3D.bin");
    #endif

    Domain<Real, Tree> domain(d.ngmin, d.ng0, d.ngmax);
    Density<Real> density(d.K);
    EquationOfState<Real> equationOfState;
    MomentumEnergy<Real> momentumEnergy(d.K);
    Timestep<Real> timestep(d.K, d.maxDtIncrease);
    UpdateQuantities<Real> updateQuantities(d.stabilizationTimesteps);
    EnergyConservation<Real> energyConservation;

    for(int iteration = 0; iteration < 200; iteration++)
    {
        timer::TimePoint start = timer::Clock::now();

        if(d.rank == 0) cout << "Iteration: " << iteration << endl;

        REPORT_TIME(domain.buildTree(d.x, d.y, d.z, d.h, d.bbox), "BuildTree");
        // REPORT_TIME(domain.reorder(d.data), "ReorderParticles");
        REPORT_TIME(domain.findNeighbors(d.x, d.y, d.z, d.h, d.neighbors), "FindNeighbors");
        // REPORT_TIME(density.compute(d.neighbors, d.x, d.y, d.z, d.h, d.m, d.ro), "Density");
        // REPORT_TIME(equationOfState.compute(d.ro, d.mui, d.temp, d.u, d.p, d.c, d.cv), "EquationOfState");
        // REPORT_TIME(momentumEnergy.compute(d.neighbors, d.x, d.y, d.z, d.h, d.vx, d.vy, d.vz, d.ro, d.p, d.c, d.m, d.grad_P_x, d.grad_P_y, d.grad_P_z, d.du), "MomentumEnergy");
        // REPORT_TIME(timestep.compute(d.h, d.c, d.dt_m1, d.dt), "Timestep");
        // REPORT_TIME(updateQuantities.compute(iteration, d.grad_P_x, d.grad_P_y, d.grad_P_z, d.dt, d.du, d.bbox, d.x, d.y, d.z, d.vx, d.vy, d.vz, d.x_m1, d.y_m1, d.z_m1, d.u, d.du_m1, d.dt_m1), "UpdateQuantities");
        // REPORT_TIME(energyConservation.compute(d.u, d.vx, d.vy, d.vz, d.m, d.etot, d.ecin, d.eint), "EnergyConservation");

        int totalNeighbors = neighbors_sum(d.neighbors);

        if(d.rank == 0)
        {
            cout << "### Check ### Computational domain: ";
            cout << d.bbox.xmin << " " << d.bbox.xmax << " ";
            cout << d.bbox.ymin << " " << d.bbox.ymax << " ";
            cout << d.bbox.zmin << " " << d.bbox.zmax << endl;
            cout << "### Check ### Avg. number of neighbours: " << totalNeighbors/d.n << "/" << d.ng0 << endl;
            cout << "### Check ### New Time-step: " << d.dt[0] << endl;
            cout << "### Check ### Total energy: " << d.etot << ", (internal: " << d.eint << ", cinetic: " << d.ecin << ")" << endl;

            if(iteration % 10 == 0)
            {
                std::ofstream outputFile("output" + to_string(iteration) + ".txt");
                REPORT_TIME(d.writeFile(outputFile), "writeFile");
                outputFile.close();
            }

            timer::TimePoint stop = timer::Clock::now();
            cout << "=== Total time for iteration " << timer::duration(start, stop) << "s" << endl << endl;
        }
    }

    #ifdef USE_MPI
        MPI_Finalize();
    #endif

    return 0;
}

