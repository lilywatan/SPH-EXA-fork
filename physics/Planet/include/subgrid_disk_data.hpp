//
// Created by Noah Kubli on 07.03.2024.
//

#pragma once

#include <array>
#include <iostream>

struct DiskData
{
    // std::array<double, 3> position{}; // position of the center of the disk
    // std::array<double, 3> position_m1{};
    double                m{0.};
    double                r_in{1e-5}; // inner radius of the disk
    double                r_out{1.}; // outer radius of the disk -> should be 2 or 3 times h 
    std::vector<double>   surface_density_profile{};
    std::vector<double>   temperature_profile{};
    size_t                num_radial_bins{100};
    double r_step         {(r_out - r_in) / num_radial_bins};
    double                rho_boundary; // density at outer boundary of disk
    double                T_boundary; // temperature at outer boundary of disk
    double                sigma_boundary; // surface density at outer boundary of disk


    // Q: can use this as such even if nothing will be in input file? Or how else to initialize? 
    template<typename Archive>
    void loadOrStoreAttributes(Archive* ar)
    {
        //! @brief load or store an attribute, skips non-existing attributes on load.
        auto optionalIO = [ar](const std::string& attribute, auto* location, size_t attrSize)
        {
            try
            {
                ar->stepAttribute(attribute, location, attrSize);
            }
            catch (std::out_of_range&)
            {
                if (ar->rank() == 0)
                {
                    std::cout << "Attribute " << attribute
                              << " not set in file or initializer, setting to default value " << *location << std::endl;
                }
            }
        };
        //optionalIO("disk::x", &position[0], 1);
        //optionalIO("disk::y", &position[1], 1);
        //optionalIO("disk::z", &position[2], 1);
        //optionalIO("disk::x_m1", &position_m1[0], 1);
        //optionalIO("disk::y_m1", &position_m1[1], 1);
        //optionalIO("disk::z_m1", &position_m1[2], 1);
        optionalIO("disk::num_radial_bins", &num_radial_bins, 1);
        optionalIO("disk::r_out", &r_out, 1);
        optionalIO("disk::surface_density_profile", surface_density_profile.data(), num_radial_bins);
        optionalIO("disk::temperature_profile", temperature_profile.data(), num_radial_bins); 
        
    }; 

    // Local to Rank

    size_t                n_accreted_local_subdisk{};
    size_t                n_removed_local_subdisk{};
    double                m_accreted_local_subdisk{};
    double                m_removed_local_subdisk{};
    std::array<double, 3> p_accreted_local_subdisk{};
    std::array<double, 3> p_removed_local_subdisk{};
    double                boundary_density{};
    double                boundary_temperature{}; 
    double                boundary_sigma{};
    double                boundary_mass{};
    size_t                n_boundary_particles{};
};