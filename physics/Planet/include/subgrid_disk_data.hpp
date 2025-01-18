//
// Created by Noah Kubli on 07.03.2024.
//

#pragma once

#include <array>
#include <iostream>

struct DiskData
{
    double                m{0.};
    double                c{0.}; // sound speed of the disk
    double                H_r{0.5}; // disk H/r ratio
    //double                r_in{1e-5}; // inner radius of the disk
    //double                r_out{1.}; // outer radius of the disk -> should be 2 or 3 times h 
    std::vector<double>   surface_density_profile{};
    std::vector<double>   temperature_profile{};
    size_t                num_radial_bins{100};
    //double r_step         {(r_out - r_in) / num_radial_bins};
    double                rho_boundary{0.}; // density at outer boundary of disk
    double                T_boundary{0.}; // temperature at outer boundary of disk
    double                sigma_boundary{0.}; // surface density at outer boundary of disk


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
        optionalIO("disk::m", &m, 1);
        optionalIO("disk::c", &c 1);
        optionalIO("disk::H", &H, 1);
        optionalIO("disk::num_radial_bins", &num_radial_bins, 1);
        //optionalIO("disk::r_out", &r_out, 1);
        optionalIO("disk::surface_density_profile", surface_density_profile.data(), num_radial_bins);
        optionalIO("disk::temperature_profile", temperature_profile.data(), num_radial_bins); 
        optionalIO("disk::rho_boundary", &rho_boundary, 1);
        optionalIO("disk::T_boundary", &T_boundary, 1);
        optionalIO("disk::sigma_boundary", &sigma_boundary, 1);
        
    }; 

    // Local to Rank

    double                m_accreted_local_subdisk{};
    double                rho_boundary_local{};
    double                temp_boundary_local{}; 
    double                sigma_boundary_local{};
    double                sound_speed_boundary_local{}; 
    size_t                n_boundary_particles{};
};