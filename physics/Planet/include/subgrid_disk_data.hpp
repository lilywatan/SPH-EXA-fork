//
// Created by Noah Kubli on 07.03.2024.
//

#pragma once

#include <array>
#include <iostream>

struct DiskData
{
    // set disk parameters 
    double                alpha{0.01}; // disk alpha parameter

    // variable disk parameters
    double                H_r{0.25}; // disk H/r ratio
    double                m{0.};
    double                r{0.}; // outer radius of the disk 
    double                h{0.}; // smoothing length of recently accreted particles
    double                c_boundary{0.}; // sound speed of the disk
    double                r0{1.}; // radius of the boundary layer -> should be 2 or 3 times h 
    double                rho_boundary{0.}; // density at outer boundary of disk
    double                T_boundary{0.}; // temperature at outer boundary of disk
    double                sigma0{0.}; // surface density at outer boundary of disk (sigma naught)
    double                Hr_boundary{0.}; // kinematic viscosity of the disk


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
        optionalIO("disk::H_r", &H_r, 1);
        optionalIO("disk::alpha", &alpha, 1);
        optionalIO("disk::m", &m, 1);
        optionalIO("disk::c_boundary", &c_boundary, 1);
        optionalIO("disk::r0", &r0, 1);
        optionalIO("disk::r", &r, 1);
        optionalIO("disk::h", &h, 1);
        optionalIO("disk::rho_boundary", &rho_boundary, 1);
        optionalIO("disk::T_boundary", &T_boundary, 1);
        optionalIO("disk::sigma0", &sigma0, 1);
        optionalIO("disk::Hr_boundary", &Hr_boundary, 1);
        
    }; 

    // Local to Rank

    double                m_accreted_local{};
    size_t                n_accreted_local{};  
    double                h_accreted_local{};
    double                r_accreted_local{}; 
    double                c_boundary_local{};
    double                rho_boundary_local{};
    double                T_boundary_local{}; 
    double                sigma0_local{}; 
    size_t                n_boundary_local{};
    double                r0_local{};
    double                H2_boundary_local{};
    
};