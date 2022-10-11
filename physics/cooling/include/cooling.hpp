/*
 * MIT License
 *
 * Copyright (c) 2022 CSCS, ETH Zurich
 *               2022 University of Basel
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/*! @file
 * @brief Interface to the GRACKLE library for radiative cooling
 *
 * @author Noah Kubli <noah.kubli@uzh.ch>
 */

#pragma once

#include <cmath>
#include <optional>
#include <cstring>
#include "cstone/util/util.hpp"

#define CONFIG_BFLOAT_8

extern "C"
{
#include "grackle.h"
}

namespace cooling
{

//Struct containing all cooling information, the same for all particles
template<typename T>
struct cooling_data
{
    //These values need to be passed to GRACKLE at each call
    struct global_values
    {
        code_units             units;
        chemistry_data         data;
        chemistry_data_storage rates;
    } global_values;

    cooling_data(const std::string &grackle_options_file_path,
                 const double ms_sim = 1e16,
                 const double kp_sim = 46400.,
                 const int comoving_coordinates = 0,
                 const std::optional<double> t_sim = std::nullopt);
    
    cooling_data(const cooling_data&) = delete;

    //Solar mass in g
    constexpr static T      ms_g = 1.989e33;
    //kpc in cm
    constexpr static T      kp_cm = 3.086e21;
    //Gravitational constant in cgs units
    constexpr static T      G_newton = 6.674e-8;
    //code unit mass in solar masses
    T                       ms = 1e16;
    //code unit length in kpc
    T                       kpc = 46400.;
    //Path to Grackle data file
    std::string             grackle_data_file_path;

    //Indices to access fields
    struct gr_data
    {
        enum IndexNames
        {
            HI_fraction = 0,
            HII_fraction  = 1,
            HM_fraction  = 3,
            HeI_fraction  = 4,
            HeII_fraction  = 5,
            HeIII_fraction  = 6,
            H2I_fraction = 7,
            H2II_fraction  = 8,
            DI_fraction  = 9,
            DII_fraction  = 10,
            HDI_fraction  = 11,
            e_fraction  = 12,
            metal_fraction  = 13,
            volumetric_heating_rate = 14,
            specific_heating_rate  = 15,
            RT_heating_rate  = 16,
            RT_HI_ionization_rate  = 17,
            RT_HeI_ionization_rate  = 18,
            RT_HeII_ionization_rate  = 19,
            RT_H2_dissociation_rate  = 20,
            H2_self_shielding_length = 21
        };
    };
};


template<typename T>
cooling_data<T>::cooling_data(const std::string &grackle_options_file_path,
                              const double ms_sim, //Solar mass in simulation units
                              const double kp_sim, //Kpc in simulation units
                              const int comoving_coordinates,
                              const std::optional<double> t_sim) //sec in simulation units, just used for the unit test, otherwise G=1)
    : ms(ms_sim), kpc(kp_sim)
{

    grackle_verbose = 1;

    //Density
    const double density_unit = ms * ms_g / std::pow(kpc * kp_cm, 3);
    //Time
    const double time_unit = t_sim.value_or(std::sqrt(1. / (density_unit * G_newton)));
    //Length
    const double length_unit = kpc * kp_cm;
    //Velocity
    const double velocity_unit = length_unit / time_unit;


    global_values.units.density_units        = density_unit; // m_sun / (pc * pc * pc);
    global_values.units.time_units           = time_unit;    // code_time;
    global_values.units.length_units         = length_unit;  // pc;
    global_values.units.velocity_units       = velocity_unit;
    global_values.units.a_units              = 1.0;
    global_values.units.a_value              = 1.0;
    global_values.units.comoving_coordinates = comoving_coordinates;



    std::cout <<"debug\n";
    std::cout << ms << "\t" << ms_g << "\t" <<  kpc << "\n";
    std::cout <<"code units\n";
    std::cout << global_values.units.density_units <<"\t" <<  global_values.units.time_units
              << "\t" <<
              global_values.units.length_units << "\n";

    global_values.data = _set_default_chemistry_parameters();

    grackle_data_file_path = PROJECT_SOURCE_DIR "/extern/grackle/grackle_repo/input/CloudyData_UVB=HM2012.h5";

    global_values.data.grackle_data_file = &grackle_data_file_path[0];

    auto setGrackleOption = [&](chemistry_data &data, const std::string &key, const std::string &value)
    {
      if (key.find("grackle_data_file_path") != std::string::npos) {
          grackle_data_file_path = std::string(value);
          data.grackle_data_file = &grackle_data_file_path[0];
      }
      if (key.find("with_radiative_cooling") != std::string::npos)
          data.with_radiative_cooling = std::stoi(value);
      if (key.find("primordial_chemistry") != std::string::npos)
          data.primordial_chemistry = std::stoi(value);
      if (key.find("h2_on_dust") != std::string::npos)
          data.h2_on_dust = std::stoi(value);
      if (key.find("metal_cooling") != std::string::npos)
          data.metal_cooling = std::stoi(value);
      if (key.find("cmb_temperature_floor") != std::string::npos)
          data.cmb_temperature_floor = std::stoi(value);
      if (key.find("UVbackground") != std::string::npos)
          data.UVbackground = std::stoi(value);
      if (key.find("UVbackground_redshift_on") != std::string::npos)
          data.UVbackground_redshift_on = std::stoi(value);
      if (key.find("UVbackground_redshift_fullon") != std::string::npos)
          data.UVbackground_redshift_fullon = std::stoi(value);
      if (key.find("UVbackground_redshift_drop") != std::string::npos)
          data.UVbackground_redshift_drop = std::stoi(value);
      if (key.find("UVbackground_redshift_off") != std::string::npos)
          data.UVbackground_redshift_off = std::stoi(value);
      if (key.find("Gamma") != std::string::npos)
          data.Gamma = std::stoi(value);
      if (key.find("three_body_rate") != std::string::npos)
          data.three_body_rate = std::stoi(value);
      if (key.find("cie_cooling") != std::string::npos)
          data.cie_cooling = std::stoi(value);
      if (key.find("h2_optical_depth_approximation") != std::string::npos)
          data.h2_optical_depth_approximation = std::stoi(value);
      if (key.find("photoelectric_heating_rate") != std::string::npos)
          data.photoelectric_heating_rate = std::stoi(value);
      if (key.find("Compton_xray_heating") != std::string::npos)
          data.Compton_xray_heating = std::stoi(value);
      if (key.find("LWbackground_intensity") != std::string::npos)
          data.LWbackground_intensity = std::stoi(value);
      if (key.find("LWbackground_sawtooth_suppression") != std::string::npos)
          data.LWbackground_sawtooth_suppression = std::stoi(value);
      if (key.find("use_volumetric_heating_rate") != std::string::npos)
          data.use_volumetric_heating_rate = std::stoi(value);
      if (key.find("use_specific_heating_rate") != std::string::npos)
          data.use_specific_heating_rate = std::stoi(value);
      if (key.find("use_radiative_transfer") != std::string::npos)
          data.use_radiative_transfer = std::stoi(value);
      if (key.find("radiative_transfer_coupled_rate_solver") != std::string::npos)
          data.radiative_transfer_coupled_rate_solver = std::stoi(value);
      if (key.find("radiative_transfer_intermediate_step") != std::string::npos)
          data.radiative_transfer_intermediate_step = std::stoi(value);
      if (key.find("radiative_transfer_hydrogen_only") != std::string::npos)
          data.radiative_transfer_hydrogen_only = std::stoi(value);
      if (key.find("H2_self_shielding") != std::string::npos)
          data.H2_self_shielding = std::stoi(value);
      if (key.find("dust_chemistry") != std::string::npos)
          data.dust_chemistry = std::stoi(value);
    };



    FILE* file = fopen(grackle_options_file_path.c_str(), "r");
    char key[32], value[64];

    while (fscanf(file, "%31s = %63s", key, value) == 2) {
        setGrackleOption(global_values.data, std::string(key), std::string(value));
    }
    std::cout << global_values.data.grackle_data_file << std::endl;
    fclose(file);

    global_values.data.use_grackle = 1;

    if (0 == _initialize_chemistry_data(&global_values.data,
                                    &global_values.rates,
                                    &global_values.units)) {
        throw std::runtime_error("Grackle: Error in _initialize_chemistry_data");
    };
}



// Initialize Grackle chemistry arrays with default data
template<typename Dataset>
void initGrackleData(Dataset& d, typename cooling_data<typename Dataset::RealType>::global_values& gv)
{

    using T = typename Dataset::RealType;
    // This is done so in the sample implementation from GRACKLE – don't know if really needed
    constexpr T tiny_number = 1.e-20;

    constexpr T metal_fraction     = 0.0;
    const T     non_metal_fraction = 1. - metal_fraction;

    util::array<T, 21> gr_test;
    std::fill(gr_test.begin(), gr_test.end(), 0.);
    using gr_data                    = typename cooling_data<T>::gr_data;
    gr_test[gr_data::HI_fraction]    = non_metal_fraction * gv.data.HydrogenFractionByMass;
    gr_test[gr_data::HeI_fraction]   = non_metal_fraction * (1. - gv.data.HydrogenFractionByMass);
    gr_test[gr_data::DI_fraction]    = 2.0 * 3.4e-5;
    gr_test[gr_data::metal_fraction] = metal_fraction;
    gr_test[gr_data::HII_fraction]   = tiny_number;
    gr_test[gr_data::HeII_fraction]  = tiny_number;
    gr_test[gr_data::HeIII_fraction] = tiny_number;
    gr_test[gr_data::e_fraction]     = tiny_number;
    gr_test[gr_data::HM_fraction]    = tiny_number;
    gr_test[gr_data::H2I_fraction]   = tiny_number;
    gr_test[gr_data::H2II_fraction]  = tiny_number;
    gr_test[gr_data::DII_fraction]   = tiny_number;
    gr_test[gr_data::HDI_fraction]   = tiny_number;

    // size_t first = domain.startIndex();
    // size_t last  = domain.endIndex();
    for (size_t i = 0; i < d.grackleData.size(); i++)
    {
        d.grackleData[i] = gr_test;
    }
}

//
template<typename T>
void cool_particle(typename cooling_data<T>::global_values& gv, const T& dt, T& rho, T& u, T& HI_fraction,
                   T& HII_fraction, T& HM_fraction, T& HeI_fraction, T& HeII_fraction, T& HeIII_fraction,
                   T& H2I_fraction, T& H2II_fraction, T& DI_fraction, T& DII_fraction, T& HDI_fraction, T& e_fraction,
                   T& metal_fraction, T& volumetric_heating_rate, T& specific_heating_rate, T& RT_heating_rate,
                   T& RT_HI_ionization_rate, T& RT_HeI_ionization_rate, T& RT_HeII_ionization_rate,
                   T& RT_H2_dissociation_rate, T& H2_self_shielding_length)
{
    grackle_field_data grackle_fields;
    grackle_fields.grid_rank      = 3;
    int zero[]                    = {0, 0, 0};
    int one[]                     = {1, 1, 1};
    grackle_fields.grid_dimension = one;
    grackle_fields.grid_start     = zero;
    grackle_fields.grid_end       = zero;
    grackle_fields.grid_dx        = 0.0;

    gr_float gr_rho                      = (gr_float)rho;
    grackle_fields.density               = &gr_rho;
    gr_float gr_u                        = (gr_float)u;
    grackle_fields.internal_energy       = &gr_u;
    gr_float x_velocity                  = 0.;
    grackle_fields.x_velocity            = &x_velocity;
    gr_float y_velocity                  = 0.;
    grackle_fields.y_velocity            = &y_velocity;
    gr_float z_velocity                  = 0.;
    grackle_fields.z_velocity            = &z_velocity;
    gr_float HI_density                  = (gr_float)HI_fraction * (gr_float)rho;
    gr_float HII_density                 = (gr_float)HII_fraction * (gr_float)rho;
    gr_float HM_density                  = (gr_float)HM_fraction * (gr_float)rho;
    gr_float HeI_density                 = (gr_float)HeI_fraction * (gr_float)rho;
    gr_float HeII_density                = (gr_float)HeII_fraction * (gr_float)rho;
    gr_float HeIII_density               = (gr_float)HeIII_fraction * (gr_float)rho;
    gr_float H2I_density                 = (gr_float)H2I_fraction * (gr_float)rho;
    gr_float H2II_density                = (gr_float)H2II_fraction * (gr_float)rho;
    gr_float DI_density                  = (gr_float)DI_fraction * (gr_float)rho;
    gr_float DII_density                 = (gr_float)DII_fraction * (gr_float)rho;
    gr_float HDI_density                 = (gr_float)HDI_fraction * (gr_float)rho;
    gr_float e_density                   = (gr_float)e_fraction * (gr_float)rho;
    gr_float metal_density               = (gr_float)metal_fraction * (gr_float)rho;
    gr_float volumetric_heating_rate_gr  = (gr_float)volumetric_heating_rate;
    gr_float specific_heating_rate_gr    = (gr_float)specific_heating_rate;
    gr_float RT_heating_rate_gr          = (gr_float)RT_heating_rate;
    gr_float RT_HI_ionization_rate_gr    = (gr_float)RT_HI_ionization_rate;
    gr_float RT_HeI_ionization_rate_gr   = (gr_float)RT_HeI_ionization_rate;
    gr_float RT_HeII_ionization_rate_gr  = (gr_float)RT_HeII_ionization_rate;
    gr_float RT_H2_dissociation_rate_gr  = (gr_float)RT_H2_dissociation_rate;
    gr_float H2_self_shielding_length_gr = (gr_float)H2_self_shielding_length;

    grackle_fields.HI_density    = &HI_density;
    grackle_fields.HII_density   = &HII_density;
    grackle_fields.HeI_density   = &HeI_density;
    grackle_fields.HeII_density  = &HeII_density;
    grackle_fields.HeIII_density = &HeIII_density;
    grackle_fields.e_density     = &e_density;
    grackle_fields.HM_density    = &HM_density;
    grackle_fields.H2I_density   = &H2I_density;
    grackle_fields.H2II_density  = &H2II_density;
    grackle_fields.DI_density    = &DI_density;
    grackle_fields.DII_density   = &DII_density;
    grackle_fields.HDI_density   = &HDI_density;
    grackle_fields.metal_density = &metal_density;

    grackle_fields.volumetric_heating_rate  = &volumetric_heating_rate_gr;
    grackle_fields.specific_heating_rate    = &specific_heating_rate_gr;
    grackle_fields.RT_heating_rate          = &RT_heating_rate_gr;
    grackle_fields.RT_HI_ionization_rate    = &RT_HI_ionization_rate_gr;
    grackle_fields.RT_HeI_ionization_rate   = &RT_HeI_ionization_rate_gr;
    grackle_fields.RT_HeII_ionization_rate  = &RT_HeII_ionization_rate_gr;
    grackle_fields.RT_H2_dissociation_rate  = &RT_H2_dissociation_rate_gr;
    grackle_fields.H2_self_shielding_length = &H2_self_shielding_length_gr;

    //Grackle uses 0 as a return code to indicate failure
    if (0 == local_solve_chemistry(&gv.data, &gv.rates, &gv.units, &grackle_fields, dt)) {
        throw std::runtime_error("Grackle: Error in local_solve_chemistry");
    }
    rho                      = gr_rho;
    u                        = gr_u;
    HI_fraction              = HI_density / gr_rho;
    HII_fraction             = HII_density / gr_rho;
    HM_fraction              = HM_density / gr_rho;
    HeI_fraction             = HeI_density / gr_rho;
    HeII_fraction            = HeII_density / gr_rho;
    HeIII_fraction           = HeIII_density / gr_rho;
    H2I_fraction             = H2I_density / gr_rho;
    H2II_fraction            = H2II_density / gr_rho;
    DI_fraction              = DI_density / gr_rho;
    DII_fraction             = DII_density / gr_rho;
    HDI_fraction             = HDI_density / gr_rho;
    e_fraction               = e_density / gr_rho;
    metal_fraction           = metal_density / gr_rho;
    volumetric_heating_rate  = volumetric_heating_rate_gr;
    specific_heating_rate    = specific_heating_rate_gr;
    RT_heating_rate          = RT_heating_rate_gr;
    RT_HI_ionization_rate    = RT_HI_ionization_rate_gr;
    RT_HeI_ionization_rate   = RT_HeI_ionization_rate_gr;
    RT_HeII_ionization_rate  = RT_HeII_ionization_rate_gr;
    RT_H2_dissociation_rate  = RT_H2_dissociation_rate_gr;
    H2_self_shielding_length = H2_self_shielding_length_gr;
}

// Function not used now but may be needed for initializing the internal energy
template<typename T>
T energy_to_temperature(typename cooling_data<T>::global_values& gv, const T& dt, T& rho, T& u, T& HI_fraction,
                        T& HII_fraction, T& HM_fraction, T& HeI_fraction, T& HeII_fraction, T& HeIII_fraction,
                        T& H2I_fraction, T& H2II_fraction, T& DI_fraction, T& DII_fraction, T& HDI_fraction,
                        T& e_fraction, T& metal_fraction, T& volumetric_heating_rate, T& specific_heating_rate,
                        T& RT_heating_rate, T& RT_HI_ionization_rate, T& RT_HeI_ionization_rate,
                        T& RT_HeII_ionization_rate, T& RT_H2_dissociation_rate, T& H2_self_shielding_length)
{
    grackle_field_data grackle_fields;
    grackle_fields.grid_rank      = 3;
    int zero[]                    = {0, 0, 0};
    int one[]                     = {1, 1, 1};
    grackle_fields.grid_dimension = one;
    grackle_fields.grid_start     = zero;
    grackle_fields.grid_end       = zero;
    grackle_fields.grid_dx        = 0.0;

    gr_float gr_rho                      = (gr_float)rho;
    grackle_fields.density               = &gr_rho;
    gr_float gr_u                        = (gr_float)u;
    grackle_fields.internal_energy       = &gr_u;
    gr_float x_velocity                  = 0.;
    grackle_fields.x_velocity            = &x_velocity;
    gr_float y_velocity                  = 0.;
    grackle_fields.y_velocity            = &y_velocity;
    gr_float z_velocity                  = 0.;
    grackle_fields.z_velocity            = &z_velocity;
    gr_float HI_density                  = (gr_float)HI_fraction * (gr_float)rho;
    gr_float HII_density                 = (gr_float)HII_fraction * (gr_float)rho;
    gr_float HM_density                  = (gr_float)HM_fraction * (gr_float)rho;
    gr_float HeI_density                 = (gr_float)HeI_fraction * (gr_float)rho;
    gr_float HeII_density                = (gr_float)HeII_fraction * (gr_float)rho;
    gr_float HeIII_density               = (gr_float)HeIII_fraction * (gr_float)rho;
    gr_float H2I_density                 = (gr_float)H2I_fraction * (gr_float)rho;
    gr_float H2II_density                = (gr_float)H2II_fraction * (gr_float)rho;
    gr_float DI_density                  = (gr_float)DI_fraction * (gr_float)rho;
    gr_float DII_density                 = (gr_float)DII_fraction * (gr_float)rho;
    gr_float HDI_density                 = (gr_float)HDI_fraction * (gr_float)rho;
    gr_float e_density                   = (gr_float)e_fraction * (gr_float)rho;
    gr_float metal_density               = (gr_float)metal_fraction * (gr_float)rho;
    gr_float volumetric_heating_rate_gr  = (gr_float)volumetric_heating_rate;
    gr_float specific_heating_rate_gr    = (gr_float)specific_heating_rate;
    gr_float RT_heating_rate_gr          = (gr_float)RT_heating_rate;
    gr_float RT_HI_ionization_rate_gr    = (gr_float)RT_HI_ionization_rate;
    gr_float RT_HeI_ionization_rate_gr   = (gr_float)RT_HeI_ionization_rate;
    gr_float RT_HeII_ionization_rate_gr  = (gr_float)RT_HeII_ionization_rate;
    gr_float RT_H2_dissociation_rate_gr  = (gr_float)RT_H2_dissociation_rate;
    gr_float H2_self_shielding_length_gr = (gr_float)H2_self_shielding_length;

    grackle_fields.HI_density    = &HI_density;
    grackle_fields.HII_density   = &HII_density;
    grackle_fields.HeI_density   = &HeI_density;
    grackle_fields.HeII_density  = &HeII_density;
    grackle_fields.HeIII_density = &HeIII_density;
    grackle_fields.e_density     = &e_density;
    grackle_fields.HM_density    = &HM_density;
    grackle_fields.H2I_density   = &H2I_density;
    grackle_fields.H2II_density  = &H2II_density;
    grackle_fields.DI_density    = &DI_density;
    grackle_fields.DII_density   = &DII_density;
    grackle_fields.HDI_density   = &HDI_density;
    grackle_fields.metal_density = &metal_density;

    grackle_fields.volumetric_heating_rate  = &volumetric_heating_rate_gr;
    grackle_fields.specific_heating_rate    = &specific_heating_rate_gr;
    grackle_fields.RT_heating_rate          = &RT_heating_rate_gr;
    grackle_fields.RT_HI_ionization_rate    = &RT_HI_ionization_rate_gr;
    grackle_fields.RT_HeI_ionization_rate   = &RT_HeI_ionization_rate_gr;
    grackle_fields.RT_HeII_ionization_rate  = &RT_HeII_ionization_rate_gr;
    grackle_fields.RT_H2_dissociation_rate  = &RT_H2_dissociation_rate_gr;
    grackle_fields.H2_self_shielding_length = &H2_self_shielding_length_gr;

    gr_float temp;

    if (0 == local_calculate_temperature(&gv.data, &gv.rates, &gv.units, &grackle_fields, &temp) == 0)
    {
        throw std::runtime_error("Grackle: Error in local_calculate_temperature");
    }
    return temp;
}

}; //namespace cooling










