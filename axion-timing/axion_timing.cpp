///
// Created by David Seery on 08/08/2017.
// --@@
// Copyright (c) 2016 University of Sussex. All rights reserved.
//
// This file is part of the CppTransport platform.
//
// CppTransport is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// CppTransport is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with CppTransport.  If not, see <http://www.gnu.org/licenses/>.
//
// @license: GPL-2
// @contributor: David Seery <D.Seery@sussex.ac.uk>
// --@@
//

#include "axion-canonical_mpi.h"
#include "axion-nontrivial_mpi.h"

using DataType = double;
using StateType = std::vector<DataType>;


template <typename RepositoryType, typename ICType>
transport::zeta_threepf_task<DataType>
make_task(RepositoryType& repo, ICType& ics, const std::string suffix, double N_init, double N_end, double efolds_ics)
  {
    transport::basic_range<double> ts(N_init, N_end, 300, transport::spacing::linear);
    
    const double kt_lo = std::exp(3.0);
    const double kt_hi = std::exp(8.0);
    
    transport::basic_range<double> kts(kt_lo, kt_hi, 1, transport::spacing::linear);
    transport::basic_range<double> alphas(0.0, 0.0, 0, transport::spacing::linear);
    transport::basic_range<double> betas(1.0/3.0, 0.99, 1, transport::spacing::linear);
    
    std::ostringstream task_name;
    task_name << "axion-timing.threepf-" << static_cast<unsigned int>(efolds_ics*100) << "." << suffix;
    
    std::ostringstream zeta_task_name;
    zeta_task_name << "axion-timing.zeta-threepf-" << static_cast<unsigned int>(efolds_ics*100) << "." << suffix;
    
    std::ostringstream task_desc;
    task_desc << "Study timing as a function of massless e-folds for different triangle shapes - " << efolds_ics << " e-folds";
    
    std::ostringstream zeta_task_desc;
    zeta_task_desc << "Zeta 3pf task associated with " << task_name.str() << " - " << efolds_ics << " e-folds";
    
    transport::threepf_alphabeta_task<> tk3(task_name.str(), ics, ts, kts, alphas, betas);
    tk3.set_adaptive_ics_efolds(efolds_ics).set_description(task_desc.str());
    
    transport::zeta_threepf_task<> ztk3(zeta_task_name.str(), tk3);
    ztk3.set_description(zeta_task_desc.str());
    
    return ztk3;
  };


template <typename TaskType, typename SQLQuery>
vis_toolkit::cost_wavenumber<DataType>
make_cost(TaskType& tk, SQLQuery& query, const std::string label)
  {
    vis_toolkit::cost_wavenumber<DataType> cost{
      *dynamic_cast< transport::threepf_task<DataType>* >(tk.get_parent_task()),
      query
    };
    cost.set_label_text(label, label);
    
    return cost;
  }


template <typename RepositoryType, typename Model>
vis_toolkit::wavenumber_series_table<DataType>
make_cost_table(RepositoryType& repo, Model m, std::string suffix, std::string table_name, std::string table_file)
  {
    const double M_Planck = 1.0;
    const double g        = 1e-10;
    const double f        = M_Planck;
    const double Lambda   = std::pow(g, 1.0/4.0) * std::pow(25.0/(2*M_PI), 1.0/2.0) * M_Planck;
    
    const double phi_init = 23.5 * M_Planck;
    const double chi_init = f/2.0 - 0.001*M_Planck;
    
    const double N_init = 0.0;
    const double N_pre  = 17.0;
    const double N_max  = 68.0;
    
    transport::parameters<DataType> params(M_Planck, { g, Lambda, f, M_PI }, m);
    
    std::ostringstream ics_name;
    ics_name << "axion-timing-" << suffix;
    
    transport::initial_conditions<DataType> ics(ics_name.str(), params, { phi_init, chi_init }, N_init, N_pre);
    
    
    auto ztk300 = make_task(repo, ics, suffix, N_init, N_max, 3.00);
    auto ztk325 = make_task(repo, ics, suffix, N_init, N_max, 3.25);
    auto ztk350 = make_task(repo, ics, suffix, N_init, N_max, 3.50);
    auto ztk375 = make_task(repo, ics, suffix, N_init, N_max, 3.75);
    
    auto ztk400 = make_task(repo, ics, suffix, N_init, N_max, 4.00);
    auto ztk425 = make_task(repo, ics, suffix, N_init, N_max, 4.25);
    auto ztk450 = make_task(repo, ics, suffix, N_init, N_max, 4.50);
    auto ztk475 = make_task(repo, ics, suffix, N_init, N_max, 4.75);
    
    auto ztk500 = make_task(repo, ics, suffix, N_init, N_max, 5.00);
    auto ztk525 = make_task(repo, ics, suffix, N_init, N_max, 5.25);
    auto ztk550 = make_task(repo, ics, suffix, N_init, N_max, 5.50);
    auto ztk575 = make_task(repo, ics, suffix, N_init, N_max, 5.75);
    
    auto ztk600 = make_task(repo, ics, suffix, N_init, N_max, 6.00);
    auto ztk625 = make_task(repo, ics, suffix, N_init, N_max, 6.25);
    auto ztk650 = make_task(repo, ics, suffix, N_init, N_max, 6.50);
    auto ztk675 = make_task(repo, ics, suffix, N_init, N_max, 6.75);
    
    auto ztk700 = make_task(repo, ics, suffix, N_init, N_max, 7.00);
    auto ztk725 = make_task(repo, ics, suffix, N_init, N_max, 7.25);
    auto ztk750 = make_task(repo, ics, suffix, N_init, N_max, 7.50);
    auto ztk775 = make_task(repo, ics, suffix, N_init, N_max, 7.75);
    
    auto ztk800 = make_task(repo, ics, suffix, N_init, N_max, 8.00);
    
    
    vis_toolkit::SQL_time_query last_time("serial IN (SELECT MAX(serial) FROM time_samples)");
    vis_toolkit::SQL_threepf_query all_threepfs("1=1");
    
    
    auto cost300 = make_cost(ztk300, all_threepfs, "3.00");
    auto cost325 = make_cost(ztk325, all_threepfs, "3.25");
    auto cost350 = make_cost(ztk350, all_threepfs, "3.50");
    auto cost375 = make_cost(ztk375, all_threepfs, "3.75");
    
    auto cost400 = make_cost(ztk400, all_threepfs, "4.00");
    auto cost425 = make_cost(ztk425, all_threepfs, "4.25");
    auto cost450 = make_cost(ztk450, all_threepfs, "4.50");
    auto cost475 = make_cost(ztk475, all_threepfs, "4.75");
    
    auto cost500 = make_cost(ztk500, all_threepfs, "5.00");
    auto cost525 = make_cost(ztk525, all_threepfs, "5.25");
    auto cost550 = make_cost(ztk550, all_threepfs, "5.50");
    auto cost575 = make_cost(ztk575, all_threepfs, "5.75");
    
    auto cost600 = make_cost(ztk600, all_threepfs, "6.00");
    auto cost625 = make_cost(ztk625, all_threepfs, "6.25");
    auto cost650 = make_cost(ztk650, all_threepfs, "6.50");
    auto cost675 = make_cost(ztk675, all_threepfs, "6.75");
    
    auto cost700 = make_cost(ztk700, all_threepfs, "7.00");
    auto cost725 = make_cost(ztk725, all_threepfs, "7.25");
    auto cost750 = make_cost(ztk750, all_threepfs, "7.50");
    auto cost775 = make_cost(ztk775, all_threepfs, "7.75");
    
    auto cost800 = make_cost(ztk800, all_threepfs, "8.00");
    
    
    std::ostringstream tab_name;
    tab_name << table_name << "." << suffix;
    
    vis_toolkit::wavenumber_series_table<DataType> table(tab_name.str(), table_file);
    table +=
      cost300 + cost325 + cost350 + cost375 +
      cost400 + cost425 + cost450 + cost475 +
      cost500 + cost525 + cost550 + cost575 +
      cost600 + cost625 + cost650 + cost675 +
      cost700 + cost725 + cost750 + cost775 +
      cost800;
    
    return table;
  }


template <typename RepositoryType, typename CanonicalModel, typename NontrivialModel>
void write_tasks(RepositoryType& repo, CanonicalModel cmod, NontrivialModel nmod)
  {
    auto canonical_table = make_cost_table(repo, cmod, "can", "canonical-table", "canonical.csv");
    auto nontrivial_table = make_cost_table(repo, nmod, "non", "nontrivial-table", "nontrivial.csv");

    transport::output_task<DataType> out_tk("axion-timing.output");
    out_tk += canonical_table + nontrivial_table;
    
    repo.commit(out_tk);
  };


int main(int argc, char* argv[])
  {
    // create task manager instance
    transport::task_manager<DataType> mgr(argc, argv);
    
    // create model instances
    auto canonical_model = mgr.create_model< transport::axion_canonical_mpi<DataType, StateType> >();
    auto nontrivial_model = mgr.create_model< transport::axion_nontrivial_mpi<DataType, StateType> >();
    
    // write tasks to repository
    mgr.add_generator([=](auto& repo) -> void { write_tasks(repo, canonical_model.get(), nontrivial_model.get()); });
    
    // hand off control
    mgr.process();
    
    return EXIT_SUCCESS;
  }
