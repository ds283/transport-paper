//
// Created by David Seery on 25/06/2016.
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

#include "transport-runtime/transport.h"
#include "new_axion_mpi.h"


void write_tasks(transport::repository<>& repo, transport::new_axion_mpi<>* model);


void write_tasks(transport::repository<>& repo, transport::new_axion_mpi<>* model)
  {
    const double M_Planck = 1.0;
    const double g        = 1e-10;
    const double f        = M_Planck;
    const double Lambda   = std::pow(g, 1.0/4.0) * std::pow(25.0/(2*M_PI), 1.0/2.0) * M_Planck;

    const double phi_init = 23.5 * M_Planck;
    const double chi_init = f/2.0 - 0.001*M_Planck;

    const double N_init = 0.0;
    const double N_pre  = 10.0;
    const double N_max  = 68.0;

    transport::parameters<> params(M_Planck, { g, Lambda, f, M_PI }, model);
    transport::initial_conditions<> ics("axion_shape", params, { phi_init, chi_init }, N_init, N_pre);

    transport::basic_range<> times(N_init, N_max, 100, transport::spacing::linear);


    transport::basic_range<> ks(exp(5.0), exp(5.0), 0, transport::spacing::linear);
    transport::basic_range<> alphas(-0.98, 0.98, 98, transport::spacing::linear);
    transport::basic_range<> betas(0.0, 0.99, 99, transport::spacing::linear);


    struct StoragePolicy
      {
      public:
        transport::storage_outcome operator()(const transport::threepf_kconfig& data) { return transport::storage_outcome::accept; }
      };

    struct TrianglePolicy
      {
      public:
        bool operator()(double alpha, double beta)
          {
            // require only triangle condition, but *not* ordering condition

            // beta should lie between 0 and 1
            if(beta < 0.0) return false;
            if(beta > 1.0) return false;

            // alpha should lie between 1-beta and beta-1
            if(beta - 1.0 - alpha > 1E-6) return false;
            if(alpha - (1.0 - beta) > 1E-6) return false;

            // demand that squeezing should not be too small
            if(std::abs(1.0 - beta) < 1E-6) return false;
            if(std::abs(1.0 + alpha + beta) < 1E-6) return false;
            if(std::abs(1.0 - alpha + beta) < 1E-6) return false;

            return true;
          }
      };

    // construct a threepf task
    transport::threepf_alphabeta_task<> tk3("axion_shape.threepf", ics, times, ks, alphas, betas, false, StoragePolicy(), TrianglePolicy());
    tk3.set_adaptive_ics_efolds(4.0);
    tk3.set_collect_initial_conditions(true);

    // construct a zeta threepf task, paired with the primary integration task
    transport::zeta_threepf_task<> ztk3("axion_shape.threepf-zeta", tk3);
    ztk3.set_paired(true);

    repo.commit(ztk3);
  }


int main(int argc, char* argv[])
  {
    // create task manager instance
    transport::task_manager<> mgr(argc, argv);

    // create model instance
    std::shared_ptr< transport::new_axion_mpi<> > model = mgr.create_model< transport::new_axion_mpi<> >();

    // write tasks to repository
    mgr.add_generator([=](transport::repository<>& repo) -> void { write_tasks(repo, model.get()); });

    // hand off control to task manager
    mgr.process();

    return(EXIT_SUCCESS);
  }