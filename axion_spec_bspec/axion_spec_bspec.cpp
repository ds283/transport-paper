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
    transport::initial_conditions<> ics("axion_cfs", params, { phi_init, chi_init }, N_init, N_pre);

    transport::basic_range<> times(N_init, N_max, 500, transport::spacing::linear);


    transport::basic_range<> ks(exp(5.0), exp(5.0), 0, transport::spacing::linear);
    transport::basic_range<> alphas(0.0, 0.0, 0, transport::spacing::linear);

    transport::basic_range<> beta_equi(1.0/3.0, 1.0/3.0, 0, transport::spacing::linear);
    transport::basic_range<> beta_sq(0.95, 0.95, 0, transport::spacing::linear);

    // construct a threepf task for the equilateral mode
    transport::threepf_alphabeta_task<> tk3_equi("axion_cfs.equi.threepf", ics, times, ks, alphas, beta_equi);
    tk3_equi.set_adaptive_ics_efolds(4.0);
    tk3_equi.set_collect_initial_conditions(true);

    transport::zeta_threepf_task<> ztk3_equi("axion_cfs.equi.threepf-zeta", tk3_equi);
    ztk3_equi.set_paired(true);

    // construct a threepf task for the squeezed mode
    transport::threepf_alphabeta_task<> tk3_sq("axion_cfs.sq.threepf", ics, times, ks, alphas, beta_sq);
    tk3_sq.set_adaptive_ics_efolds(4.0);
    tk3_sq.set_collect_initial_conditions(true);

    transport::zeta_threepf_task<> ztk3_sq("axion_cfs.sq.threepf-zeta", tk3_sq);
    ztk3_sq.set_paired(true);

    repo.commit(ztk3_equi);
    repo.commit(ztk3_sq);
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