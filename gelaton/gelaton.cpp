//
// Created by David Seery on 27/06/2016.
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
#include "gelaton_mpi.h"


void write_tasks(transport::repository<>& repo, transport::gelaton_mpi<>* model);


void write_tasks(transport::repository<>& repo, transport::gelaton_mpi<>* model)
  {
    const double M_P       = 1.0;
    const double M_chi     = std::sqrt(10.0) * M_P;
    const double epsilon_s = 10.0;

    const double x_init    = -2.0 * M_P;
    const double y_init    = 1E-4 * M_P;

    const double N_init    = 0.0;
    const double N_pre     = 8.0;
    const double N_max     = 29.1;

    transport::parameters<> params(M_P, { M_chi, epsilon_s }, model);
    transport::initial_conditions<> ics("gelaton", params, { x_init, y_init, 0.0, 0.0 }, N_init, N_pre);

    transport::basic_range<> times(N_init, N_max, 100, transport::spacing::linear);


    transport::basic_range<> ks(exp(0.0), exp(20.5), 500, transport::spacing::log_bottom);
    transport::basic_range<> alphas(0.0, 0.0, 0, transport::spacing::linear);
    transport::basic_range<> betas(1.0/3.0, 1.0/3.0, 0, transport::spacing::linear);

    // construct a threepf task
    transport::threepf_alphabeta_task<> tk3("gelaton.threepf", ics, times, ks, alphas, betas);
    tk3.set_collect_initial_conditions(true).set_adaptive_ics_efolds(6.0);

    transport::zeta_threepf_task<> ztk3("gelaton.threepf-zeta", tk3);

    repo.commit(ztk3);
  }


int main(int argc, char* argv[])
  {
    // create task manager instance
    transport::task_manager<> mgr(argc, argv);

    // create model instance
    std::shared_ptr< transport::gelaton_mpi<> > model = mgr.create_model< transport::gelaton_mpi<> >();

    // write tasks to repository
    mgr.add_generator([=](transport::repository<>& repo) -> void { write_tasks(repo, model.get()); });

    // hand off control to task manager
    mgr.process();

    return(EXIT_SUCCESS);
  }
