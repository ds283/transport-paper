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
    const double eta_chi   = 10.0;
    const double g_chi     = 3.2E5;
    const double epsilon_s = 10.0;

    const double x_init    = -2.0 * M_P;
    const double y_init    = 1E-4 * M_P;

    const double N_init    = 0.0;
    const double N_pre     = 4.0;
    const double N_max     = 29.1;
    
    transport::parameters<> params(M_P, { eta_chi, g_chi, epsilon_s }, model);
    transport::initial_conditions<> ics("gelaton", params, { x_init, y_init, 0.0, 0.0 }, N_init, N_pre);

    transport::basic_range<> times(N_init, N_max, 500, transport::spacing::linear);

    // set up a 2pf job consisting of a single k-mode; this lets us produce the background
    // evolution over the whole range of e-folds
    transport::basic_range<> k(1.0, 1.0, 0, transport::spacing::linear);

    // construct a two task
    transport::twopf_task<> tk2("gelaton.background-history", ics, times, k);
    tk2.set_adaptive_ics(false);

    repo.commit(tk2);
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
