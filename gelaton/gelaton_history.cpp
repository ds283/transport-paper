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
    const double M_P           = 1.0;
    const double P_zeta        = 1E-9;                                                          // desired amplitude of fluctuations
    const double omega         = M_PI / 30.0;                                                   // try to get round pi radians in ~ 30 e-folds
    
    const double V0            = 0.1 * P_zeta * M_P*M_P*M_P*M_P;                                // adjust uplift to get sufficient inflation
    const double eta_R         = 1.0/std::sqrt(3.0);                                            // adjust radial mass to be of order Hubble
    const double g_R           = M_P*M_P / std::sqrt(V0);                                       // adjust radial cubic coupling to be of order Hubble
    const double lambda_R      = 0.5 * M_P*M_P*M_P / std::pow(V0, 3.0/4.0) / std::sqrt(omega);  // adjust radial quartic coupling to dominate the displacement
    
    const double alpha         = 7.25*omega;                                                    // adjust angular tilt to get desired omega
    
    const double R0            = std::sqrt(V0/3.0) / (M_P * omega * std::sqrt(P_zeta));         // adjust radial minimum to give desired P_zeta normalization
    
    const double x_init        = -R0;
    const double y_init        = (1E-2)*R0;
    
    const double N_init        = 0.0;
    const double N_pre         = 8.0;
    const double N_max         = 28.0;
    
    transport::parameters<> params(M_P, { R0, V0, eta_R, g_R, lambda_R, alpha }, model);
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
