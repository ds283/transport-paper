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
#include "Gao_mpi.h"


void write_tasks(transport::repository<>& repo, transport::Gao_mpi<>* model);


void write_tasks(transport::repository<>& repo, transport::Gao_mpi<>* model)
  {
    constexpr double M_P        = 1.0;

    constexpr double shift1     = 231 * M_P;
    constexpr double shift2     = 231 * M_P;

    constexpr double Mass       = 1e-4 * M_P;
    constexpr double mphi       = 1e-7 * M_P;
    constexpr double deltaTheta = M_PI / 10.0;
    const     double phi0       = (-100.0 * std::sqrt(6.0)) * M_P + shift1;
    const     double s          = (1000.0 * std::sqrt(3.0)) * M_P;
    constexpr double pi         = M_PI;

    const     double phi_init   = (-2.0-100.0*std::sqrt(6.0)) * M_P + shift2;
    const     double chi_init   = (2.0*std::tan(M_PI/20.0)) * M_P;
    constexpr double dphi_init  = 0.0;
    constexpr double dchi_init  = 0.0;

    constexpr double N_init  = 0.0;
    constexpr double N_pre   = 14.0;
    constexpr double N_max   = 52.0;


    transport::parameters<> params(M_P, { Mass, mphi, deltaTheta, phi0, s, pi }, model);
    transport::initial_conditions<> ics("nonadiabatic", params, { phi_init, chi_init, dphi_init, dchi_init }, N_init, N_pre);

    transport::basic_range<> times(N_init, N_max, 100, transport::spacing::linear);


    transport::basic_range<> k(1.0, 1.0, 0, transport::spacing::linear);

    // construct a twopf task for this single mode
    transport::twopf_task<> tk2("nonadiabatic.history", ics, times, k);
    tk2.set_adaptive_ics_efolds(4.0);

    repo.commit(tk2);
  }


int main(int argc, char* argv[])
  {
    // create task manager instance
    transport::task_manager<> mgr(argc, argv);

    // create model instance
    std::shared_ptr< transport::Gao_mpi<> > model = mgr.create_model< transport::Gao_mpi<> >();

    // write tasks to repository
    mgr.add_generator([=](transport::repository<>& repo) -> void { write_tasks(repo, model.get()); });

    // hand off control to task manager
    mgr.process();

    return(EXIT_SUCCESS);
  }
