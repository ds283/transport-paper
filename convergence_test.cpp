//
// Created by David Seery on 07/01/2014.
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
#include "dq_mpi.h"


void write_tasks(transport::repository<>& repo, transport::dquad_mpi<>* model);

transport::zeta_threepf_task<> make_task(transport::repository<>& repo, transport::initial_conditions<>& ics,
                                         double N_init, double N_end, double efolds_ics);


void write_tasks(transport::repository<>& repo, transport::dquad_mpi<>* model)
  {
    const double Mp = 1.0;
    const double Mphi = 9E-5 * Mp;
    const double Mchi = 1E-5 * Mp;

    transport::parameters<> params(Mp, { Mphi, Mchi }, model);

    const double phi_init = 10.0 * Mp;
    const double chi_init = 12.9 * Mp;

    const double N_init = 0.0;
    const double N_pre = 17.0;
    const double N_end = 60.0;

    transport::initial_conditions<> ics("dquad", params, { phi_init, chi_init }, N_init, N_pre);

    transport::zeta_threepf_task<> ztk300 = make_task(repo, ics, N_init, N_end, 3.00);
    transport::zeta_threepf_task<> ztk325 = make_task(repo, ics, N_init, N_end, 3.25);
    transport::zeta_threepf_task<> ztk350 = make_task(repo, ics, N_init, N_end, 3.50);
    transport::zeta_threepf_task<> ztk375 = make_task(repo, ics, N_init, N_end, 3.75);
    transport::zeta_threepf_task<> ztk400 = make_task(repo, ics, N_init, N_end, 4.00);
    transport::zeta_threepf_task<> ztk425 = make_task(repo, ics, N_init, N_end, 4.25);
    transport::zeta_threepf_task<> ztk450 = make_task(repo, ics, N_init, N_end, 4.50);
    transport::zeta_threepf_task<> ztk475 = make_task(repo, ics, N_init, N_end, 4.75);
    transport::zeta_threepf_task<> ztk500 = make_task(repo, ics, N_init, N_end, 5.00);
    transport::zeta_threepf_task<> ztk525 = make_task(repo, ics, N_init, N_end, 5.25);
    transport::zeta_threepf_task<> ztk550 = make_task(repo, ics, N_init, N_end, 5.50);
    transport::zeta_threepf_task<> ztk575 = make_task(repo, ics, N_init, N_end, 5.75);
    transport::zeta_threepf_task<> ztk600 = make_task(repo, ics, N_init, N_end, 6.00);
    transport::zeta_threepf_task<> ztk625 = make_task(repo, ics, N_init, N_end, 6.25);
    transport::zeta_threepf_task<> ztk650 = make_task(repo, ics, N_init, N_end, 6.50);
    transport::zeta_threepf_task<> ztk675 = make_task(repo, ics, N_init, N_end, 6.75);
    transport::zeta_threepf_task<> ztk700 = make_task(repo, ics, N_init, N_end, 7.00);
    transport::zeta_threepf_task<> ztk725 = make_task(repo, ics, N_init, N_end, 7.25);
    transport::zeta_threepf_task<> ztk750 = make_task(repo, ics, N_init, N_end, 7.50);
    transport::zeta_threepf_task<> ztk775 = make_task(repo, ics, N_init, N_end, 7.75);
    transport::zeta_threepf_task<> ztk800 = make_task(repo, ics, N_init, N_end, 8.00);
    transport::zeta_threepf_task<> ztk825 = make_task(repo, ics, N_init, N_end, 8.25);
    transport::zeta_threepf_task<> ztk850 = make_task(repo, ics, N_init, N_end, 8.50);
    transport::zeta_threepf_task<> ztk875 = make_task(repo, ics, N_init, N_end, 8.75);
    transport::zeta_threepf_task<> ztk900 = make_task(repo, ics, N_init, N_end, 9.00);
    transport::zeta_threepf_task<> ztk925 = make_task(repo, ics, N_init, N_end, 9.25);
    transport::zeta_threepf_task<> ztk950 = make_task(repo, ics, N_init, N_end, 9.50);
    transport::zeta_threepf_task<> ztk975 = make_task(repo, ics, N_init, N_end, 9.75);
    transport::zeta_threepf_task<> ztk1000 = make_task(repo, ics, N_init, N_end, 10.00);

    vis_toolkit::SQL_time_query last_time("tserial IN (SELECT MAX(teserial) FROM time_samples)");
    vis_toolkit::SQL_threepf_query all_threepfs("1=1");

    vis_toolkit::zeta_threepf_time_series<> z300(ztk300, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z325(ztk325, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z350(ztk350, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z375(ztk375, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z400(ztk400, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z425(ztk425, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z450(ztk450, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z475(ztk475, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z500(ztk500, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z525(ztk525, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z550(ztk550, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z575(ztk575, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z600(ztk600, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z625(ztk625, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z650(ztk650, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z675(ztk675, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z700(ztk700, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z725(ztk725, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z750(ztk750, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z775(ztk775, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z800(ztk800, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z825(ztk825, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z850(ztk850, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z875(ztk875, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z900(ztk900, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z925(ztk925, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z950(ztk950, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z975(ztk975, last_time, all_threepfs);
    vis_toolkit::zeta_threepf_time_series<> z1000(ztk1000, last_time, all_threepfs);

    vis_toolkit::time_series_table<> table("table", "table.txt");
    table += z300 + z325 + z350 + z375
      + z400 + z425 + z450 + z475
      + z500 + z525 + z550 + z575
      + z600 + z625 + z650 + z675
      + z700 + z725 + z750 + z775
      + z800 + z825 + z850 + z875
      + z900 + z925 + z950 + z975
      + z1000;

    transport::output_task<> out_tk("convergence.output", table);

    repo.commit(out_tk);
  }


transport::zeta_threepf_task<> make_task(transport::repository<>& repo, transport::initial_conditions<>& ics,
                               double N_init, double N_end, double efolds_ics)
  {
    transport::basic_range<> ts(N_init, N_end, 300, transport::spacing::linear);

    const double kt_lo = std::exp(3.0);
    const double kt_hi = std::exp(8.0);

    transport::basic_range<> kts(kt_lo, kt_hi, 1, transport::spacing::linear);
    transport::basic_range<> alphas(0.0, 0.0, 0, transport::spacing::linear);
    transport::basic_range<> betas(1.0 / 3.0, 0.99, 1, transport::spacing::linear);

    std::ostringstream task_name;
    task_name << "convergence.threepf-" << static_cast<unsigned int>(efolds_ics*100);

    std::ostringstream zeta_task_name;
    zeta_task_name << "convergence.threepf-zeta-" << static_cast<unsigned int>(efolds_ics*100);

    std::ostringstream task_desc;
    task_desc << "Study convergence with number of massless e-folds for different triangle shapes - " << efolds_ics << " e-folds";

    std::ostringstream zeta_task_desc;
    zeta_task_desc << "Zeta 3pf task associated with " << task_name.str() << " - " << efolds_ics << " e-folds";

    transport::threepf_alphabeta_task<> tk3(task_name.str(), ics, ts, kts, alphas, betas);
    tk3.set_adaptive_ics_efolds(efolds_ics);
    tk3.set_description(task_desc.str());

    transport::zeta_threepf_task<> ztk3(zeta_task_name.str(), tk3);
    ztk3.set_description(zeta_task_desc.str());

    return ztk3;
  }


int main(int argc, char* argv[])
  {
    // create task manager instance
    transport::task_manager<> mgr(argc, argv);

    // create model instance
    std::shared_ptr< transport::dquad_mpi<> > model = mgr.create_model< transport::dquad_mpi<> >();

    // write tasks to repository
    mgr.add_generator([=](transport::repository<>& repo) -> void { write_tasks(repo, model.get()); });

    // hand off control to task manager
    mgr.process();

    return(EXIT_SUCCESS);
  }
