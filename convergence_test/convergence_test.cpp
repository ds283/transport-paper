//
// Created by David Seery on 24/06/2016.
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

    const double phi_init = 10.0 * Mp;
    const double chi_init = 12.9 * Mp;

    const double N_init = 0.0;
    const double N_pre = 17.0;
    const double N_end = 60.0;

    transport::parameters<> params(Mp, { Mphi, Mchi }, model);
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

    vis_toolkit::SQL_time_query last_time("serial IN (SELECT MAX(serial) FROM time_samples)");
    vis_toolkit::SQL_threepf_query all_threepfs("1=1");

    vis_toolkit::zeta_threepf_wavenumber_series<> z300(ztk300, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost300(*dynamic_cast< transport::threepf_task<>* >(ztk300.get_parent_task()), all_threepfs);
    z300.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("300", "300");
    cost300.set_label_text("300", "300");

    vis_toolkit::zeta_threepf_wavenumber_series<> z325(ztk325, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost325(*dynamic_cast< transport::threepf_task<>* >(ztk325.get_parent_task()), all_threepfs);
    z325.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("325", "325");
    cost325.set_label_text("325", "325");

    vis_toolkit::zeta_threepf_wavenumber_series<> z350(ztk350, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost350(*dynamic_cast< transport::threepf_task<>* >(ztk350.get_parent_task()), all_threepfs);
    z350.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("350", "350");
    cost350.set_label_text("350", "350");

    vis_toolkit::zeta_threepf_wavenumber_series<> z375(ztk375, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost375(*dynamic_cast< transport::threepf_task<>* >(ztk375.get_parent_task()), all_threepfs);
    z375.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("375", "375");
    cost375.set_label_text("375", "375");

    vis_toolkit::zeta_threepf_wavenumber_series<> z400(ztk400, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost400(*dynamic_cast< transport::threepf_task<>* >(ztk400.get_parent_task()), all_threepfs);
    z400.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("400", "400");
    cost400.set_label_text("400", "400");

    vis_toolkit::zeta_threepf_wavenumber_series<> z425(ztk425, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost425(*dynamic_cast< transport::threepf_task<>* >(ztk425.get_parent_task()), all_threepfs);
    z425.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("425", "425");
    cost425.set_label_text("425", "425");

    vis_toolkit::zeta_threepf_wavenumber_series<> z450(ztk450, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost450(*dynamic_cast< transport::threepf_task<>* >(ztk450.get_parent_task()), all_threepfs);
    z450.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("450", "450");
    cost450.set_label_text("450", "450");

    vis_toolkit::zeta_threepf_wavenumber_series<> z475(ztk475, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost475(*dynamic_cast< transport::threepf_task<>* >(ztk475.get_parent_task()), all_threepfs);
    z475.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("475", "475");
    cost475.set_label_text("474", "474");

    vis_toolkit::zeta_threepf_wavenumber_series<> z500(ztk500, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost500(*dynamic_cast< transport::threepf_task<>* >(ztk500.get_parent_task()), all_threepfs);
    z500.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("500", "500");
    cost500.set_label_text("500", "500");

    vis_toolkit::zeta_threepf_wavenumber_series<> z525(ztk525, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost525(*dynamic_cast< transport::threepf_task<>* >(ztk525.get_parent_task()), all_threepfs);
    z525.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("525", "525");
    cost525.set_label_text("525", "525");

    vis_toolkit::zeta_threepf_wavenumber_series<> z550(ztk550, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost550(*dynamic_cast< transport::threepf_task<>* >(ztk550.get_parent_task()), all_threepfs);
    z550.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("550", "550");
    cost550.set_label_text("550", "550");

    vis_toolkit::zeta_threepf_wavenumber_series<> z575(ztk575, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost575(*dynamic_cast< transport::threepf_task<>* >(ztk575.get_parent_task()), all_threepfs);
    z575.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("575", "575");
    cost575.set_label_text("575", "575");

    vis_toolkit::zeta_threepf_wavenumber_series<> z600(ztk600, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost600(*dynamic_cast< transport::threepf_task<>* >(ztk600.get_parent_task()), all_threepfs);
    z600.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("600", "600");
    cost600.set_label_text("600", "600");

    vis_toolkit::zeta_threepf_wavenumber_series<> z625(ztk625, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost625(*dynamic_cast< transport::threepf_task<>* >(ztk625.get_parent_task()), all_threepfs);
    z625.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("625", "625");
    cost625.set_label_text("625", "625");

    vis_toolkit::zeta_threepf_wavenumber_series<> z650(ztk650, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost650(*dynamic_cast< transport::threepf_task<>* >(ztk650.get_parent_task()), all_threepfs);
    z650.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("650", "650");
    cost650.set_label_text("650", "650");

    vis_toolkit::zeta_threepf_wavenumber_series<> z675(ztk675, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost675(*dynamic_cast< transport::threepf_task<>* >(ztk675.get_parent_task()), all_threepfs);
    z675.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("675", "675");
    cost675.set_label_text("675", "675");

    vis_toolkit::zeta_threepf_wavenumber_series<> z700(ztk700, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost700(*dynamic_cast< transport::threepf_task<>* >(ztk700.get_parent_task()), all_threepfs);
    z700.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("700", "700");
    cost700.set_label_text("700", "700");

    vis_toolkit::zeta_threepf_wavenumber_series<> z725(ztk725, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost725(*dynamic_cast< transport::threepf_task<>* >(ztk725.get_parent_task()), all_threepfs);
    z725.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("725", "725");
    cost725.set_label_text("725", "725");

    vis_toolkit::zeta_threepf_wavenumber_series<> z750(ztk750, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost750(*dynamic_cast< transport::threepf_task<>* >(ztk750.get_parent_task()), all_threepfs);
    z750.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("750", "750");
    cost750.set_label_text("750", "750");

    vis_toolkit::zeta_threepf_wavenumber_series<> z775(ztk775, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost775(*dynamic_cast< transport::threepf_task<>* >(ztk775.get_parent_task()), all_threepfs);
    z775.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("775", "775");
    cost775.set_label_text("775", "775");

    vis_toolkit::zeta_threepf_wavenumber_series<> z800(ztk800, last_time, all_threepfs);
    vis_toolkit::cost_wavenumber<> cost800(*dynamic_cast< transport::threepf_task<>* >(ztk800.get_parent_task()), all_threepfs);
    z800.set_use_kt_label(true).set_use_alpha_label(true).set_use_beta_label(true).set_label_text("800", "800");
    cost800.set_label_text("800", "800");

    vis_toolkit::wavenumber_series_table<> convergence_table("convergence_table", "convergence_table.csv");
    convergence_table += z300 + z325 + z350 + z375
      + z400 + z425 + z450 + z475
      + z500 + z525 + z550 + z575
      + z600 + z625 + z650 + z675
      + z700 + z725 + z750 + z775
      + z800;

    vis_toolkit::wavenumber_series_table<> timing_table("timing_table", "timing_table.csv");
    timing_table += cost300 + cost325 + cost350 + cost375
      + cost400 + cost425 + cost450 + cost475
      + cost500 + cost525 + cost550 + cost575
      + cost600 + cost625 + cost650 + cost675
      + cost700 + cost725 + cost750 + cost775
      + cost800;

    transport::output_task<> out_tk("convergence.output");
    out_tk += convergence_table + timing_table;

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
    tk3.set_adaptive_ics_efolds(efolds_ics).set_description(task_desc.str());

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
