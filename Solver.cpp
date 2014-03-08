//
//  Solver.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/17/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Solver.h"

void Solver (Input input_param, Problem problem_setup) {

	/* Retrieve input parameters */
	int max_si_cycles        = input_param.get_si_cycles();
	int max_cycles_per_adapt = input_param.get_max_cycles_per_adapt();
	int si_per_adapt         = input_param.get_si_per_adapt();
	double psi_left          = input_param.get_psi_left();
	double psi_right         = input_param.get_psi_right();
	double si_tol            = input_param.get_si_tol();
	double alpha             = input_param.get_alpha();
	double beta              = input_param.get_beta();
	double adapt_tol         = input_param.get_adapt_tol();
	double eps_thresh        = input_param.get_eps_thresh();
	Col<int> num_cells       = input_param.get_num_cells();
	Col<int> scat_order      = input_param.get_scat_order();
	vector<Region*> Regions  = problem_setup.get_Regions();
	vec scat_xs              = input_param.get_scat_xs();
	vec abs_xs               = input_param.get_abs_xs();
	vec ext_source           = input_param.get_ext_source();
	vec reg_size             = input_param.get_reg_size();
	string map_type          = input_param.get_map_type();
	string fixup_type        = input_param.get_fixup_type();
	string spat_type         = input_param.get_spat_type();
	string bc_type           = input_param.get_bc_type();
	string quad_type         = input_param.get_quad_type();
	string adapt_type        = input_param.get_adapt_type();

	/* Calculate basic parameters */
	int num_reg = int(Regions.size());

	/* Initialize parameters */
	int si_cycles;
	double si_error;
	double sigma;
	double sum_top;
	double sum_bot;
	double leg;
	double adapt_cycles;
	bool adapt_finished;
	bool exit_si;
	string dir;
	vec dirs_use;
	vec wgts_use;
	vec tot_source_lft(num_reg-1);
	vec tot_source_rgt(num_reg-1);
	vec scat_xs_reg;
	vector<vec> fc_source;
	vector<vector<mat> > phi_old;
	vector<vector<mat> > phi_old_old;
	vector<vector<mat> > phi;
	vector<vector<LDFE_reg*> > quad_pos;
	vector<vector<LDFE_reg*> > quad_neg;
	string input_pause;

	/* Phi vector */
	phi.resize(num_reg);
	for (int i=0; i<num_reg; i++) {
		phi[i].resize(scat_order(i)+1);
		for (int j=0; j<=scat_order(i); j++){
			phi[i][j] = zeros(2, num_cells(i));
		}
	}

	/* Initialize first collision source to zero */
	fc_source.resize(num_reg);
	for (int i=0; i<num_reg; i++) {
		fc_source[i] = zeros(num_cells(i));
	}

	/* Outer adaptive loop */
	adapt_cycles = 0;
	adapt_finished = false;
	do {

		cout << "starting adaptive cycle: " << adapt_cycles << endl;

		/* Gather bottom nodes of each quadrature tree*/
		quad_pos = find_quad (Regions, "pos");
		quad_neg = find_quad (Regions, "neg");

		/* Apply boundary conditions */
		apply_bc(bc_type, quad_pos, quad_neg, psi_left, psi_right, num_cells, num_reg, Regions, reg_size, fc_source);

		/* Inner SI loop */
		si_cycles = 0;
		si_error = si_tol*2.;
		exit_si = false;
		do {

			/* Fix-up terms */
			if (fixup_type=="FIXUP1") {
				fixup_source(fixup_type, scat_xs, num_reg, Regions, tot_source_lft, tot_source_rgt, phi, scat_order, num_cells, ext_source);
			}

			/* Positive sweep */
			for (int i=0; i<num_reg; i++) {
				// LDFE sweep
				if (spat_type == "LDFE") {
					LDFE_spat(Regions[i], phi[i], "pos", fc_source[i], quad_pos[i]);
				}
				else {
					cout << "Invalid spatial discretization type. Exiting program." << endl;
					exit (EXIT_FAILURE);
				}
				// Mapping and fixup
				if (i<num_reg-1) {
					if (map_type == "MAPPING1") {
						mapping1(Regions[i], Regions[i+1], "pos");
					}
					else if (map_type == "MAPPING2") {
						mapping2(Regions[i], Regions[i+1], "pos");
					}
					else {
						cout << "Invalid mapping type. Exiting program." << endl;
						exit (EXIT_FAILURE);
					}
					if (fixup_type == "FIXUP1") {
						fixup1(num_cells(i), num_cells(i+1), "pos", tot_source_lft(i), quad_pos[i], quad_pos[i+1]);
					}
				}
			}

			/* Negative sweep */
			for (int i=num_reg-1; i>=0; i--) {
				// LDFE sweep
				if (spat_type == "LDFE") {
					LDFE_spat(Regions[i], phi[i], "neg", fc_source[i], quad_neg[i]);
				}
				else {
					cout << "Invalid spatial discretization type. Exiting program." << endl;
					exit (EXIT_FAILURE);
				}
				// Mapping and fixup
				if (i>0) {
					if (map_type == "MAPPING1") {
						mapping1(Regions[i], Regions[i-1], "neg");
					}
					else if (map_type == "MAPPING2") {
						mapping2(Regions[i], Regions[i-1], "neg");
					}
					else {
						cout << "Invalid mapping type. Exiting program." << endl;
						exit (EXIT_FAILURE);
					}
					if (fixup_type == "FIXUP1") {
						fixup1(num_cells(i), num_cells(i-1), "neg", tot_source_rgt(i-1), quad_neg[i], quad_neg[i-1]);
					}
				}
			}

			/* Update moments, scattering source and SI error*/
			if (scat_xs.max()>0) {
				// Old moments
				if (si_cycles>0) {
					phi_old_old = phi_old;
				}
				phi_old = phi;
				// New moments
				for (int i=0; i<num_reg; i++) {
					for (int j=0; j<=scat_order[i]; j++) {
						phi[i][j]=zeros(2,num_cells(i));
						for (int k=0; k<num_cells(i); k++) {
							for (int r=0; r<quad_pos[i].size(); r++) {
								dirs_use = quad_pos[i][r]->get_dirs();
								wgts_use = quad_pos[i][r]->get_wgts();
								for (int s=0; s<2; s++) {
									leg=legendre(double(j),dirs_use(s));
									phi[i][j](0,k)=phi[i][j](0,k)+quad_pos[i][r]->psi[s](1,k)*leg*wgts_use(s);
									phi[i][j](1,k)=phi[i][j](1,k)+quad_pos[i][r]->psi[s](0,k+1)*leg*wgts_use(s);
								}
							}
							for (int r=0; r<quad_neg[i].size(); r++) {
								dirs_use = quad_neg[i][r]->get_dirs();
								wgts_use = quad_neg[i][r]->get_wgts();
								for (int s=0; s<2; s++) {
									leg=legendre(double(j),dirs_use(s));
									phi[i][j](0,k)=phi[i][j](0,k)+quad_neg[i][r]->psi[s](0,k)*leg*wgts_use(s);
									phi[i][j](1,k)=phi[i][j](1,k)+quad_neg[i][r]->psi[s](1,k+1)*leg*wgts_use(s);
								}
							}
						}
					}
				}
				// SI error
				if (si_cycles>0) {
					sum_top=0;
					sum_bot=0;
					for (int i=0; i<num_reg; i++) {
						for (int j=0; j<num_cells(i); j++) {
							sum_top=sum_top+pow((phi[i][0](0,j)+phi[i][0](1,j))/2-\
									(phi_old[i][0](0,j)+phi_old[i][0](1,j))/2,2);
							sum_bot=sum_bot+pow((phi_old[i][0](0,j)+phi_old[i][0](1,j))/2-\
									(phi_old_old[i][0](0,j)+phi_old_old[i][0](1,j))/2,2);
						}
					}
					sigma=sqrt(sum_top)/sqrt(sum_bot);
					si_error=sqrt(sum_bot)/(1.-sigma);
				}
				cout << "current SI cycle: " << si_cycles << " current SI error: "<< si_error << endl;
				si_cycles++;
			}

			/* Decide whether or not to leave SI loop */
			if (adapt_type == "JMP") {
				if (scat_xs.max()>0) {
					if ((si_cycles>0 && si_cycles%si_per_adapt<1e-12) || abs(si_error)<si_tol) {
						exit_si = true;
					}
					else if (si_cycles>max_si_cycles) {
						cout << "maximum number of SI cycles have been reached. Exiting program!" << endl;
						EXIT_FAILURE;
					}
				}
				else {
					exit_si = true;
				}
			}
			else {
				if (scat_xs.max()>0) {
					adapt_finished = true;
					if (abs(si_error)<si_tol){
						exit_si = true;
					}
					else if (si_cycles>max_si_cycles) {
						cout << "Maximum number of SI cycles have been reached. Exiting!" << endl;
						EXIT_FAILURE;
					}
				}
				else {
					adapt_finished = true;
					exit_si = true;
				}
			}

		} while (exit_si==false);

		/* *** TEST *** Write angular flux into HDF5 to see plots with adaptivity */
		cout << "Writing angular flux of adaptive cycle " << adapt_cycles << "to file." << endl;
		Output(problem_setup);
		cout << "Enter anything to continue." << endl;
		cin >> input_pause;
		cout << endl;
		/* *** END TEST *** */

		/* Apply adaptivity */
		if (adapt_cycles>max_cycles_per_adapt) {
			cout << "Maximum number of adaptive cycles reached. Exiting!" << endl;
			EXIT_FAILURE;
		}
		if (adapt_type == "JMP") {
			adapt(quad_pos, quad_neg, num_cells, quad_type, alpha, beta, Regions, adapt_finished, adapt_tol, eps_thresh);
			adapt_cycles++;
		}

	} while (adapt_finished==false);

}

vector<vector<LDFE_reg*> > find_quad (vector<Region*> Regions, string dir) {

	/* Initialize parameters */
	bool reached_bottom;
	LDFE_quad* quad_region;
	LDFE_reg* quad_selected;
	vector<LDFE_reg*> quad_ind;
	vector<LDFE_reg*> quad_children;
	vector<LDFE_reg*> quad_children_tmp_ind;
	vector<LDFE_reg*> quad_children_tmp_all;
	vector<vector<LDFE_reg*> > quad_return;

	/* Calculate needed parameters */
	int num_reg = Regions.size();
	int quad_ind_num;
	if (dir == "pos") {
		quad_ind_num = 0;
	}
	else {
		quad_ind_num = 1;
	}
	quad_return.resize(num_reg);

	/* Traverse quadrature tree for chosen set of directions */
	for (int z=0; z<num_reg; z++) {
		quad_region   = Regions[z]->get_quad();  // Region quadrature pointer
		quad_ind      = quad_region->get_quad(); // Vector of positive and negative LDFE regions
		quad_selected = quad_ind[quad_ind_num];  // Chosen LDFE region
		if (quad_selected->get_has_children() == true) {   // Check if the coarsest quadrature was chosen
			quad_children = quad_selected->get_children(); // Root node's first level children
			// Traverse the tree
			reached_bottom = false;
			while (reached_bottom == false) {
				// Go through the children nodes in the current level
				for (int i=0; i<quad_children.size(); i++) {
					// If the child has children add them to a temporary children vector
					if (quad_children[i]->get_has_children() == true) {
						quad_children_tmp_ind = quad_children[i]->get_children();
						quad_children_tmp_all.push_back(quad_children_tmp_ind[0]);
						quad_children_tmp_all.push_back(quad_children_tmp_ind[1]);
					}
					// If the child does not have children add pointer to LDFE regions in return vector
					else {
						quad_return[z].push_back(quad_children[i]);
					}
				}
				// If there are no children in the temporary children vector the bottom is reached
				if (quad_children_tmp_all.size() == 0) {
					reached_bottom = true;
				}
				// Set the children vector as the temporary children vector for the next level calculation
				quad_children = quad_children_tmp_all;
				quad_children_tmp_all.clear();
				quad_children_tmp_ind.clear();
			}
		}
		// Directly add LDFE regions for coarsest quadrature
		else {
			quad_return[z].push_back(quad_selected);
		}
		sort(quad_return[z].begin(), quad_return[z].end(), CompareByDirs());
	}
	return quad_return;
}
