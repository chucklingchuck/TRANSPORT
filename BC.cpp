/*
 * BC.cpp
 *
 *  Created on: Feb 27, 2014
 *      Author: cheuklau
 */

#include "bc.h"

void apply_bc(string bc_type, vector<vector<LDFE_reg*> > quad_pos, vector<vector<LDFE_reg*> > quad_neg, double psi_left, double psi_right, Col<int> num_cells, int num_reg, vector<Region*> Regions, vec reg_size, vector<vec>& fc_source) {

	/* Define parameters */
	vec wgts_use;
	vec dirs_use;

	/* Apply isotropic boundary condition */
	if (bc_type == "ISOTROPIC") {
		// Left boundary condition
		double sum = 0;
		for (int i=0; i<quad_pos[0].size(); i++) {
			wgts_use = quad_pos[0][i]->get_wgts();
			dirs_use = quad_pos[0][i]->get_dirs();
			for (int j=0; j<2; j++) {
				sum = sum+wgts_use(j)*dirs_use(j);
			}
		}
		for (int i=0; i<quad_pos[0].size(); i++) {
			quad_pos[0][i]->psi.resize(2);
			for (int j=0; j<2; j++) {
				quad_pos[0][i]->psi[j].resize(2, num_cells(0)+1);
				quad_pos[0][i]->psi[j](0,0)=psi_left/sum;
			}
		}
		// Right boundary condition
		sum = 0;
		for (int i=0; i<quad_neg[0].size(); i++) {
			wgts_use = quad_neg[0][i]->get_wgts();
			dirs_use = quad_neg[0][i]->get_dirs();
			for (int j=0; j<2; j++) {
				sum = sum+wgts_use(j)*abs(dirs_use(j));
			}
		}
		for (int i=0; i<quad_neg[num_reg-1].size(); i++) {
			quad_neg[num_reg-1][i]->psi.resize(2);
			for (int j=0; j<2; j++) {
				quad_neg[num_reg-1][i]->psi[j].resize(2, num_cells[num_reg-1]+1);
				quad_neg[num_reg-1][i]->psi[j](0,num_cells[num_reg-1])=psi_right/sum;
			}
		}
	}

	/* Apply beam boundary condition */
	else if (bc_type == "BEAM") {
		// First collision source
		double psi_last = psi_left;
		for (int i=0; i<num_reg; i++) {
			vec scat_xs_reg = Regions[i]->get_scat_xs();
			double abs_xs = Regions[i]->get_abs_xs();
			double tot_xs = scat_xs_reg(0)+abs_xs;
			double cell_size = reg_size(i)/num_cells(i);
			vec fc_analytic;
			fc_analytic.resize(num_cells(i)+1);
			for (int j=0; j<num_cells(i)+1; j++) {
				fc_analytic(j) = scat_xs_reg(0)*psi_last*exp(-1.*tot_xs*j*cell_size);
			}
			for (int j=0; j<num_cells(i); j++){
				fc_source[i](j) = ((fc_analytic(j)+fc_analytic(j+1))/2.)*cell_size;
			}
			psi_last = fc_analytic(num_cells(i));
		}
		// Vacuum left boundary condition
		for (int i=0; i<quad_pos[0].size(); i++) {
			quad_pos[0][i]->psi.resize(2);
			for (int j=0; j<2; j++) {
				quad_pos[0][i]->psi[j].resize(2,num_cells[0]+1);
				quad_pos[0][i]->psi[j](0,0)=0;
			}
		}
		// Vacuum right boundary condition
		for (int i=0; i<quad_neg[num_reg-1].size(); i++) {
			quad_neg[num_reg-1][i]->psi.resize(2);
			for (int j=0; j<2; j++) {
				quad_neg[num_reg-1][i]->psi[j].resize(2,num_cells[num_reg-1]+1);
				quad_neg[num_reg-1][i]->psi[j](0,num_cells[num_reg-1])=0;
			}
		}
	}

	/* Invalid boundary condition definition */
	else {
		cout << "Invalid boundary condition type. Exiting program." << endl;
		exit (EXIT_FAILURE);
	}
}
