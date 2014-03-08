/*
 * Adapt.cpp
 *
 *  Created on: Feb 28, 2014
 *      Author: cheuklau
 */

#include "adapt.h"

void adapt (vector<vector<LDFE_reg*> > quad_pos, vector<vector<LDFE_reg*> > quad_neg, Col<int> num_cells, string quad_type, double alpha, double beta, vector<Region*> Regions, bool& adapt_finished, double adapt_tol, double eps_thresh) {

	/* Initialize parameters */
	int num_reg = quad_pos.size();
	int num_LDFE_reg;
	double err_global;
	double mu_min;
	double mu_center;
	double mu_max;
	vec dirs_tmp;
	vector<vector<vec> > err_ind_pos(num_reg);
	vector<vector<vec> > err_ind_neg(num_reg);
	vector<vector<vector<bool> > > refine_pos(num_reg);
	vector<vector<vector<bool> > > coarsen_pos(num_reg);
	vector<vector<vector<bool> > > refine_neg(num_reg);
	vector<vector<vector<bool> > > coarsen_neg(num_reg);
	vector<vector<LDFE_reg*> > quad_test;
	string input_pause;

	/* Calculate error indicators */
	for (int i=0; i<num_reg; i++) {
		err_ind_pos[i].resize(2);
		err_ind_neg[i].resize(2);
		for (int j=0; j<2; j++) {
			err_ind_pos[i][j].resize(quad_pos[i].size());
			err_ind_neg[i][j].resize(quad_neg[i].size());
		}
		apply_adapt_ind(quad_pos[i], err_ind_pos[i], num_cells(i), quad_type, eps_thresh);
		apply_adapt_ind(quad_neg[i], err_ind_neg[i], num_cells(i), quad_type, eps_thresh);
	}

	/* Calculate global error */
	err_global = 0;
	num_LDFE_reg = 0;
	for (int i=0; i<num_reg; i++) {
		for (int j=0; j<2; j++) {
			for (int k=0; k<quad_pos[i].size(); k++){
				err_global = err_global+pow(err_ind_pos[i][j](k),2);
				num_LDFE_reg++;
			}
			for (int k=0; k<quad_neg[i].size(); k++){
				err_global = err_global+pow(err_ind_neg[i][j](k),2);
				num_LDFE_reg++;
			}
		}
	}
	err_global = pow(err_global/num_LDFE_reg, 0.5);
	cout << "Global error is: " << err_global << endl;

	/* Check if adaptive tolerance is satisfied */
	if (err_global < adapt_tol) {
		adapt_finished = true;
	}

	/* Resize refine and coarsen vectors */
	else {
		// Flag each LDFE regions for refinement or coarsenment
		for (int i=0; i<num_reg; i++) {
			refine_pos[i].resize(2);
			refine_neg[i].resize(2);
			coarsen_pos[i].resize(2);
			coarsen_neg[i].resize(2);
			for (int j=0; j<2; j++) {
				refine_pos[i][j].resize(quad_pos[i].size());
				refine_neg[i][j].resize(quad_neg[i].size());
				coarsen_pos[i][j].resize(quad_pos[i].size());
				coarsen_neg[i][j].resize(quad_neg[i].size());
				for (int k=0; k<quad_pos[i].size(); k++) {
					if (err_ind_pos[i][j](k)>alpha*err_global) {
						refine_pos[i][j][k] = true;
						coarsen_pos[i][j][k] = false;
					}
					else if (err_ind_pos[i][j](k)<beta*err_global) {
						refine_pos[i][j][k] = false;
						coarsen_pos[i][j][k]= false;
					}
					else {
						refine_pos[i][j][k] = false;
						coarsen_pos[i][j][k]= false;
					}
				}
				for (int k=0; k<quad_neg[i].size(); k++) {
					if (err_ind_neg[i][j](k)>alpha*err_global) {
						refine_neg[i][j][k] = true;
						coarsen_neg[i][j][k] = false;
					}
					else if (err_ind_neg[i][j](k)<beta*err_global) {
						refine_neg[i][j][k] = false;
						coarsen_neg[i][j][k]= false;
					}
					else {
						refine_neg[i][j][k] = false;
						coarsen_neg[i][j][k]= false;
					}
				}
			}
		}

		// Refine positive LDFE regions
		cout << "refining positive regions" << endl;
		for (int i=0; i<num_reg; i++) {
			for (int j=0; j<quad_pos[i].size(); j++)
			{
				if (refine_pos[i][0][j] == true || refine_pos[i][1][j] == true) {
					quad_pos[i][j]->has_children = true;
					dirs_tmp = quad_pos[i][j]->get_dirs();
					mu_min    = dirs_tmp(0)-(dirs_tmp(1)-dirs_tmp(0))/2.;
					mu_center = dirs_tmp(0)+(dirs_tmp(1)-dirs_tmp(0))/2.;
					mu_max    = dirs_tmp(1)+(dirs_tmp(1)-dirs_tmp(0))/2.;
					quad_pos[i][j]->children[0] = new LDFE_reg(mu_min, mu_center, 0, 0, quad_type); // Left positive LDFE regions
					quad_pos[i][j]->children[1] = new LDFE_reg(mu_center, mu_max, 0, 0, quad_type); // Right positive LDFE regions
				}
			}
		}

		// Refine negative LDFE regions
		cout << "refining negative regions" << endl;
		for (int i=0; i<num_reg; i++) {
			for (int j=0; j<quad_neg[i].size(); j++)
			{
				if (refine_neg[i][0][j] == true || refine_neg[i][1][j] == true) {
					quad_neg[i][j]->has_children = true;
					dirs_tmp = quad_neg[i][j]->get_dirs();
					mu_min    = dirs_tmp(0)-(dirs_tmp(1)-dirs_tmp(0))/2.;
					mu_center = dirs_tmp(0)+(dirs_tmp(1)-dirs_tmp(0))/2.;
					mu_max    = dirs_tmp(1)+(dirs_tmp(1)-dirs_tmp(0))/2.;
					quad_neg[i][j]->children[1] = new LDFE_reg(mu_min, mu_center, 0, 0, quad_type); // Positive LDFE regions
					quad_neg[i][j]->children[0] = new LDFE_reg(mu_center, mu_max, 0, 0, quad_type); // Negative LDFE regions
				}
			}
		}

		/* Coarsen positive LDFE regions
			cout << "coarsening positive regions" << endl;
			int counter = 0;
			do {
				if (coarsen_pos[i][counter][0] == true && coarsen_pos[i][counter][1] == true) {
					// Check if adjacent LDFE range also needs to be coarsened
					if (coarsen_pos[i][counter+1][0] == true && coarsen_pos[i][counter][1] == true) {
						// Check if they are the same refinement level
						vec dirs_tmp_1 = quad_pos[i][counter]->get_dirs();
						double dist_tmp_1 = dirs_tmp_1(1) - dirs_tmp_1(0);
						vec dirs_tmp_2 = quad_pos[i][counter+1]->get_dirs();
						double dist_tmp_2 = dirs_tmp_2(1) - dirs_tmp_2(0);
						if (abs(dist_tmp_1 - dist_tmp_2)<1e-12) {
							// Find parent and delete the children nodes
							coarsen(Regions, "pos");
							counter = counter+2;
						}
						else {
							counter++;
						}
					}
					else {
						counter++;
					}
				}
				else {
					counter++;
				}
			} while (counter < quad_pos[i].size());

			// Coarsen negative LDFE regions
			cout << "coarsening negative regions" << endl;
			counter = 0;
			do {
				if (coarsen_neg[i][counter][0] == true && coarsen_neg[i][counter][1] == true) {
					// Check if adjacent LDFE range also needs to be coarsened
					if (coarsen_neg[i][counter+1][0] == true && coarsen_neg[i][counter][1] == true) {
						// Check if they are the same refinement level
						vec dirs_tmp_1 = quad_neg[i][counter]->get_dirs();
						double dist_tmp_1 = dirs_tmp_1(1) - dirs_tmp_1(0);
						vec dirs_tmp_2 = quad_neg[i][counter+1]->get_dirs();
						double dist_tmp_2 = dirs_tmp_2(1) - dirs_tmp_2(0);
						if (abs(dist_tmp_1 - dist_tmp_2)<1e-12) {
							// Coarsen it!
							counter = counter+2;
						}
						else {
							counter++;
						}
					}
					else {
						counter++;
					}
				}
				else {
					counter++;
				}

			} while (counter < quad_neg[i].size());
		 */

	}

}

void apply_adapt_ind (vector<LDFE_reg*> quad_use, vector<vec>& err_ind, int num_cells, string quad_type, double eps_thresh) {

	/* Initialize parameters */
	Col<int> spat_pos(2);
	double left_edge;
	double right_edge;
	double norm_ind;
	double norm_half_space;
	vec dirs;
	mat basis;
	mat proj_store;

	/* Spatial indices */
	spat_pos(0) = 0;
	spat_pos(1) = num_cells;

	/* Calculate error indicators for each spatial side */
	for (int r=0; r<2; r++) {

		/* Projected solutions for each LDFE region */
		proj_store.resize(2,quad_use.size());
		for (int k=0; k<quad_use.size(); k++) {
			basis = quad_use[k]->get_basis();
			dirs = quad_use[k]->get_dirs();
			if (quad_type == "LDFE_MU") {
				left_edge = dirs(0)-(dirs(1)-dirs(0))/2.;
				right_edge = dirs(1)+(dirs(1)-dirs(0))/2.;
			}
			else if (quad_type == "LDFE_THETA") {
				left_edge = cos(acos(dirs(0))-(acos(dirs(1))-acos(dirs(0)))/2.);
				right_edge = cos(acos(dirs(1))+(acos(dirs(1))-acos(dirs(0)))/2.);
			}
			else {
				cout << "Error with quadrature definition." << endl;
				EXIT_FAILURE;
			}
			proj_store(0,k) = quad_use[k]->psi[0](0,spat_pos(r))*(basis(0,0)+basis(1,0)*left_edge)+\
					quad_use[k]->psi[1](0,spat_pos(r))*(basis(0,1)+basis(1,1)*left_edge);
			proj_store(1,k) = quad_use[k]->psi[0](0,spat_pos(r))*(basis(0,0)+basis(1,0)*right_edge)+\
					quad_use[k]->psi[1](0,spat_pos(r))*(basis(0,1)+basis(1,1)*right_edge);
		}

		/* Norm of the half-space */
		norm_half_space = 0;
		for (int k=0; k<quad_use.size(); k++) {
			norm_half_space = norm_half_space+pow(quad_use[k]->psi[0](0,spat_pos(r)),2)+\
					pow(quad_use[k]->psi[1](0,spat_pos(r)),2);
		}
		norm_half_space = pow(norm_half_space/double(quad_use.size()), 0.5);

		/* Error indicator for each LDFE range*/
		if (dirs(0)>0) {
			norm_ind = pow((pow(quad_use[0]->psi[0](0,spat_pos(r)),2)+\
					pow(quad_use[0]->psi[1](0,spat_pos(r)),2))/2.,0.5);
			err_ind[r](0) = abs(proj_store(1,0)-proj_store(0,1))/(2.*(norm_ind+eps_thresh*norm_half_space));
			for (int k=1; k<quad_use.size()-1; k++) {
				norm_ind = pow((pow(quad_use[k]->psi[0](0,spat_pos(r)),2)+\
						pow(quad_use[k]->psi[1](0,spat_pos(r)),2))/2.,0.5);
				err_ind[r](k) = (abs(proj_store(1,k-1)-proj_store(0,k))+\
						abs(proj_store(1,k)-proj_store(0,k+1)))/(2.*(norm_ind+eps_thresh*norm_half_space));
			}
			norm_ind = pow((pow(quad_use[quad_use.size()-1]->psi[0](0,spat_pos(r)),2)+\
					pow(quad_use[quad_use.size()-1]->psi[1](0,spat_pos(r)),2))/2., 0.5);
			err_ind[r](quad_use.size()-1) = abs(proj_store(1,quad_use.size()-2)-proj_store(0,quad_use.size()-1))/(2.*(norm_ind+eps_thresh*norm_half_space));
		}
		else {
			norm_ind = pow((pow(quad_use[0]->psi[0](0,spat_pos(r)),2)+\
					pow(quad_use[0]->psi[1](0,spat_pos(r)),2))/2.,0.5);
			err_ind[r](0) = abs(proj_store(1,1)-proj_store(0,0))/(2.*(norm_ind+eps_thresh*norm_half_space));
			for (int k=1; k<quad_use.size()-1; k++) {
				norm_ind = pow((pow(quad_use[k]->psi[0](0,spat_pos(r)),2)+\
						pow(quad_use[k]->psi[1](0,spat_pos(r)),2))/2.,0.5);
				err_ind[r](k) = (abs(proj_store(1,k+1)-proj_store(0,k))+\
						abs(proj_store(1,k)-proj_store(0,k-1)))/(2.*(norm_ind+eps_thresh*norm_half_space));
			}
			norm_ind = pow((pow(quad_use[quad_use.size()-1]->psi[0](0,spat_pos(r)),2)+\
					pow(quad_use[quad_use.size()-1]->psi[1](0,spat_pos(r)),2))/2.,0.5);
			err_ind[r](quad_use.size()-1) = abs(proj_store(1,quad_use.size()-1)-proj_store(0,quad_use.size()-2))/(2.*(norm_ind+eps_thresh*norm_half_space));
		}

	}
}

void coarsen (vector<Region*> Regions, string dir) {
	cout << "coarsening called!" << endl;
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
		quad_region   = Regions[z]->get_quad();  // Pointer region quadrature
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
					// If the child does not have children check if this is the parent, if so coarsen
					else {
						quad_children[i]->has_children = false;
						quad_children[i]->children.clear();
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
		// Coarsen for coarsest quadrature
		else {
			quad_selected->has_children = false;
			quad_selected->children.clear();
		}
	}
}
