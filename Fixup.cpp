//
//  Fixup.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 2/14/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Fixup.h"

void fixup1 (int num_cells_1, int num_cells_2, string dir, double tot_source, vector<LDFE_reg*> quad_inc, vector<LDFE_reg*> quad_out) {

	/* Initialize parameters */
	int spat_pos_1;
	int spat_pos_2;
	int counter;
	double left_bound;
	double left_dir;
	double right_bound;
	double right_dir;
	double m;
	double b;
	vec dir_tmp;
	vec dir_tmp_order;
	uvec sorted_index;
	mat row_tmp_1(1,2);
	mat row_tmp_2(1,2);
	mat interp_table;
	mat interp_table_sorted;
	vector<mat> psi_tmp;

	/* Spatial edge index */
	if (dir == "pos") {
		spat_pos_1 = num_cells_1;
		spat_pos_2 = 0;
	}
	else {
		spat_pos_1 = 0;
		spat_pos_2 = num_cells_2;
	}

	/* Create interpolation table from quadrature you are mapping from */
	counter = 0;
	for (int i=0; i<quad_inc.size(); i++) {
		psi_tmp = quad_inc[i]->psi;
		dir_tmp = quad_inc[i]->get_dirs();
		row_tmp_1(0,0) = dir_tmp(0);
		row_tmp_1(0,1) = psi_tmp[0](0,spat_pos_1);
		row_tmp_2(0,0) = dir_tmp(1);
		row_tmp_2(0,1) = psi_tmp[1](0,spat_pos_1);
		interp_table.insert_rows(counter, row_tmp_1);
		interp_table.insert_rows(counter+1, row_tmp_2);
		counter = counter+2;
	}

	/* Reorder interpolation table */
	dir_tmp_order.resize(interp_table.n_rows);
	for (int i=0; i<interp_table.n_rows; i++) {
		dir_tmp_order(i) = interp_table(i,0);
	}
	sorted_index = sort_index(dir_tmp_order);
	interp_table_sorted.resize(interp_table.n_rows, 2);
	for (int i=0; i<interp_table.n_rows; i++) {
		interp_table_sorted(i,0) = interp_table(sorted_index(i),0);
		interp_table_sorted(i,1) = interp_table(sorted_index(i),1);
	}

	/* Fixup mapped solution */
	for (int i=0; i<quad_out.size(); i++) {
		psi_tmp = quad_out[i]->psi;
		dir_tmp = quad_out[i]->get_dirs();
		// Go through each direction
		for (int j=0; j<2; j++){
			// Check if direction is between the leftmost edge and the first direction
			if (dir_tmp(j) < interp_table_sorted(0,0)) {
				// For positive directions, leftmost edge is at mu = 0
				if (dir == "pos") {
					left_bound = tot_source;
					left_dir = 0;
					right_bound = interp_table_sorted(0,1);
					right_dir = interp_table_sorted(0,0);
					// Check to see if mapped solution poses an extrema
					if ((psi_tmp[j](0,spat_pos_2) < right_bound && psi_tmp[j](0,spat_pos_2) < left_bound) || \
							(psi_tmp[j](0,spat_pos_2) > right_bound && psi_tmp[j](0,spat_pos_2) > left_bound)) {
						m = (right_bound-left_bound)/(right_dir-left_dir);
						b = (right_dir*left_bound-left_dir*right_bound)/(right_dir-left_dir);
						quad_out[i]->psi[j](0,spat_pos_2) = m*dir_tmp(j)+b;
					}
				}
				// For negative directions, leftmost edge is at mu = -1
				else {
					quad_out[i]->psi[j](0,spat_pos_2) = interp_table_sorted(0,1);
				}
			}
			// Check if direction is between the rightmost edge and the last direction
			else if (dir_tmp(j) > interp_table_sorted(interp_table_sorted.n_rows-1, 0)) {
				// For positive directions, rightmost edge is at mu = 1
				if (dir == "pos") {
					quad_out[i]->psi[j](0,spat_pos_2) = interp_table_sorted(interp_table_sorted.n_rows-1,1);
				}
				// For negative directions, rightmost edge is at mu = 0
				else {
					right_bound = tot_source;
					right_dir = 0;
					left_bound = interp_table_sorted(interp_table_sorted.n_rows-1,1);
					left_dir = interp_table_sorted(interp_table_sorted.n_rows-1,0);
					// Check to see if mapped solution poses an extrema
					if ((psi_tmp[j](0,spat_pos_2) < right_bound && psi_tmp[j](0,spat_pos_2) < left_bound) || \
							(psi_tmp[j](0,spat_pos_2) > right_bound && psi_tmp[j](0,spat_pos_2) > left_bound)) {
						m = (right_bound-left_bound)/(right_dir-left_dir);
						b = (right_dir*left_bound-left_dir*right_bound)/(right_dir-left_dir);
						quad_out[i]->psi[j](0,spat_pos_2) = m*dir_tmp(j)+b;
					}
				}
			}
			// Otherwise direction is between two interior directions
			else {
				// Go through each row of the sorted interpolation table
				for (int r=1; r<interp_table_sorted.n_rows; r++) {
					// Check to see if the correct bounds are found
					if (dir_tmp(j) < interp_table_sorted(r,0) && dir_tmp(j) > interp_table_sorted(r-1,0)) {
						left_dir    = interp_table_sorted(r-1,0);
						right_dir   = interp_table_sorted(r,0);
						left_bound  = interp_table_sorted(r-1,1);
						right_bound = interp_table_sorted(r,1);
						// Check to see if mapped solution poses an extrema
						if ((psi_tmp[j](0,spat_pos_2) < right_bound && psi_tmp[j](0,spat_pos_2) < left_bound) || \
								(psi_tmp[j](0,spat_pos_2) > right_bound && psi_tmp[j](0,spat_pos_2) > left_bound)) {
							m = (right_bound - left_bound) / (right_dir - left_dir);
							b = (right_dir*left_bound-left_dir*right_bound)/(right_dir-left_dir);
							quad_out[i]->psi[j](0,spat_pos_2) = m*dir_tmp(j)+b;
							r=interp_table_sorted.n_rows;
						}
						// Mapped solution does not need fixup
						else {
							r=interp_table_sorted.n_rows;
						}
					}
				}
			}
		}
	}
}

void fixup_source(string fixup_type, vec scat_xs, int num_reg, vector<Region*> Regions, vec& tot_source_lft, vec& tot_source_rgt, vector<vector<mat> > phi, Col<int> scat_order, Col<int> num_cells, vec ext_source) {

	/* Define parameters */
	vec scat_xs_reg;
	double abs_xs_reg;
	double leg;

	/* Calculate left and right total sources at spatial interfaces for fix-up */
	if (fixup_type == "FIXUP1" && scat_xs.max() > 0.) {
		// Go through each interface
		for (int i=0; i<num_reg-1; i++) {
			// Calculate left scattering source
			tot_source_lft(i) = 0;
			scat_xs_reg = Regions[i]->get_scat_xs();
			abs_xs_reg = Regions[i]->get_abs_xs();
			// Add the contributions from each moment
			for (int s=0; s<=scat_order[i]; s++) {
				leg = legendre(double(s), 0.);
				tot_source_lft(i)=tot_source_lft(i)+((2.*s+1.)/2.)*scat_xs_reg(s)*phi[i][s](0,num_cells(i)-1)*leg;
			}
			// Add the external source
			tot_source_lft(i) = (tot_source_lft(i)+ext_source(i))/(abs_xs_reg+scat_xs_reg(0));
			// Calculate right scattering source
			tot_source_rgt(i)=0;
			scat_xs_reg = Regions[i+1]->get_scat_xs();
			abs_xs_reg = Regions[i+1]->get_abs_xs();
			// Add the contributions from each moment
			for (int s=0; s<=scat_order[i+1]; s++) {
				leg = legendre(double(s), 0.);
				tot_source_rgt(i)=tot_source_rgt(i)+((2.*s+1.)/2.)*scat_xs_reg(s)*phi[i+1][s](0,0)*leg;
			}
			// Add the external source
			tot_source_rgt(i) = (tot_source_rgt(i)+ext_source(i+1))/(abs_xs_reg+scat_xs_reg(0));
		}
	}
	else if (fixup_type == "FIXUP1") {
		for (int i=0; i<num_reg-1; i++) {
			tot_source_lft(i) = 0.;
			tot_source_rgt(i) = 0.;
		}
	}
}
