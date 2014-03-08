//
//  Sweep.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/17/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Sweep.h"

void LDFE_spat(Region* current_reg, vector<mat> phi, string current_dir, vec fc_source, vector<LDFE_reg*> quad_use) {

	/* Retrieve needed input parameters */
	int num_cells          = current_reg->get_num_cells();
	int scat_order         = current_reg->get_scat_order();
	double abs_xs          = current_reg->get_abs_xs();
	double reg_size        = current_reg->get_reg_size();
	double ext_source      = current_reg->get_ext_source();
	vec scat_xs            = current_reg->get_scat_xs();

	/* Calculate basic parameters */
	double tot_xs    = abs_xs+scat_xs(0);
	double mesh_size = reg_size/num_cells;
	int num_edges    = num_cells+1;

	/* Initialize needed variables */
	double tau;
	double leg;
	double scat_source_lft;
	double scat_source_rgt;
	vec dirs_use;

	/* Positive sweep */
	if (current_dir == "pos") {
		for (int i=0; i<quad_use.size(); i++) {
			dirs_use = quad_use[i]->get_dirs();
			for (int j=0; j<2; j++) {
				tau = tot_xs*mesh_size/dirs_use(j);
				for (int k=1; k<num_edges; k++) {
					scat_source_lft=0;
					scat_source_rgt=0;
					for (int r=0; r<=scat_order; r++) {
						leg = legendre(double(r), dirs_use(j));
						scat_source_lft=scat_source_lft+((2.*r+1.)/2.)*scat_xs(r)*phi[r](0,k-1)*leg;
						scat_source_rgt=scat_source_rgt+((2.*r+1.)/2.)*scat_xs(r)*phi[r](1,k-1)*leg;
					}
					quad_use[i]->psi[j](0,k)=\
							quad_use[i]->psi[j](0,k-1)*(6.-2.*tau)/(6.+4.*tau+pow(tau,2.))+\
							((scat_source_lft+ext_source+fc_source(k-1))/tot_xs)*(3.*tau/(6.+4.*tau+pow(tau,2.)))+\
							((scat_source_rgt+ext_source+fc_source(k-1))/tot_xs)*(3.*tau+pow(tau,2.))/(6.+4.*tau+pow(tau,2.));
					quad_use[i]->psi[j](1,k-1)=\
							quad_use[i]->psi[j](0,k-1)*(6.+4.*tau)/(6.+4.*tau+pow(tau,2.))+\
							((scat_source_lft+ext_source+fc_source(k-1))/tot_xs)*(tau+pow(tau,2.))/(6.+4.*tau+pow(tau,2.))-\
							((scat_source_rgt+ext_source+fc_source(k-1))/tot_xs)*(tau/(6.+4.*tau+pow(tau,2.)));
				}
			}
		}
	}

	/* Negative sweep */
	else if (current_dir == "neg") {
		for (int i=0; i<quad_use.size(); i++) {
			dirs_use = quad_use[i]->get_dirs();
			for (int j=0; j<2; j++) {
				tau = tot_xs*mesh_size/dirs_use(j);
				for (int k=num_edges-2; k>=0; k--) {
					scat_source_lft=0;
					scat_source_rgt=0;
					for (int r=0; r<=scat_order; r++) {
						leg = legendre(double(r), dirs_use(j));
						scat_source_lft=scat_source_lft+((2.*r+1.)/2.)*scat_xs(r)*phi[r](0,k)*leg;
						scat_source_rgt=scat_source_rgt+((2.*r+1.)/2.)*scat_xs(r)*phi[r](1,k)*leg;
					}
					quad_use[i]->psi[j](0,k)=\
							quad_use[i]->psi[j](0,k+1)*(6.+2.*tau)/(6.-4.*tau+pow(tau,2.))+\
							((scat_source_lft+ext_source+fc_source(k))/tot_xs)*(-3.*tau+pow(tau,2.))/(6.-4.*tau+pow(tau,2.))-\
							((scat_source_rgt+ext_source+fc_source(k))/tot_xs)*(3.*tau)/(6.-4.*tau+pow(tau,2.));
					quad_use[i]->psi[j](1,k+1)=\
							quad_use[i]->psi[j](0,k+1)*(6.-4.*tau)/(6.-4.*tau+pow(tau,2.))+\
							((scat_source_lft+ext_source+fc_source(k))/tot_xs)*(tau/(6.-4.*tau+pow(tau,2.)))+\
							((scat_source_rgt+ext_source+fc_source(k))/tot_xs)*(-tau+pow(tau,2.))/(6.-4.*tau+pow(tau,2.));
				}
			}

		}
	}
}

double legendre(double order, double ang) {

	if (order==0) {
		return 1.;
	}
	else if (order==1) {
		return ang;
	}
	else {
		return ((2.*order-1.)/order)*ang*legendre(order-1.,ang)-((order-1.)/order)*legendre(order-2.,ang);
	}

}
