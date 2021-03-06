//
//  Input.cpp
//  TRANSPORT
//
//  Created by Cheuk Lau on 1/20/14.
//  Copyright (c) 2014 Cheuk Lau. All rights reserved.
//

#include "Input.h"

Input::Input () {
    
    /* Define XML document */
    xml_document<> doc;
    
    /* Define XML pointers used */
    xml_node<> * root_node;
    xml_node<> * grid_node;
    xml_node<> * quad_node;
    xml_node<> * mat_node;
    xml_node<> * bc_node;
    
    /* Read the XML file into a vector  */
    ifstream input_file ("input_xml");
    vector<char> buffer((istreambuf_iterator<char>(input_file)), istreambuf_iterator<char>());
    buffer.push_back('\0');
    
    /* Parse the buffer using the XML file  */
    doc.parse<0> (&buffer[0]);
    
    /* Set the XML pointers */
    root_node = doc.first_node("prototype");
    grid_node = root_node->first_node("grid");
    quad_node = root_node->first_node("quad");
    mat_node = root_node->first_node("mat");
    bc_node = root_node->first_node("bc");
    
    /* Read grid values */
    // num_reg
    spat_type = grid_node->first_node("spat_type")->first_attribute("entry")->value();
    string num_reg_char = grid_node->first_node("num_reg")->first_attribute("entry")->value();
    istringstream iss_num_reg(num_reg_char); // Convert string to integer
    iss_num_reg >> num_reg;
    // num_cells
    string num_cells_char = grid_node->first_node("num_cells")->first_attribute("entry")->value();
    istringstream iss_num_cells(num_cells_char); // Convert string to vector of integers
    int conversion_num_cells;
    int counter = 0;
    num_cells.resize(num_reg);
    while (iss_num_cells >> conversion_num_cells) {
        num_cells(counter) = conversion_num_cells;
        counter++;
    }
    // reg_size
    string reg_size_char = grid_node->first_node("reg_size")->first_attribute("entry")->value();
    istringstream iss_reg_size(reg_size_char); // Convert string to vector of doubles
    double conversion_reg_size;
    counter = 0;
    reg_size.resize(num_reg);
    while (iss_reg_size >> conversion_reg_size) {
        reg_size(counter) = conversion_reg_size;
        counter++;
    }
    
    /* Read quad values */
    // quad_type
    quad_type = quad_node->first_node("type")->first_attribute("entry")->value();
    // quad_order
    string quad_order_char = quad_node->first_node("order")->first_attribute("entry")->value();
    istringstream iss_quad_order(quad_order_char); // Convert string to vector of integers
    int conversion_quad_order;
    counter = 0;
    quad_order.resize(num_reg);
    while (iss_quad_order >> conversion_quad_order) {
        quad_order(counter) = conversion_quad_order;
        counter++;
    }
    
    /* Read mat values */
    // abs_xs
    string abs_xs_char = mat_node->first_node("abs_xs")->first_attribute("entry")->value();
    istringstream iss_abs_xs(abs_xs_char); // Convert string to vector of doubles
    double conversion_abs_xs;
    counter = 0;
    abs_xs.resize(num_reg);
    while (iss_abs_xs >> conversion_abs_xs) {
        abs_xs(counter) = conversion_abs_xs;
        counter++;
    }
    // scat_xs
    string scat_xs_char = mat_node->first_node("scat_xs")->first_attribute("entry")->value();
    istringstream iss_scat_xs(scat_xs_char); // Convert string to vector of doubles
    double conversion_scat_xs;
    counter = 0;
    scat_xs.resize(num_reg);
    while (iss_scat_xs >> conversion_scat_xs) {
        scat_xs(counter) = conversion_scat_xs;
        counter++;
    }
    // ext_source
    string ext_source_char = mat_node->first_node("source")->first_attribute("entry")->value();
    istringstream iss_ext_source(ext_source_char); // Convert string to vector of doubles
    double conversion_ext_source;
    counter = 0;
    ext_source.resize(num_reg);
    while (iss_ext_source >> conversion_ext_source) {
        ext_source(counter) = conversion_ext_source;
        counter++;
    }
    
    /* Read bc values */
    // psi_left
    string psi_left_char = bc_node->first_node("left")->first_attribute("entry")->value();
    istringstream iss_psi_left(psi_left_char); // Convert string to double
    iss_psi_left >> psi_left;
    // psi_right
    string psi_right_char = bc_node->first_node("right")->first_attribute("entry")->value();
    istringstream iss_psi_right(psi_right_char); // Convert string to double
    iss_psi_right >> psi_right;

}