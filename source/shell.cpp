

#include "shell.h"
//#include "parameters.h"
#include <vector>
#include <math.h>
#include <fstream>
#include <iostream>
using namespace std;
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_errno.h>

void shell::initialize(){
    //when point is not on the edge the edge is 0, if it is on the edge it is 1
    vertices point0(0., 0., 0., 1, 1);
    vertices point1(0.5, R3_2, 0.0, 1, 1);
    vertices point2(1., 0., 0., 1, 1);

    // vertices point0(0.2, 0., 0.4, 1, 1);
    // vertices point1(0.5, R3_2, 0.03, 1, 1);
    // vertices point2(0., 0.5, 0.3, 1, 1);

    //nt num_v = points.size();

    linevertices l0(0, 1, 0, 1000, 1, 2, 1);
    linevertices l1(1, 2, 0, 1000, 1, 0, 2);
    linevertices l2(2, 0, 0, 1000, 1, 1, 0);
    //cout<<"initialize\n"<<l0.ll<<" "<<l0.lr<<endl;

    lines.push_back(l0);
    lines.push_back(l1);
    lines.push_back(l2);
    
    cout<<lines[2].ll<<lines[0].lr<<endl;

    trianglevertices t0(0, 1, 2, 0, 1, 2, 0);

    triangles.push_back(t0);

   
    point1.neighbor_l.push_back(1);
    point2.neighbor_l.push_back(1);
    point2.neighbor_l.push_back(2);

    point0.neighbor_t.push_back(0);
    point1.neighbor_t.push_back(0);
    point2.neighbor_t.push_back(0);

    points.push_back(point0);
    points.push_back(point1);
    points.push_back(point2);

    int num_v = points.size();
    int num_l = lines.size();
    int num_t = triangles.size();
    cout<<"Yinan"<<endl;

    grow();

}
