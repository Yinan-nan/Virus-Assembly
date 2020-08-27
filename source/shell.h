//
//  Created by Yinan on 10/11/2019
//

#ifndef Shell_H
#define Shell_H


#include <vector>
#include "parameters.h"
#include "Arrow.h"
#include "Tensor.h"
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_errno.h>
using namespace std;

class vertices {

public:
    double x, y, z;
    vector<int> neighbor_l;
    vector<int> neighbor_t;
    int nn, edge;  //nn number of the triangles to this point
    vertices(double x1, double x2, double x3, int x4, int x5) : x(x1), y(x2), z(x3), nn(x4), edge(x5) {}

};

class linevertices {

public:
    int lv1, lv2, lt1, lt2, ltn, ll, lr;  // ltn: line traingle number,lr:line right NO.,ll:line left NO.
    linevertices(int v1, int v2, int t1, int t2, int ltn1, int ll1, int lr1 ) : lv1(v1), lv2(v2), lt1(t1), lt2(t2), ltn(ltn1), ll(ll1), lr(lr1) {}

};


class vectors {

public:
    double x, y, z, mag;
    vectors(double vx1, double vx2, double vx3, double mv) : x(vx1), y(vx2), z(vx3), mag(mv) {}

};

class openangles {

public:
    double th;
    int li, type; //li is the No. of the line or the vertex, type = 1 is lineangle,type = 0 is vertexangle
    openangles(int i1, double angle1, int type1 ) : li(i1), th(angle1), type(type1) {}

};
class matrix9 {

public:
    double dnx1,dny1,dnz1,dnx2,dny2,dnz2,dnx3,dny3,dnz3;

    matrix9 (double nx1, double ny1, double nz1, double nx2, double ny2, double nz2, double nx3, double ny3, double nz3) : dnx1(nx1), dny1(ny1), dnz1(nz1),dnx2(nx2),dny2(ny2),dnz2(nz2),dnx3(nx3),dny3(ny3),dnz3(nz3) {}

};



class shell{
public:

    double th0 = 2 * asin(1 / sqrt((12 * r0 * r0) - 3));

    vertices *v;
    linevertices *l;
    trianglevertices *t;
    int num_t;
    int num_v;
    int num_l;
    vector<vertices> points;
    vector<linevertices> lines;
    vector<trianglevertices> triangles;
    vector<vertices> pointsbackup;
    vector<linevertices> linesbackup;
    vector<trianglevertices> trianglesbackup;
    gsl_vector *pointlist;
    gsl_vector *dElist;
    int number();
    int chooseangle();

    vector<vectors> lvector;
    vector<vectors> nvector;
    vector<vectors> nhat;
    vector<Arrow> larrow;
    vector<Arrow> rvlarrow;
    vector<Arrow> narrow;
    vector<Arrow> nhatarrow;

    vector<openangles> theta;
    vector<matrix9> dnhat;
    vector<denergy> dEs;
    vector<denergy> dEb1;
    vector<denergy> dEb2;
    vector<denergy> dEt;
    vector<vertices> oldpoints;
    vector<denergy> dEtold;

    void initialize();
    void linevec();
    void normalvec();
    int openangle();
    double lineangle(int i, int j, int sign, int h);
    int nextmove();
    void nextvertex(int ln);
    void grow();
    void insert_tri(int vertex);
    void merge_lines(int vertex);
    double Energy(const gsl_vector *pointlist);
    void update(const gsl_vector *pointlist);
    double relax();
    void backup();
    void restore();
    int check_closeshell(int v,int w);
    void outputfile(int step);

    inline Arrow dmag_darrow(Arrow a);
    inline Arrow ddot_darrow(Arrow &a,Arrow &b);
    inline Tensor dcross_darrow(Arrow &a,Arrow &b);
    inline Tensor dhat_darrow(Arrow &a);
    inline Tensor dnorm_dvi(const int &t, const int &vi);
    inline Arrow dhs_dvi(const int &l, const int  &vi);
    inline Arrow dhb_dvi(const int &tm, const int &tn, Arrow &n_m, Arrow &n_n, const int &vi, const int &sin_sign);

    //-----------------------------------------------------------
    static double Energy_wrapper(const gsl_vector *pointlist, void *params){
        return static_cast<shell*>(params)->Energy(pointlist);
    }
    void gradE (const gsl_vector *pointlist, gsl_vector *dElist);
    static void gradE_wrapper(const gsl_vector *pointlist, void *params, gsl_vector *dElist){
        return static_cast<shell*>(params)->gradE(pointlist,dElist);
    }
    void EdE (const gsl_vector *pointlist, double *h, gsl_vector *dElist);
    static void EdE_wrapper(const gsl_vector *pointlist,void *params, double *h,
                            gsl_vector *dElist){
        return static_cast<shell*>(params)->EdE(pointlist,h,dElist);
    }

};

#endif