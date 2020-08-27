//
// Created by Yinan on 2/20/2020.
//
#include <iostream>
#include "shell.h"
#include <fstream>
#include <math.h>
#include "parameters.h"
#include <time.h>
// #include <bits/stdc++.h> 
using namespace std;


void shell::update(const gsl_vector *pointlist) {
	num_v = points.size();
	for (int i = 0; i < num_v; i++) {
		points[i].x = gsl_vector_get(pointlist,3*i);
		points[i].y = gsl_vector_get(pointlist,3*i+1);
		points[i].z = gsl_vector_get(pointlist,3*i+2);
	}
	linevec();
	normalvec();
}

inline Arrow shell::dmag_darrow(Arrow a) {					  //d|a|/da
	return a / a.mag();
}

inline Arrow shell::ddot_darrow(Arrow &a, Arrow &b) {         //d(a.b)/da
	return b;
}

inline Tensor shell::dcross_darrow(Arrow &a, Arrow &b) { 	//d(a*b)/da
	Arrow a0(0, b.cor[2], -b.cor[1]);
	Arrow a1(-b.cor[2], 0, b.cor[0]);
	Arrow a2(b.cor[1], -b.cor[0], 0);
	Tensor t(a0,a1,a2);
	return t;
}

inline Tensor shell::dhat_darrow(Arrow &a) {				//d(n.hat())/dn
	double mag_3 = pow(a.mag(),3);
	double x = a.cor[0], y = a.cor[1], z = a.cor[2];
	Arrow a0(y*y + z*z, -x*y, -x*z);
	Arrow a1(-y*x, x*x + z*z, -y*z);
	Arrow a2(-x*z, -y*z, x*x + y*y);
	Tensor t(a0/mag_3, a1/mag_3, a2/mag_3);
	return t;
}

inline Tensor shell::dnorm_dvi(const int &t, const int &vi) {           //dn/dvi
	int v1 = triangles[t].tv1;
	int v2 = triangles[t].tv2;
	int v3 = triangles[t].tv3;
	// cout<<"\tbegin dnorm_dvi\n";
	// cout<<"\tv1="<<v1<<" v2="<<v2<<" v3="<<v3<<endl;
	// cout<<"\tt="<<t<<" vi="<<vi<<endl;
	Arrow a1(points[v1].x, points[v1].y, points[v1].z);
	Arrow a2(points[v2].x, points[v2].y, points[v2].z);
	Arrow a3(points[v3].x, points[v3].y, points[v3].z);
	// cout<<"points[v2"
	Tensor tensor;
	if (vi == v1) {
		tensor = dcross_darrow(a1,a2) - dcross_darrow(a1,a3);
	}
	else if (vi == v2) {
		tensor = dcross_darrow(a2,a3) - dcross_darrow(a2,a1);
	}
	else if (vi == v3) {
		tensor = dcross_darrow(a3,a1) - dcross_darrow(a3,a2) ;
	}
	
	 //tensor.print();
	return tensor;
}

inline Arrow shell::dhs_dvi(const int &l, const int &vi) {
	Arrow dhs_dvi_a;
	double b;
	int vj;
	if (vi == lines[l].lv1) {
		vj = lines[l].lv2;
		Arrow arrowl(points[vj].x - points[vi].x, points[vj].y - points[vi].y, points[vj].z - points[vi].z);
		b = arrowl.mag();
		dhs_dvi_a = arrowl * (ks * (b - b0) / b);
	}
	else if (vi == lines[l].lv2) {
		vj = lines[l].lv1;
		Arrow arrowl(points[vj].x - points[vi].x, points[vj].y - points[vi].y, points[vj].z - points[vi].z);
		b = arrowl.mag();
		dhs_dvi_a = arrowl * (ks * (b - b0) / b);
	}
	return dhs_dvi_a;

}

inline Arrow shell::dhb_dvi(const int &t1, const int &t2, Arrow &n1, Arrow &n2, const int &vi, const int &sign) {
	// int t1,t2, sign;
	Arrow n1_hat, n2_hat, cross;
	
	n1_hat = n1.hat();
	n2_hat = n2.hat();

	// cross = n1_hat * n2_hat;
	// double cdotl = cross.dot(larrow[l]);
	// if (cdotl < 0)
	// 	sign = 1;
	// else
	// 	sign = -1;

	Arrow dcos_dvi, dsin_dvi, dhb_dv_a;
	dcos_dvi = ddot_darrow(n1_hat,n2_hat) * dhat_darrow(n1) * dnorm_dvi(t1,vi) +
			   ddot_darrow(n2_hat,n1_hat) * dhat_darrow(n2) * dnorm_dvi(t2,vi);
	dsin_dvi = dmag_darrow(n1_hat * n2_hat) * (dcross_darrow(n1_hat,n2_hat) * dhat_darrow(n1) * dnorm_dvi(t1,vi) 
			   								  - dcross_darrow(n2_hat,n1_hat) * dhat_darrow(n2) * dnorm_dvi(t2,vi));	
    dhb_dv_a = dcos_dvi * (-kb * cos(th0)) + dsin_dvi * (-kb * sin(th0) * sign);


	// Arrow dcos_dvi1, dsin_dvi1, dcos_dvi2, dsin_dvi2;
 //    dcos_dvi1 = ddot_darrow(n1_hat,n2_hat) * dhat_darrow(n1);
 //    dcos_dvi2 = ddot_darrow(n2_hat,n1_hat) * dhat_darrow(n2);
 //    cout<<"dnormt1_dv"<<vi<<"=";
	// dnorm_dvi(t1,vi).print();
	// cout<<"dnormt2_dv"<<vi<<"=";
	// dnorm_dvi(t2,vi).print();
	// cout<<"dcos_dv"<<vi<<"=";
	// (dcos_dvi * (-kb*cos(th0))).print();
	// cout<<"dsin_dv"<<vi<<"=";
	// (dsin_dvi * (-kb*sin(th0) * sign)).print();
	// cout<<"dhb_dv_a=";
	// dhb_dv_a.print();

    return dhb_dv_a;  
}

double shell::Energy(const gsl_vector *pointlist) {
	//cout<<"start Energy\n";

	double Es = 0, Eb = 0, Et = 0, co, si, cdotl;
	int t1, t2,sign;
	update(pointlist);
	// linevec();
	// normalvec();
	num_l = lines.size();


	for (int i = 0; i < num_l; i++) {
		if (lines[i].ltn == 2) {
			//cout<<"lines[i].ltn == 2 and the length of line:"<<i<<" = "<<larrow[i].mag()<<endl;
			Es += ks * (larrow[i].mag() - b0) * (larrow[i].mag() - b0);

			t1 = lines[i].lt1;
			t2 = lines[i].lt2;
			Arrow n1 = nhatarrow[t1];
			Arrow n2 = nhatarrow[t2];
			co = n1.dot(n2);
			Arrow cross = n1 * n2;
			cdotl = cross.dot(larrow[i]);
			
			if (cdotl > 0)
				sign = 1;
			else
				sign = -1;
			si = cross.mag() * sign;
			Eb += kb * (1 - co * cos(th0) - si * sin(th0));
			//cout<<"co="<<co<<" si="<<cross.mag()<<" sign="<<sign<<" th0="<<th0<<endl;


		}
		else if (lines[i].ltn == 1) {
			//cout<<"lines[i].ltn == 1 and the length of line:"<<i<<" = "<<larrow[i].mag()<<endl;
			Es += 0.5 * ks * (larrow[i].mag() - b0) * (larrow[i].mag() - b0);
		}
	}
	Et = Es + Eb;
	// Et = Es;
	//cout<<"\tEt="<<Et<<endl;
	return Et;
}


void shell::gradE(const gsl_vector *pointlist, gsl_vector *dElist) {
	//cout<<"begin gradE"<<endl;
	// clock_t start, end;
	// start = clock();

	int v0, v1, v2, v3, t1, t2, sign;
	Arrow n1, n2, cross, dhs_dv0, dhs_dv1, dhb_dv0, dhb_dv1, dhb_dv2, dhb_dv3,n1_hat,n2_hat;
	update(pointlist);
	// linevec();
	// normalvec();
	num_v = points.size();
    num_l = lines.size();
    num_t = triangles.size();
    vector<vector<double>> dh_dv(num_v,vector<double>(3,0.));
    for (int l = 0; l < num_l; l++) {
    	v0 = lines[l].lv1;
    	v1 = lines[l].lv2;
    	dhs_dv0 = dhs_dvi(l,v0);
    	dhs_dv1 = dhs_dvi(l,v1);
    	double coef = (lines[l].ltn == 2 ? 2.: 1.);

    	


    	for (int alpha = 0; alpha < 3; alpha++) {
    		dh_dv[v0][alpha] -= coef * dhs_dv0[alpha];
    		dh_dv[v1][alpha] -= coef * dhs_dv1[alpha];
    	}


    	if (lines[l].ltn == 2) {
    		t1 = lines[l].lt1;
    		t2 = lines[l].lt2;
    		n1 = narrow[t1];
    		n2 = narrow[t2];
    		n1_hat = nhatarrow[t1];
    		n2_hat = nhatarrow[t2];


    		cross = n1_hat * n2_hat;
			double cdotl = cross.dot(larrow[l]);
			
			if (cdotl > 0)
				sign = 1;
			else
				sign = -1;

    		if (triangles[t1].tv1 != v0 && triangles[t1].tv1 != v1)
    			v2 = triangles[t1].tv1;
    		else if (triangles[t1].tv2 != v0 && triangles[t1].tv2 != v1)
    			v2 = triangles[t1].tv2;
    		else if (triangles[t1].tv3 != v0 && triangles[t1].tv3 != v1)
    			v2 = triangles[t1].tv3;

    		if (triangles[t2].tv1 != v0 && triangles[t2].tv1 != v1)
    			v3 = triangles[t2].tv1;
    		else if (triangles[t2].tv2 != v0 && triangles[t2].tv2 != v1)
    			v3 = triangles[t2].tv2;
    		
    		
    		dhb_dv0 = dhb_dvi(t1,t2,n1,n2,v0,sign);
    		dhb_dv1 = dhb_dvi(t1,t2,n1,n2,v1,sign);
    		

    		for (int alpha = 0; alpha < 3; alpha++) {
    			dh_dv[v0][alpha] += dhb_dv0[alpha];
    			dh_dv[v1][alpha] += dhb_dv1[alpha];
    		

    			 //cout<<"dhb_dv["<<alpha<<"]="<<dhb_dv0[alpha]<<" "<<dhb_dv1[alpha]<<" "<<dhb_dv2[alpha]<<" "<<dhb_dv3[alpha]<<endl;

    		}
    	}  
    

    for (int i = 0; i < num_v; i++) {
    	for (int alpha = 0; alpha < 3; alpha++) {
    		//cout<<"dh_dv="<<dh_dv[i][alpha]<<endl;
    		gsl_vector_set(dElist,3*i + alpha, dh_dv[i][alpha]);
    	}
    }
   
}

void shell:: EdE(const gsl_vector *pointlist,double *h, gsl_vector *dElist){

    *h=Energy(pointlist);
    int v0, v1, v2, v3, t1, t2, sign;
	Arrow n1, n2, cross, dhs_dv0, dhs_dv1, dhb_dv0, dhb_dv1, dhb_dv2, dhb_dv3,n1_hat,n2_hat;
	update(pointlist);
	// linevec();
	// normalvec();
	num_v = points.size();
    num_l = lines.size();
    num_t = triangles.size();
    vector<vector<double>> dh_dv(num_v,vector<double>(3,0.));
    for (int l = 0; l < num_l; l++) {
    	v0 = lines[l].lv1;
    	v1 = lines[l].lv2;
    	dhs_dv0 = dhs_dvi(l,v0);
    	dhs_dv1 = dhs_dvi(l,v1);
    	double coef = (lines[l].ltn == 2 ? 2.: 1.);


    	for (int alpha = 0; alpha < 3; alpha++) {
    		dh_dv[v0][alpha] -= coef * dhs_dv0[alpha];
    		dh_dv[v1][alpha] -= coef * dhs_dv1[alpha];
    	}


    	if (lines[l].ltn == 2) {
    		t1 = lines[l].lt1;
    		t2 = lines[l].lt2;
    		n1 = narrow[t1];
    		n2 = narrow[t2];
    		n1_hat = nhatarrow[t1];
    		n2_hat = nhatarrow[t2];


    		cross = n1_hat * n2_hat;
			double cdotl = cross.dot(larrow[l]);
			
			if (cdotl > 0)
				sign = 1;
			else
				sign = -1;

    		if (triangles[t1].tv1 != v0 && triangles[t1].tv1 != v1)
    			v2 = triangles[t1].tv1;
    		else if (triangles[t1].tv2 != v0 && triangles[t1].tv2 != v1)
    			v2 = triangles[t1].tv2;
    		else if (triangles[t1].tv3 != v0 && triangles[t1].tv3 != v1)
    			v2 = triangles[t1].tv3;

    		if (triangles[t2].tv1 != v0 && triangles[t2].tv1 != v1)
    			v3 = triangles[t2].tv1;
    		else if (triangles[t2].tv2 != v0 && triangles[t2].tv2 != v1)
    			v3 = triangles[t2].tv2;
    		
    		// cout<<"v0="<<v0<<" v1="<<v1<<" v2="<<v2<<" v3="<<v3<<endl;
    		// cout<<"n1="<<t1<<" n2="<<t2<<endl;
    		
    		dhb_dv0 = dhb_dvi(t1,t2,n1,n2,v0,sign);
    		dhb_dv1 = dhb_dvi(t1,t2,n1,n2,v1,sign);
    		

    		for (int alpha = 0; alpha < 3; alpha++) {
    			dh_dv[v0][alpha] += dhb_dv0[alpha];
    			dh_dv[v1][alpha] += dhb_dv1[alpha];
    			

    			 //cout<<"dhb_dv["<<alpha<<"]="<<dhb_dv0[alpha]<<" "<<dhb_dv1[alpha]<<" "<<dhb_dv2[alpha]<<" "<<dhb_dv3[alpha]<<endl;

    		}
    	}  
    	
    }

    for (int i = 0; i < num_v; i++) {
    	for (int alpha = 0; alpha < 3; alpha++) {
    		//cout<<"dh_dv="<<dh_dv[i][alpha]<<endl;
    		gsl_vector_set(dElist,3*i + alpha, dh_dv[i][alpha]);
    	}
    }
    
    //gradE(pointlist, dElist);
}






double shell::relax() {
	//cout<<"start relax\n";
	int iter = 0;
	int status;
	double min_energy;
	// double stepsize=.0001;
 //    double tol=pow(10,-7);
	double stepsize=.0001;
    double tol=pow(10,-7);

    num_v = points.size();
    num_t = triangles.size();

    // cout<<"num_v="<<num_v<<" points.size()="<<points.size()<<endl;
    // cout<<"num_t="<<num_t<<" triangles.size()="<<triangles.size()<<endl;

    gsl_multimin_function_fdf my_func;
    my_func.n = 3*num_v;
    my_func.f = &shell::Energy_wrapper;
    my_func.df = &shell::gradE_wrapper;
    my_func.fdf = &shell::EdE_wrapper;
    my_func.params = this;

    pointlist = gsl_vector_alloc(3*num_v);
    dElist = gsl_vector_alloc(3*num_v);

    for (int i = 0; i < num_v; i++) {
    	gsl_vector_set(pointlist, 3 * i, points[i].x);
    	gsl_vector_set(pointlist, 3 * i + 1, points[i].y);
    	gsl_vector_set(pointlist, 3 * i + 2, points[i].z);
    }

    double E=Energy(pointlist);
    num_t = triangles.size();
    cout<<"initial E"<<E<<" "<<num_t<<" "<<E/num_t<<endl;

    const gsl_multimin_fdfminimizer_type *type_ptr = gsl_multimin_fdfminimizer_vector_bfgs;
    gsl_multimin_fdfminimizer *minimizer_ptr = gsl_multimin_fdfminimizer_alloc(type_ptr, 3*num_v);
    gsl_multimin_fdfminimizer_set (minimizer_ptr, &my_func, pointlist, stepsize, tol);
    //cout<<"begin iteration"<<endl;
    do {
    	iter++;
    	status = gsl_multimin_fdfminimizer_iterate(minimizer_ptr);
    	//cout<<"update points"<<endl;
    	update(minimizer_ptr->x);
    	if (status) {
            cout<<"Error:"<<status<<endl;
            break;
        }
        //cout<<"check status"<<endl;
        status = gsl_multimin_test_gradient(minimizer_ptr->gradient, tol);
        //cout << "iteration number" << iter << endl;
        //cout<<"E="<<E<<endl;
        if (status == GSL_SUCCESS) {
            // cout << "iteration number" << iter << endl;
            min_energy = gsl_multimin_fdfminimizer_minimum(minimizer_ptr);
            //cout << "energy min" << min_energy << endl;
        }
    } while(status == GSL_CONTINUE && iter<max_iter);
    /*The below part is just added on 11/15/2017*/
    // if I get "error27", it means the gsl has reached a local min and can't progress
    //further. So I will accept that point and below I will get the min_energy
    if (status){
        min_energy=gsl_multimin_fdfminimizer_minimum(minimizer_ptr);
        //cout<<"gsl_stopped, best min occured at"<<iter<<" "<<min_energy<<endl;
    }
     // update(pointlist);
    gsl_multimin_fdfminimizer_free (minimizer_ptr);
    gsl_vector_free (pointlist);
    cout<<"end relax\n";
    return min_energy;
    
}


