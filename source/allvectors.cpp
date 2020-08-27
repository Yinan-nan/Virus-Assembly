// #include "allvectors.h"
#include <iostream>
#include<vector>
#include"shell.h"
#include <math.h>
/*#include "parameters.h"*/
#include <fstream>
using namespace std;

void shell::linevec() {
    double vx, vy, vz, mag;
    int   left, right, num_l;

    lvector.clear();
    larrow.clear();
    rvlarrow.clear();
    num_l  =lines.size();
    for (int i = 0; i < num_l; ++i) {
        left = lines[i].lv1;
        right = lines[i].lv2;
        vx = points[right].x - points[left].x;
        vy = points[right].y - points[left].y;
        vz = points[right].z - points[left].z;
        mag = sqrt((vx * vx) + (vy * vy) + (vz * vz));

        vectors vec(vx,vy,vz,mag);
        lvector.push_back(vec);

        Arrow arr(vx, vy, vz);
        larrow.push_back(arr);

        Arrow rvarr(-vx,-vy,-vz);
        rvlarrow.push_back(rvarr);
    }
}

void shell::normalvec() {    //point to outside
    double nvec_x,nvec_y,nvec_z,nvec_mag;
    double nx,ny,nz;

    nhat.clear();
    nhatarrow.clear();
    narrow.clear();
    num_t = triangles.size();


        //arrow mehtod
        Arrow arrv1(points[v1].x, points[v1].y, points[v1].z);
        Arrow arrv2(points[v2].x, points[v2].y, points[v2].z);
        Arrow arrv3(points[v3].x, points[v3].y, points[v3].z);
        Arrow arrnh = arrv1 * arrv2 + arrv2 * arrv3 + arrv3 * arrv1;
        narrow.push_back(arrnh);
        arrnh = arrnh / arrnh.mag();
        nhatarrow.push_back(arrnh);
    
    // cout<<"normal vector"<<endl;
    // cout<<"vectors\n";
    // for (int i = 0;i < nhat.size();i++) {
    //     cout<<nhat[i].x<<","<<nhat[i].y<<","<<nhat[i].z<<endl;
    // }
    // cout<<"Arrow"<<"\n";
    // for (int i = 0;i < nhatarrow.size();i++) {
    //     cout<<nhatarrow[i].cor[0]<<","<<nhatarrow[i].cor[1]<<","<<nhatarrow[i].cor[2]<<endl;
    // }
}

int shell::openangle() {
    double theta1,theta2, angle;
    int lv, rv, j, k,lline,rline;
    num_l=lines.size();
    num_t=triangles.size();
    num_v=points.size();

    linevec();
    normalvec();
    //vertex-angle
    for (int m=0; m<num_v; m++){
        if(points[m].nn > 3 and points[m].edge == 1){  // if points[].nn==4, calculated both vertexangle and lineangle
            for (int i = 0;i < points[m].neighbor_l.size();i++) {
                int b = points[m].neighbor_l[i];
                // cout<<"lin:"<<b<<" ";
                // cout<<"ll:"<<lines[b].lv1<<" "<<"lr:"<<lines[b].lv2<<endl;
                if (lines[b].ltn == 1) {
                    if (lines[b].lv2 == m)  {
                        lline = b;
                        // cout<<"lr:"<<lines[b].lr<<endl;
                    }
                    if (lines[b].lv1 == m) {
                        rline = b;
                        // cout<<"ll:"<<lines[b].ll<<endl;
                    }
                }
            }
            // cout<<"lline:"<<lline<<" "<<"rline:"<<rline<<endl;
            angle=lineangle(lline,rline,0,m);
            openangles nangle(m, angle,0);
            theta.push_back(nangle);
        }
    }
    //line-openangle
    for (int i = 0; i < num_l; ++i) { 
        if (lines[i].ltn ==1) {
            // cout<<"i:"<<i<<endl;
            // cout<<"v1:"<<lines[i].lv1<<" "<<"v2:"<<lines[i].lv2<<endl;
            

            if(points[lv].nn<=4 and points[rv].nn<=4) {
                // cout<<"j:"<<j<<" "<<"i:"<<i<<endl;
                theta1 = lineangle(j, i, 0, lv); 
                theta2 = lineangle(i, k, 0, rv); // sign is 0, means the direction of cross and triangle is the same，1 means opposite direction.
             
                angle = min(theta1, theta2);
                openangles nangle(i, angle, 1);
                theta.push_back(nangle);
            }

        } else {
            // cout<<i<<" the line is not on the edge\n";
        }
    }
    int t = theta.size();
    for(int p=0; p<t; p++){
        //cout<<"openangle:"<<theta[p].li<<" "<<theta[p].th<<" "<<theta[p].type<<endl;
    }
    int angle_index = chooseangle();
    cout<<"angle was chosed:"<<angle_index<<endl;
    return angle_index;
}



int shell::chooseangle() {
    int v_on_edge=0,closeshell,index = -1;
    double mintheta = 2 * 3.14;
    for (int i = 0; i < theta.size(); ++i) {
        // mintheta = min(mintheta, theta[i].th);
        if (theta[i].th < mintheta) {
            index = i;
            mintheta = theta[i].th;
        }
    }
    for(int k=0; k<points.size(); k++) {
        if (points[k].edge==1)
            v_on_edge++;
    }
    if(v_on_edge < 4) {
        closeshell = check_closeshell(v_on_edge,index);
        if (closeshell == 1000)
            return 1000;
        else
            return index;
    }
    else
        return index;
}














