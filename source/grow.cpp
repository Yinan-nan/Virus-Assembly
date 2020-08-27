

#include <iostream>
//#include "grow.h"
#include "shell.h"
#include <fstream>
#include <math.h>
#include <sys/stat.h>
/*#include "energy.cpp"*/
using namespace std;

void shell::grow() {
    clock_t start, end;
    start = clock();
	int con;
	int step = 1;
	while (num_t < MAX_num_tri) {
        cout<<"num of triangles="<<num_t<<endl;
		theta.clear();
        con = nextmove();      
        outputfile(step);
        step++;
        relax();
        outputfile(step);
        step++;
        cout<<"step="<<step<<endl;
        cout<<"con="<<con<<endl;
        if (con == 1000)
        	break;
        
	}
	cout<<"finish"<<endl;
    end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC); 
     
    cout << "\tTime taken by program is : " << time_taken <<" sec\n";
	//outputfile();
}

int shell::nextmove() {
	int thi = openangle();
	int v;
	double Epentamer, Ehexamer;
    if (thi == 1000)
    	return 1000;
    else {
    	begin:
    	if (theta[thi].type == 0) {
    		v = theta[thi].li;
    		switch (points[v].nn) {

    			case 5:
    			cout << "Case 5:";
    			if (theta[thi].th <= 0.5) { 
    				
    				merge_lines(v);
    			}
    			else if (theta[thi].th >= 0.6) {
    			
    				
    				insert_tri(v);
    		    }
    		    else if (theta[thi].th < 0.6 && theta[thi].th > 0.5) {
    		    	
                    //cout<<"before backup"<<endl;
    		    	backup();
                    //cout<<"before merge line"<<endl;
    		    	merge_lines(v);
    		    	Epentamer = relax() / triangles.size();
                   
                    cout<<Epentamer<<endl;
                    
    		    	restore();
    		    	insert_tri(v);
    		    	Ehexamer = relax() / triangles.size();
                    
                    cout<<Ehexamer<<endl;
    		    	if (Ehexamer < Epentamer) {
    		    		
    		    	}
    		    	else {
    		    		restore();
    		    		merge_lines(v);
    		    		
    		    	}
    		    }
    		    break;

    		    case 6:
    		    
    		    merge_lines(v);
    		    break;
    		}
    	}
    	else if (theta[thi].type == 1) {
            cout << "attach a new triangle in front of the line" << endl;
            nextvertex(theta[thi].li);
        }
        theta.clear();
        return 1;
    }

}

void shell::nextvertex(int ln) {
	num_v = points.size();
	num_l = lines.size();
	num_t = triangles.size();
    int tn = lines[ln].lt1;
    int lline = lines[ln].ll;
    int rline = lines[ln].lr;
    int lpt = lines[ln].lv1;;
    int rpt = lines[ln].lv2;

    Arrow arrnh = nhatarrow[tn];
    Arrow arrl = larrow[ln].hat();;
    Arrow drt = arrl * arrnh;
    double x, y, z,l;
    l = sqrt(1. - 1. / 4. * larrow[ln].mag() * larrow[ln].mag());
    x = (points[lpt].x + points[rpt].x) / 2. + l * drt.cor[0] * cos(th0) - l * arrnh.cor[0] * sin(th0);
    y = (points[lpt].y + points[rpt].y) / 2. + l * drt.cor[1] * cos(th0) - l * arrnh.cor[1] * sin(th0);
    z = (points[lpt].z + points[rpt].z) / 2. + l * drt.cor[2] * cos(th0) - l * arrnh.cor[2] * sin(th0);
    
    vertices newv(x,y,z,1,1);
    points.push_back(newv);

    linevertices la(lpt, num_v, num_t, 1000, 1, lline, num_l+1);
    linevertices lb(num_v, rpt, num_t, 1000, 1, num_l, rline);
    lines.push_back(la);
    lines.push_back(lb);

    trianglevertices newtri(lpt,num_v,rpt,num_l,num_l + 1,ln,0);
    triangles.push_back(newtri);
    
    points[lpt].neighbor_l.push_back(num_l);
    points[rpt].neighbor_l.push_back(num_l + 1);
    points[lpt].neighbor_t.push_back(num_t);
    points[rpt].neighbor_t.push_back(num_t);
    points[lpt].nn += 1;
    points[rpt].nn += 1;

  
    lines[ln].lt2 = num_t;
    lines[ln].ltn += 1;
    lines[lline].lr = num_l;
    lines[rline].ll = num_l + 1;

    num_v = points.size();
	num_l = lines.size();
	num_t = triangles.size();
    cout<<"\tpoint: ";
    for (int i=0; i<points.size(); i++) {
        cout << points[i].x << " " << points[i].y << " " << points[i].z <<"\n";
    }
        
}
void shell::insert_tri(int vertex) {
	int vertex_lline, vertex_rline,lline,rline;
	int lv, rv;
	num_l = lines.size();
	num_t = triangles.size();
	for (int i = 0;i < points[vertex].neighbor_l.size();i++) {
		int b = points[vertex].neighbor_l[i];

		if (lines[b].ltn == 1) {
			if (lines[b].lv2 == vertex) {
				vertex_lline = b;
			}
			else if (lines[b].lv1 == vertex) {
				vertex_rline = b;
			}
		}
	}
    lv = lines[vertex_lline].lv1;
    rv = lines[vertex_rline].lv2;
    lline = lines[vertex_lline].ll; //
    rline = lines[vertex_rline].lr;
    
    linevertices newl (lv,rv,num_t,1000,1,lline,rline);
    lines.push_back(newl);

    lines[lline].lr = num_l; //I don't know the meaning of this procedure. ll and lr are only valid when they are on the edge;
    lines[rline].ll = num_l;

  

    points[vertex].nn += 1;
    points[lv].nn += 1;
    points[rv].nn += 1;
    points[vertex].edge = 0;
    points[vertex].neighbor_t.push_back(num_t);
  

    trianglevertices newtri(vertex,lv,rv,vertex_lline,num_l,vertex_rline,0);
    triangles.push_back(newtri);

    num_v = points.size();
	num_l = lines.size();
	num_t = triangles.size();
    // cout<<"\tpoint: ";
    // for (int i=0; i<points.size(); i++) {
    //     cout << points[i].x << " " << points[i].y << " " << points[i].z <<"\n";
    // }
}

void shell::merge_lines(int vertex) {
	num_v = points.size();
	num_l = lines.size();
	num_t = triangles.size();
	points[vertex].edge = 0;
    cout<<"vertex="<<vertex<<endl;
	int vertex_lline,vertex_rline;
    cout<<"\t neighbor_l num:"<<points[vertex].neighbor_l.size()<<endl;
    int label;
	for (int i = 0;i < points[vertex].neighbor_l.size();i++) {
		int b = points[vertex].neighbor_l[i];
        // cout<<"b="<<b<<endl;
        // cout<<"adjacent triangle num"<<lines[b].ltn<<endl;
		if (lines[b].ltn == 1) {
            // cout<<"b.lv2="<<lines[b].lv2<<endl;
            // cout<<"b.lv1="<<lines[b].lv1<<endl;
			if (lines[b].lv2 == vertex) {
				vertex_lline = b;
			}
			else if (lines[b].lv1 == vertex) {
				vertex_rline = b;
                label = i;
			}
		}
	}
    points[vertex].neighbor_l.erase(points[vertex].neighbor_l.begin() + label);

    // cout<<"\t beginning\n";
	int ltri = lines[vertex_lline].lt1;
	int rtri = lines[vertex_rline].lt1;
    // cout<<"\t first\n";
	int a = lines[vertex_lline].lv1;
	int b = lines[vertex_rline].lv2;
     // cout<<"\t second\n";
	int lline = lines[vertex_lline].ll;
	int rline = lines[vertex_rline].lr;

  
	if (lines[l1].lv1 == vertex && lines[l1].lv2 == b) {  
		triangles[rtri].tl1 = vertex_lline;
	}
	else if (lines[l2].lv1 == vertex && lines[l2].lv2 == b) {
		triangles[rtri].tl2 = vertex_lline;
	}
	else if (lines[l3].lv1 == vertex && lines[l3].lv2 == b){
		triangles[rtri].tl2 = vertex_lline;
	}
	lines[vertex_lline].lt2 = rtri;
	lines[vertex_lline].ltn = 2;
    lines[rline].ll = lline;
    lines[lline].lr = rline;
    
    // cout<<"\t second\n";
    //changing the b->a for all the triangles having b
    for (int i = 0;i < num_t;i++) {
    	if (triangles[i].tv1 == b) {  
    		triangles[i].tv1 = a;
    	}
    	else if (triangles[i].tv2 == b) {
    		triangles[i].tv2 = a;
    	}
    	else if (triangles[i].tv3 == b) {
    		triangles[i].tv3 = a;
    	}
	}
    // cout<<"\t third\n";

      
    points[a].nn += points[b].nn;
    for (int i = 0;i < points[b].neighbor_l.size();i++) {
    	int nbl = points[b].neighbor_l[i];
    	if (lines[nbl].lv1 != vertex && lines[nbl].lv2 != vertex)
    		points[a].neighbor_l.push_back(nbl);
    }
    for (int i = 0;i < points[b].neighbor_t.size();i++) {
    	points[a].neighbor_t.push_back(points[b].neighbor_t[i]);
    }
    

    points.erase(points.begin() + b);
    lines.erase(lines.begin() + vertex_rline);
    
    // change the index of lines and poitns because of erasing 
    num_v = points.size();
	num_l = lines.size();
	num_t = triangles.size();
	for (int i = 0; i < num_v; i++) {
		for (int j = 0; j < points[i].neighbor_l.size(); j++) {
			if (points[i].neighbor_l[j] > vertex_rline) 
				points[i].neighbor_l[j] -= 1;
		}
	}
	for (int i = 0; i < num_l; i++) {
		if (lines[i].ll > vertex_rline)
			lines[i].ll -= 1;
		if (lines[i].lr > vertex_rline)
			lines[i].lr -= 1;
	
	}
	
	num_v = points.size();
	num_l = lines.size();
	num_t = triangles.size();
    // cout<<"\tpoint: ";
    // for (int i=0; i<points.size(); i++) {
    //     cout << points[i].x << " " << points[i].y << " " << points[i].z <<"\n";
    // }
}



void shell::restore() {
	points.clear();
	lines.clear();
	triangles.clear();
	points = pointsbackup;
	lines = linesbackup;
	triangles = trianglesbackup;
}

int shell::check_closeshell (int v_on_edge, int index) {
    int l1 = -1,l2 = -1,l3 = -1,v1 = -1,v2 = -1,v3 = -1;
    num_t = triangles.size();
    if (v_on_edge == 3 && theta[index].th < 1) {
        cout<<"close the shell with a triangle"<<endl;
        if (theta[index].type == 1) {
            l1 = theta[index].li;
            l2 = lines[l1].lr;
            l3 = lines[l1].ll;
            v1 = lines[l1].lv1;
            v2 = lines[l1].lv2;
            v3 = lines[l2].lv2;
            if (lines[l2].lv2 != lines[l3].lv1) {
                cout<<"there is a mistake in closing the shell-code_1"<<endl;
                return 0;
            }
        }
        if (theta[index].type == 0) {
            v1 = theta[index].li;
            for(int i = 0; i < points[v1].neighbor_l.size();i++) {
                int l = points[v1].neighbor_l[i];
                if(lines[l].ltn == 1) {
                    if (lines[l].lv1 == v1)
                        l1 = l;
                    else if(lines[l].lv2 == v1)
                        l3 = l;
              }
            }
            if (l1 != -1 && l3 != -1){
                l2 = lines[l1].lr;
                v2 = lines[l1].lv2;
                v3 = lines[l3].lv1;
            }
        }
        
       
        lines[l2].ltn = 2;
        lines[l2].lt2 = num_t;
        lines[l3].ltn = 2;
        lines[l3].lt2 = num_t;
       

        trianglevertices tri(v1,v3,v2,l1,l3,l2,0);
        triangles.push_back(tri);
        return 1000;
    }
    else if (v_on_edge==2) {
        cout<<"close the shell with merging the last two lines"<<endl;
        if (theta[index].type == 1) {
            cout<<"there is a problem in closing-code 2"<<endl;
            return 0;
        }
        else {
            v1 = theta[index].li;
            for (int i = 0;i < points[v1].neighbor_l.size();i++) {
                int l = points[v1].neighbor_l[i];
                if (lines[l].ltn == 1) {
                    if (lines[l].lv1 == v1)
                        l1 = l;
                    else if (lines[l].lv2 == v1)
                        l2 = l;
                }
            }
            if (lines[l1].lv2 != lines[l2].lv1)
                cout<<"there is a problem in closing-code 2"<<endl;
            v2 = lines[l1].lv2;
            int l, lmax;
            l = min(l1,l2);
            lmax = max(l1,l2);
            lines[l].ltn = 2;
            int t1 = lines[lmax].lt1;
            lines[l].lt2 = t1;
           
            for (int i = 0;i < points[v1].neighbor_l.size();i++) {
                if (points[v1].neighbor_l[i] == lmax)
                    points[v1].neighbor_l.erase(points[v1].neighbor_l.begin() + i);
            }
            for (int i = 0;i < points[v2].neighbor_l.size();i++) {
                if (points[v2].neighbor_l[i] == lmax)
                    points[v2].neighbor_l.erase(points[v2].neighbor_l.begin() + i);
            }
            lines.erase(lines.begin() + lmax);
            if (lmax != lines.size()) {
                cout<<"the erased line in closing is not the last line and we need reindexing"<<endl;
                for (int i = 0; i < num_t;i++) {
                    if (triangles[i].tl1 > lmax) 
                        triangles[i].tl1 -= 1;
                    if (triangles[i].tl2 > lmax)
                        triangles[i].tl2 -= 1;
                                
                }
                num_l = lines.size(); 
                for (int i = 0; i < num_l;i++) {
                    if (lines[i].ll > lmax)
                        lines[i].ll -= 1;
                    if (lines[i].lr > lmax)
                        lines[i].lr -= 1;
                }
            }

            return 1000;
        }
    }
    else  {      
        return 0;
    }
}

void shell::outputfile(int step) {

	
	ofstream file;
	file.open("output.txt");
	file << "yinan\n";

	file << "points' size: "<<points.size()<<"\n";
	for (int i=0; i<points.size(); i++) {
		file <<"vertex coordinates: \n";
        file << points[i].x << " " << points[i].y << " " << points[i].z << " " << points[i].nn <<"\n";
        file <<"neighbor_l: \n";
        for (int j = 0;j < points[i].neighbor_l.size();j++) {
        	file << points[i].neighbor_l[j] <<" ";
        }
        file << "\n" << "neighbor_t: \n";
        for (int j = 0;j < points[i].neighbor_t.size();j++) {
        	file << points[i].neighbor_t[j] <<" ";
        }
        file << "\n";

    }
   
    file << "triangles' size: " << triangles.size() << "\n";
    for (int i = 0; i < triangles.size(); ++i)
    {
        file <<triangles[i].tv1<<" "<<triangles[i].tv2<<" "<<triangles[i].tv3<<" "<<
                   triangles[i].tl1<<" "<<triangles[i].tl2<<" "<<triangles[i].tl3<<"\n";
    }
    file.close();


	string File = "Movie";
	string Path = File + "/";
	remove(File.c_str());
	mkdir(File.c_str(),S_IRWXU);
	string Name = "step=" + to_string(step);
	
    ofstream file2;
	file2.open(Path+Name+"-Parameter.txt");
	// file2 << "yinan\n";
	file2 << "points' size: "<<points.size()<<"\n";
	for (int i=0; i<points.size(); i++) {
        file2 << points[i].x << " " << points[i].y << " " << points[i].z <<"\n";
    }

    for (int i = 0; i < points.size(); i++) {
         file2 << points[i].neighbor_l.size()<<"\n";
    }

   
    
    file2 << "triangles' size: " << triangles.size() << "\n";
    for (int i = 0; i < triangles.size(); ++i)
    {
        file2 <<triangles[i].tv1<<" "<<triangles[i].tv2<<" "<<triangles[i].tv3<<"\n";
    }
    file2.close();
}













