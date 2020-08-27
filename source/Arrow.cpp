

#include <stdio.h>
#include <iostream>
#include <cmath>

#include "Arrow.h"
#include "Tensor.h"


double& Arrow::operator[](const int& index){
    return cor[index];
}
void Arrow::operator=(const double& number){
    for (unsigned int i=0; i<3; ++i) {
        cor[i]=number;
    }
}
Arrow Arrow::operator+(Arrow a){
    Arrow tmp(cor[0]+a.cor[0],cor[1]+a.cor[1],cor[2]+a.cor[2]);
    return tmp;
}
Arrow Arrow::operator-(Arrow a){
    Arrow tmp(cor[0]-a.cor[0],cor[1]-a.cor[1],cor[2]-a.cor[2]);
    return tmp;
}
double Arrow::dot(Arrow a){
    return cor[0]*a.cor[0]+cor[1]*a.cor[1]+cor[2]*a.cor[2];
}
Arrow Arrow::operator*(Arrow a){
    Arrow tmp(cor[1]*a.cor[2]-cor[2]*a.cor[1],
              cor[2]*a.cor[0]-cor[0]*a.cor[2],
              cor[0]*a.cor[1]-cor[1]*a.cor[0]);
    return tmp;
}
Arrow Arrow::operator*(const double& numble){
    Arrow tmp(numble*cor[0],numble*cor[1],numble*cor[2]);
    return tmp;
}

Arrow Arrow::operator*(Tensor t){
	Arrow tmp(dot(t.column(0)),dot(t.column(1)),dot(t.column(2)));
    return tmp;
}
Arrow Arrow::operator/(const double& numble){
    Arrow tmp(cor[0]/numble,cor[1]/numble,cor[2]/numble);
    return tmp;
}
double Arrow::mag(){
    return sqrt(pow(cor[0],2)+pow(cor[1],2)+pow(cor[2],2));
}
Arrow Arrow::hat(){
    double m=mag();
    Arrow tmp(cor[0]/m,cor[1]/m,cor[2]/m);
    return tmp;
}
void Arrow::print(){
    std::cout<<"{"<<cor[0]<<","<<cor[1]<<","<<cor[2]<<"}"<<std::endl;
}
