

#ifndef Arrow_hpp
#define Arrow_hpp

#include <stdio.h>

class Tensor;
//ARROW CLASS
class Arrow{
public:
    //coordinate:
    double cor[3];
    //initial:
    Arrow():cor{0.,0.,0.}{}
    Arrow(const double x,const double y,const double z):cor{x,y,z}{}
    Arrow(const double *p):cor{p[0],p[1],p[2]}{}
    Arrow(const double *p1,const double *p2):cor{p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2]}{}
    //functions:
    double& operator[](const int&);
    void  operator=(const double&);
    Arrow operator+(Arrow);
    Arrow operator-(Arrow);
    double dot(Arrow);
    Arrow operator*(Arrow);
    Arrow operator*(const double&);
    Arrow operator*(Tensor);
    Arrow operator/(const double&);
    double mag();
    Arrow hat();
    
    void print();
};


#endif /* Arrow_hpp */