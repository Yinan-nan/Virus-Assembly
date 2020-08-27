

#ifndef Tensor_hpp
#define Tensor_hpp

#include <stdio.h>
#include "Arrow.h"

class Tensor{
public:
    //coordinates:
    Arrow row[3];
    //functions:
    Tensor();
    Tensor(Arrow a0, Arrow a1, Arrow a2):row{a0,a1,a2}{}
    
    Arrow column(const int&);
    Tensor transpose();
    Arrow& operator[](const int&);
    Arrow  operator*(Arrow);
    Tensor operator-(Tensor);
    Tensor operator*(Tensor);

    void print();
};

#endif /* Tensor_hpp */
