
#include <stdio.h>
#include <iostream>

#include "Tensor.h"
#include "Arrow.h"

Tensor::Tensor(){
	row[0]=0;row[1]=0;row[2]=0;
}
Arrow  Tensor::column(const int& index){
    Arrow col(row[0][index],row[1][index],row[2][index]);
    return col;
}
Tensor Tensor::transpose(){
	Tensor tmp_t(column(0),column(1),column(2));
	return tmp_t;
}
Arrow& Tensor::operator[](const int& index){
    return row[index];
}
Arrow Tensor::operator*(Arrow a){
    Arrow tmp_a(row[0].dot(a),row[1].dot(a),row[2].dot(a));
    return tmp_a;
}
Tensor Tensor::operator-(Tensor t){
    Tensor tmp_t(row[0]-t[0],row[1]-t[1],row[2]-t[2]);
    return tmp_t;
}
Tensor Tensor::operator*(Tensor t){
	
    Tensor tmp_t(row[0]*t,row[1]*t,row[2]*t);
    return tmp_t;
}
void Tensor::print(){
	std::cout<<"{";
	row[0].print();
	std::cout<<",";
	row[1].print();
	std::cout<<",";
	row[2].print();
	std::cout<<"}"<<std::endl;
}