#pragma once
//class Mat2;
#include "Mat2.h"

#include <string>

class Vec2{
public:
    complex x, y;

    Vec2(){
        x = 0;
        y = 0;
    };
    Vec2(complex x, complex y){
        this->x = x;
        this->y = y;
    };
    std::string ToString();
    Vec2 operator+(Vec2 b);
    Vec2& operator+=(Vec2 b);
    Vec2 operator-(Vec2 b);
    Vec2& operator-=(Vec2 b);
    
    Vec2 Conj();
    Vec2 operator*(double b);
    complex dot(Vec2 b);

    Vec2& operator/=(double b){
        this->x/=b;
        this->y/=b;
        return *this;
    }

    Vec2& normalize();

    Vec2 dot(Mat2 b);

    complex operator[](int index);
    
};