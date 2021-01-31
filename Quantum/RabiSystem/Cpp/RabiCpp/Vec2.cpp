#pragma once

#include "Vec2.h"
#include <iostream>
#include <complex>
Vec2 Vec2::operator+(Vec2 b){

    return Vec2(this->x + b.x, this->y + b.y);

}
Vec2& Vec2::operator+=(Vec2 b){
    this->x += b.x;
    this->y += b.y;
    return *this;
}
Vec2 Vec2::operator-(Vec2 b){
    return Vec2(this->x - b.x, this->y - b.y);
}
Vec2& Vec2::operator-=(Vec2 b){
    this->x -= b.x;
    this->y -= b.y;
    return *this;
}

Vec2& Vec2::normalize(){

    double mag = real(conj(this->x)*this->x + conj(this->y)*this->y);
    mag = sqrt(mag);
    this-> x /= mag;
    this-> y /= mag;
    return *this;
}

complex Vec2::dot(Vec2 b){

    return (this->x * b.x) + (this->y * b.y);

}
Vec2 Vec2::operator*(double b){
    return Vec2(this->x*b, this->y * b);
}

Vec2 Vec2::dot(Mat2 b){

    return Vec2(this->x * b.M11() + this->y * b.M21(), this->x * b.M12() + this->y * b.M22());

}
complex Vec2::operator[](int index){
    if(index == 0){
        return this->x;
    }
    if(index == 1){
        return this->y;
    }
    std::cout << "Can only access elements 0 or 1" << std::endl;
    return 0;
}

Vec2 Vec2::Conj(){
    return Vec2(conj(this->x), conj(this->y));
}

std::string Vec2::ToString() {


    std::string x = std::to_string(real(this->x)) +" + " + std::to_string(imag(this->x)) + "i";
    std::string y = std::to_string(real(this->y)) +" + " + std::to_string(imag(this->y)) + "i";
    return x + ", " + y;


}