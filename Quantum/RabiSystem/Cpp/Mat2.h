#pragma once
#include <complex>
#define complex std::complex<double>
const complex I = complex(0., 1.);
class Vec2;


class Mat2{

private:
    complex m11, m12, m21, m22;
public:
    Mat2(){
        m11=m12=m21=m22= I * 0.0;
    }
    Mat2(complex m11, complex m12, complex m21, complex m22){
        this->m11 = m11;
        this->m12 = m12;
        this->m21 = m21;
        this->m22 = m22;
    }



    Mat2 Hermitian();
    Mat2 dot(Mat2 &b);
    Mat2 operator*(complex c);

    Mat2& operator*=(complex c);

    Mat2 operator*(double c);
    Mat2& operator*=(double c);

    Mat2 operator+(Mat2 c);

    Mat2& operator+=(Mat2 c);
    Mat2 operator-(Mat2 c);

    Vec2 dot(Vec2 b);

};
inline Mat2 DOT(Mat2 &a, Mat2 &b){
    return a.dot(b);
};
