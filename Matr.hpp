#ifndef _MATR_H_
#define _MATR_H_

#define EPS 1e-12
#define eps 1e-12

#include <SFML/Graphics.hpp> //Для работы с изображением
#include <GL/glut.h>
#include <iostream>
#include <iomanip>
#include <math.h>

class Matr {
    private:
        double *a;
        size_t n; //кол строк
        size_t m; //столбцов
    public:
        Matr();
        ~Matr();
        Matr(int r);
        Matr(int r,int d);
        Matr operator+ (Matr const & b);
        double & index(size_t i, size_t j);
        double index(size_t i, size_t j) const;
        Matr & operator= (Matr const & b);
        Matr operator* (Matr const & b);
        Matr operator* (double b);
        size_t getStr() const;
        size_t getStl() const;
        void readMatr(const char * path);
        void Print();
        double Norma();
        void obrxod1();
        void obrxod2();
        void obrez(int k);
        Matr trans();
};

void Sobstv3nach (Matr &a, Matr &lambda);
void SobstvVector (Matr &a, Matr &u); 
void MainAxisXaf(int x, int y, Matr **labclaster, Matr &uX1, Matr &uX2, double l);

class Pixel
{
    public:
    double RGB[3];
    double linRGB[3];
    double lab[3];
    Pixel();
    ~Pixel();
    Pixel(sf::Color p);
};

#endif
