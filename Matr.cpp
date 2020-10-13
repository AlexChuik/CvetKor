#include "Matr.hpp"

double formula (int i, int j, int n) //с ее помощью инициализирую матрицу 
{
    return 0;
}

Pixel::Pixel() 
{
    for(int i=0; i<3;i++) {
        RGB[i] = 0;
        linRGB[i] = 0;
        lab[i] = 0;
    }
}
Pixel::~Pixel(){}
Pixel::Pixel(sf::Color p)
{
    RGB[0] = (double)p.r/255;
    if(RGB[0] <= 0.04045) linRGB[0] = RGB[0]/(12.92);
    else linRGB[0] = pow((RGB[0]+0.055)/(1.055),2.4);
    
    RGB[1] = (double)p.g/255;
    if(RGB[1] <= 0.04045) linRGB[1] = RGB[1]/(12.92);
    else linRGB[1] = pow((RGB[1]+0.055)/(1.055),2.4);
    
    RGB[2] = (double)p.b/255;
    if(RGB[2] <= 0.04045) linRGB[2] = RGB[2]/(12.92);
    else linRGB[2] = pow((RGB[2]+0.055)/(1.055),2.4);

    lab[0] = (linRGB[0] - linRGB[1])/sqrt(2);
    lab[1] = (2*linRGB[2] - linRGB[0] - linRGB[1])/sqrt(6);
    lab[2] = (linRGB[2] + linRGB[0] + linRGB[1])/sqrt(3);
}


Matr::Matr() {
    m=0;
    n=0;
    a = new double[n*m];
}

Matr::~Matr() {
    delete[] a;
}

Matr::Matr(int r, int d) {
    n=r;
    m=d;
    a = new double[n*m];
    int i,j;
    for (i=0;i<n;i++)
    {
        for (j=0; j<m;j++) 
            this->index(i,j)=formula(i,j,n);
    }
}

Matr::Matr(int r) {
    n=r;
    m=r;
    a=new double[n*m];
    int i,j;
    for (i=0;i<n;i++)
    {
        for (j=0; j<m;j++) 
            this->index(i,j)=formula(i,j,n);
    }
}

size_t Matr::getStl () const {
    return m;
}
size_t Matr::getStr () const {
    return n;
}

Matr Matr::operator+ (/*Matr const & c,*/const Matr & b) {
    int i, j;
    if(( (*this).getStl() != b.getStl() )||( (*this).getStr() != b.getStr() )) {
        printf("Error Size1\n");
        return *this;
    }
    
    Matr* z = new Matr((*this).getStr(), (*this).getStl());
    for ( i = 0; i < (*this).getStr(); i++) {
        for ( j = 0; j < (*this).getStl(); j++) 
            (*z).index(i,j) = (*this).index(i,j) + b.index(i,j);
    
    }
    return (*z);
}

double & Matr::index (size_t i, size_t j) {
    return a[m*i+j];
}
double Matr::index(size_t i, size_t j) const {
    return a[m*i+j];
}

Matr & Matr::operator= (const Matr & c) {
    int i, j;
    if (( c.getStl() != this->getStl() )||( c.getStr() != this->getStr() )) {
        printf("Error Size2\n");
        return *this;
    }
    for ( i = 0; i < c.getStr(); i++) {
        for ( j = 0; j < c.getStl(); j++) 
            (*this).index(i,j) = c.index(i,j);
    
    }

    return *this;
}

Matr Matr::operator* ( const Matr & b) {
    int i, j, l;
    Matr* z = new Matr((*this).getStr(), b.getStl());

    for (i = 0; i < (*this).getStr(); i++) {
        for ( j = 0; j < b.getStl(); j++) {
            (*z).index(i,j) = 0;
           
            for ( l = 0; l < (*this).getStl(); l++) {
                 (*z).index(i,j) += (*this).index(i,l)*b.index(l,j);
            }
        }
    }
    return (*z);
}

Matr Matr::operator* (double b) {
    int i, j;
    
    Matr* z = new Matr((*this).getStr(), (*this).getStl());

    for (i = 0; i < (*this).getStr(); i++) {
        for ( j = 0; j < (*this).getStl(); j++) 
            (*z).index(i,j) = (*this).index(i,j)*b;
    }
    return (*z);
}


void Matr::Print()
{
    int i,j;
    if (m<21&&n<11) { 
    for(i=0;i<n;i++)
    {
       for (j=0;j<m;j++)
       {
          printf ("%10.6lf  ", (*this).index(i,j));

       }
          printf("\n");
    }
    printf ("\n");
    }
}

void Matr::readMatr(const char * path) {
    FILE * f;
    int i, j;
    double tmp;

    f = fopen (path, "r");
    if (f == NULL) {
        printf ("Can't open file\n");
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
               if(fscanf(f,"%lf",&tmp) != 1)
               {
                   printf("%d %d\n",i,j);
                   printf(" Error read\n");
                   return;
               }
               (*this).index(i,j) = tmp;
          }
     }

  fclose (f);
}

void Matr::obrez(int k) {
    double* b = new double[n*m];
    int i,j;
    
    for (i=0;i<k;i++) {
        for (j=0; j<k;j++) 
            b[i*k+j] = this->index(i,j);
    }
    n=k;
    m=k;
    a = b;
}

double Matr::Norma() {
    int i,c;
    double norm;
    norm = 0;
   // if((n != 1)&&(m != 1)) return -1;
    for (i=0;i<n;i++) {
        for (c=0; c<m; c++) 
            norm += ((*this).index(i,c))*((*this).index(i,c));
    } 
   
    return sqrt(norm);
}

void Matr::obrxod1 ()  //ne rabochii
{
    int i,j,c;
    if (n > m) return;
    for( i = 0; i < n; i++)  
    {
        for(c=0; c<m; c++) {
            if (abs((*this).index(i,i))>EPS)
                (*this).index(i,c) = (*this).index(i,c)/((*this).index(i,i));
            else printf("Shok tyt problemu v grebannom gausse\n");
        }

        for( j = i+1; j<n; j++) {
            for(c = 0; c < m; c++) 
                (*this).index(j,c) = ((*this).index(j,c))-((*this).index(i,c))*((*this).index(j,i));    
        }
    }
}
void Matr::obrxod2 ()
{
    int i,j,c;
    double prick;
    if (n > m) return;
    for( i = n-1; i >= 0; i--)  
    {
        prick = (*this).index(i,i);
        for(c=0; c<m; c++) {
            if (abs((*this).index(i,i))>EPS) 
                (*this).index(i,c) = (*this).index(i,c)/(prick);
            //printf("%lf   %lf\n",(*this).index(i,i),(*this).index(i,c));}
            else printf("Opredel. = 0\n");
        }
        
        for( j = i-1; j>=0; j--) {
            prick = (*this).index(j,i);
            for(c=0; c<m; c++) 
               // printf("%lf  %lf  %lf\n",(*this).index(j,c),(*this).index(i,c),prick);
                (*this).index(j,c) = ((*this).index(j,c))-((*this).index(i,c))*(prick);  
               // printf("%lf   %lf\n",(*this).index(j,c),(*this).index(j,i));}  
        }
    }
}

Matr Matr::trans() {
    Matr T(m,n);
    double tmp;
    int i,j;
    tmp=0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) 
               T.index(j,i)=(*this).index(i,j);
    }
    return T;
}
