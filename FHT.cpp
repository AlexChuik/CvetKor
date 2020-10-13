#include "Matr.hpp"
using namespace sf;

//параметр дискретизации
const int N = 512; //в пространстве sRGB координаты дискретны от 0 до 225, поэтому этого должно хватить 

//самая длинная ось l примерно 1.73, координаты получены преобразованием дискретных
double f (int k) {
    if (k == 0) return sqrt(2);
    if (k == 1) return 4./sqrt(6);
    return 3./sqrt(3);
}

void razmutGauss(int x, int y, int koor, double sigma, Matr** labclaster, double* razmclaster)
{
    int n = (int)(6*sigma+1);
    double* window = new double[n];
    double* temp = new double[N*N];

    for(int i = 0; i < n; i++) { //создаю окно размытия
        window[i] = pow(M_E, ((i - (n-1)/2)*(i - (n-1)/2))/(-2*sigma*sigma));
    }
    for(int i = 0; i < N*N; i++) { 
        temp[i] = 0;
        razmclaster[i] = 0;
    }
    for(int i = 0; i < x; i++) {  //размываем по горизонтали
        for(int j = 0; j < y; j++) {
            for(int k = 0; k < n; k++) {
                if (((int)(((labclaster[i*y+j]->index(0,koor))+f(koor)/2)*N/f(koor)) + (k - (n-1)/2)>=0)
                && ((int)(((labclaster[i*y+j]->index(0,koor))+f(koor)/2)*N/f(koor)) + (k - (n-1)/2) < N))
                {
                    temp[(int)(((labclaster[i*y+j]->index(0,koor))+f(koor)/2)*N/f(koor)) *N  //дискр коор по ширине
                    + (k - (n-1)/2)*N//разброс влево вправо
                    + (int)((labclaster[i*y+j]->index(0,2))*N/f(2))] //дискр коор по высоте
                    += window[k];
                }
            } 
        }
    }
    for(int i = 0; i < N; i++) {  //размываем по вертикали
        for(int j = 0; j < N; j++) {
            for(int k = 0; k < n; k++) {
                if((j + (k - (n-1)/2) >= 0) && (j + (k - (n-1)/2) < N)){
                    razmclaster[i*N + j + (k - (n-1)/2)] += 
                    temp[i*N + j]*window[k];
                }    
            } 
        }
    }
    delete[] temp;
    delete[] window;
}

void Xaf1(int x, int y, int koor, double* mass, int pos = 1) //вертикальные с наклоном вправо
{
    if(pos == N) return;
    double* mass2 = new double[N*N];
    for(int k = 0; k < N; k += 2*pos) 
    {
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < pos; j++)
            {   
                mass2[i*N+2*j+k] = mass[i*N+j+k] + mass[((i+j+2*N)%N)*N+j+k+pos];
                mass2[i*N+2*j+1+k] = mass[i*N+j+k] + mass[((i+1+j+2*N)%N)*N+j+k+pos];
            }
        }
    }

    for(int k = 0; k < N*N; k++)
    {
        mass[k] = mass2[k];
    }
    pos = pos*2;
    delete[] mass2;
    Xaf1(x, y, koor, mass, pos);
}

void Xaf2(int x, int y, int koor, double* mass, int pos = 1) //вертикальные с наклоном влево
{
    if(pos == N) return;
    double* mass2 = new double[N*N];
    for(int k = 0; k < N; k += 2*pos)
    {
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < pos; j++)
            {   
                mass2[i*N+2*j+k] = mass[i*N+j+k] + mass[((i-j+2*N)%N)*N+j+k+pos];
                mass2[i*N+2*j+1+k] = mass[i*N+j+k] + mass[((i-1-j+2*N)%N)*N+j+k+pos];
            }
        }
    }
    for(int k = 0; k < N*N; k++) 
    {
        mass[k] = mass2[k];
    }
    pos = pos*2;
    delete[] mass2;
    Xaf2(x, y, koor, mass, pos);
}


void MainAxisXaf(int x, int y, Matr **labclaster, Matr &uX1, Matr &uX2, double l) 
{
    Image heroimage; //буду рисовать размытые проекции и результаты преобразований Хафа для контроля
    heroimage.create(N,N,Color(0,0,0));
    Matr M1(1,2), M2(1,2);
    double max1(0), max2(0), sigma(l*sqrt(2)*N);
    double* Xafmas = new double[N*N];
    double* Xafmas2 = new double[N*N];
    razmutGauss(x,y,0,sigma,labclaster,Xafmas);
    for(int k = 0; k < N*N; k++) 
    {
        Xafmas2[k] = Xafmas[k];
    }
    for(int i = 0; i < N; i++) { //рисую размытие la плоскости
        for(int j = 0; j < N; j++) {
            Color ct(Xafmas[i*N+j],
            Xafmas[i*N+j],
            Xafmas[i*N+j]
            );
            heroimage.setPixel(i, j, ct);
        }
    }
    heroimage.saveToFile("TestRazmLA.jpg");


    Xaf1(x, y, 0, Xafmas2); //перобразование Хафа для вертик. вправо и рисунок
    for(int i = 0; i < N; i++) { 
        for(int j = 0; j < N; j++) {
            Color ct(Xafmas2[i*N+j]/1000,
            Xafmas2[i*N+j]/1000,
            Xafmas2[i*N+j]/1000
            );
            heroimage.setPixel(i, j, ct);
        }
    }
    heroimage.saveToFile("TestXaf1LA.jpg");

    for(int i = 0; i < N; i++) { //Поиск наилучшей прямой
        for(int j = 0; j < N; j++) {
            if (max1 < Xafmas2[i*N+j]) {
                M1.index(0,0) = i;
                M1.index(0,1) = j;
                max1 = Xafmas2[i*N+j];
            }
        }
    }

    for(int k = 0; k < N*N; k++) 
    {
        Xafmas2[k] = Xafmas[k];
    }
    
    Xaf2(x, y, 0, Xafmas2); //перобразование Хафа для вертик. влево и рисунок
    for(int i = 0; i < N; i++) { 
        for(int j = 0; j < N; j++) {
            Color ct(Xafmas2[i*N+j]/1000,
            Xafmas2[i*N+j]/1000,
            Xafmas2[i*N+j]/1000
            );
            heroimage.setPixel(i, j, ct);
        }
    }
    heroimage.saveToFile("TestXaf2LA.jpg");

    for(int i = 0; i < N; i++) {   //Поиск наилучшей прямой
        for(int j = 0; j < N; j++) {
            if (max2 < Xafmas2[i*N+j]) {
                M2.index(0,0) = i;
                M2.index(0,1) = j;
                max2 = Xafmas2[i*N+j];
            }
        }
    }
    if (max2 > max1) {M1 = M2; M1.index(0,1) = M1.index(0,1)*(-1);} //выбираю между наклоном вправо и влево
    //M1.Print();


    //аналогичные действия для плоскости lb
    Matr u1(2), u2(2);
    u1.index(1,0) = (M1.index(0,0) + M1.index(0,1)/2) * f(0)/N - f(0)/2;
    u1.index(0,0) = (M1.index(0,0)) * f(0)/N - f(0)/2;
    u1.index(0,1) = 0;
    u1.index(1,1) = (N/2) * f(2)/N;
    max1 = 0; max2 = 0;
    
    
    razmutGauss(x,y,1,sigma,labclaster,Xafmas);
    for(int k = 0; k < N*N; k++) 
    {
        Xafmas2[k] = Xafmas[k];
    }
    for(int i = 0; i < N; i++) {  //рисую размытие lb плоскости
        for(int j = 0; j < N; j++) {
            Color ct(Xafmas[i*N+j],
            Xafmas[i*N+j],
            Xafmas[i*N+j]
            );
            heroimage.setPixel(i, j, ct);
        }
    }
    heroimage.saveToFile("TestRazmLB.jpg");
    
    Xaf1(x, y, 1, Xafmas2);
    for(int i = 0; i < N; i++) { 
        for(int j = 0; j < N; j++) {
            Color ct(Xafmas2[i*N+j]/1000,
            Xafmas2[i*N+j]/1000,
            Xafmas2[i*N+j]/1000
            );
            heroimage.setPixel(i, j, ct);
        }
    }
    heroimage.saveToFile("TestXaf1LB.jpg");

    for(int i = 0; i < N; i++) { 
        for(int j = 0; j < N; j++) {
            if (max1 < Xafmas2[i*N+j]) {
                M1.index(0,0) = i;
                M1.index(0,1) = j;
                max1 = Xafmas2[i*N+j];
            }
        }
    }

    for(int k = 0; k < N*N; k++) 
    {
        Xafmas2[k] = Xafmas[k];
    }

    Xaf2(x, y, 1, Xafmas2);
    for(int i = 0; i < N; i++) { 
        for(int j = 0; j < N; j++) {
            if (max2 < Xafmas2[i*N+j]) {
                M2.index(0,0) = i;
                M2.index(0,1) = j;
                max2 = Xafmas2[i*N+j];
            }
        }
    }

    for(int i = 0; i < N; i++) { 
        for(int j = 0; j < N; j++) {
            Color ct(Xafmas2[i*N+j]/1000,
            Xafmas2[i*N+j]/1000,
            Xafmas2[i*N+j]/1000
            );
            heroimage.setPixel(i, j, ct);
        }
    }
    heroimage.saveToFile("TestXaf2LB.jpg");

    if (max2 > max1) {M1 = M2; M1.index(0,1) = M1.index(0,1)*(-1);}
    //M1.Print();
    delete[] Xafmas;
    delete[] Xafmas2;
    u2.index(1,0) = (M1.index(0,0) + M1.index(0,1)/2) * f(1)/N - f(1)/2;
    u2.index(0,0) = (M1.index(0,0)) * f(1)/N - f(1)/2;
    u2.index(0,1) = 0;
    u2.index(1,1) = (N/2) * f(2)/N;

    //склейка 2х проекций, потом трансформация в linRGB
    uX1.index(0,0) = (sqrt(3)*u1.index(0,0) - u2.index(0,0))/sqrt(6);
    uX1.index(0,1) = ((-1)*sqrt(3)*u1.index(0,0) - u2.index(0,0))/sqrt(6);
    uX1.index(0,2) = 2*(u2.index(0,0))/sqrt(6);

    uX2.index(0,0) = (sqrt(2)*u2.index(1,1) + sqrt(3)*u1.index(1,0) - u2.index(1,0))/sqrt(6);
    uX2.index(0,1) = (sqrt(2)*u2.index(1,1) + (-1)*sqrt(3)*u1.index(1,0) - u2.index(1,0))/sqrt(6);
    uX2.index(0,2) = (sqrt(2)*u2.index(1,1) + 2*u2.index(1,0))/sqrt(6);
}
