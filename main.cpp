#include "Matr.hpp"

using namespace sf; //Область имён для работы сущностей SFML
using namespace std;
static char* Imia;

double Max(double a, double b) 
{
    if (a > b) return a;
    return b;
}

void display()
{
    //рисую систему координат и плоскости
    glClear(GL_COLOR_BUFFER_BIT);
    glPushMatrix();
    
    glTranslatef(-350, -550, -900);
    glRotated(-90, 1, 0, 0);
    glRotated(-40, 0, 0, 1);
    glColor3f(1.0, 1.0, 1.0);
    for (float x = -128; x <= 128; x += 16) {
        glBegin(GL_LINE_STRIP);
        glColor4f(0.3, 0.3, 0.7,0.1);
        for (float y = -128; y <= 128; y += 16)
        {
            //if ((-x-y+255/2>=0)&&(-x-y+255/2<=1))
            glVertex3f(x+255/2, y+255/2, -x-y+255/2);
        }
        glEnd();

        glBegin(GL_LINE_STRIP);
        glColor4f(0.3, 0.3, 0.7,0.1);
        for (float y = -128; y <= 128; y += 16)
        {
            //if ((-x-y+255/2>=0)&&(-x-y+255/2<=1))
            glVertex3f(y+255/2, x+255/2, -x-y+255/2);
        }
        glEnd();
    }
    for (float x = 0; x <= 256; x += 16) {
        glPointSize(2);
        glBegin(GL_LINE_STRIP);
        glColor3f(0.3, 0.7, 0.7);
        for (float y = 0; y <= 256; y += 16)
        {
            glVertex3f(y, y, y);
        }
        glEnd();
        glPointSize(1);
        glBegin(GL_LINE_STRIP);
        glColor3f(0.3, 0.3, 0.3);
        for (float y = 0; y <= 256; y += 16)
        {
            glVertex3f(x, y, 0);
        }
        glEnd();
        glBegin(GL_LINE_STRIP);
        glColor3f(0.3, 0.3, 0.3);
        for (float y = 0; y <= 256; y += 16)
        {
            glVertex3f(0, y, x);
        }
        glEnd();
        glBegin(GL_LINE_STRIP);
        glColor3f(0.3, 0.3, 0.3);
        for (float y = 0; y <= 256; y += 16)
        {
            glVertex3f(0, x, y);
        }
        glEnd();
        glBegin(GL_LINE_STRIP);
        glColor3f(0.3, 0.3, 0.3);
        for (float y = 0; y <= 256; y += 16)
        {
            glVertex3f(y, x, 0);
        }
        glEnd();
        glBegin(GL_LINE_STRIP);
        glColor3f(0.3, 0.3, 0.3);
        for (float y = 0; y <= 256; y += 16)
        {
            glVertex3f(y, 256, x);
        }
        glEnd();
        glBegin(GL_LINE_STRIP);
        glColor3f(0.3, 0.3, 0.3);
        for (float y = 0; y <= 256; y += 16)
        {
            glVertex3f(x, 256, y);
        }
        glEnd();
    }
    //закончил рисовать
    

    Image heroimage; //Создание единичного экземпляра стороннего класса Image
    heroimage.loadFromFile(Imia); //Чтение изображения в экземпляр класса
    //Папку Image располагать рядом с Sourse.cpp текущего проекта
    int x = heroimage.getSize().x; //Выводим в консоль размер изображения ширина
    int y = heroimage.getSize().y; //Выводим в консоль размер изображения высота
    //PCA analis
    Matr m(1,3); //центр масс в линРГБ
    Matr C(3),E(3),A(3);
    E.index(0,0)=1; E.index(1,1)=1; E.index(2,2)=1;
    vector<Pixel> claster(x * y); //кластер цветов пикселей 
    Matr** linRGBclaster = new Matr*[x * y];
    Matr** labclaster = new Matr*[x * y];
    double MSE(0),MSE1(0);
    for (unsigned l = 0; l < x; ++l) //Листаем пиксели в копии изображения в экземпляре heroimage
    { 
        for (unsigned p = 0; p < y; ++p)
        { 
            linRGBclaster[l*y+p] = new Matr(1,3);
            labclaster[l*y+p] = new Matr(1,3);
            claster[l*y+p] = Pixel(heroimage.getPixel(l, p)); //Считывание цвета для пиксела с заданными координатами
            for (int i = 0; i < 3; i++) 
            {
                linRGBclaster[l*y+p]->index(0,i) = claster[l*y+p].linRGB[i];
                m.index(0,i)+=claster[l*y+p].linRGB[i];
                labclaster[l*y+p]->index(0,i) = claster[l*y+p].lab[i];
            }
            MSE += pow(labclaster[l*y+p]->index(0,0),2) + pow(labclaster[l*y+p]->index(0,1),2);
            //сразу их рисую (красным)
            *linRGBclaster[l*y+p] = (*linRGBclaster[l*y+p]) * 255;
            glBegin(GL_POINTS);
            glPointSize(1);
            glColor3f(1.0, 0, 0);
            glVertex3f(linRGBclaster[l*y+p]->index(0,0), 
            linRGBclaster[l*y+p]->index(0,1), linRGBclaster[l*y+p]->index(0,2));
            glEnd();
            *linRGBclaster[l*y+p] = (*linRGBclaster[l*y+p]) * (1./255);
            //закончил рисовать
        } 
    } 
    MSE = (sqrt(MSE))/(x*y);
    MSE1=MSE;
    cout << "Изначальное MSE: " << MSE << endl;
    MSE = 0;
    m.index(0,0)/=x*y;
    m.index(0,1)/=x*y;
    m.index(0,2)/=x*y;// m = E(c)  
    cout << "m: ";
    m.Print();

    for (unsigned l = 0; l < x; ++l) //Делаем матрицу ковариаций
    { 
        for (unsigned p = 0; p < y; ++p)
        { 
            C = C + ((*(linRGBclaster[l*y+p])) + m*(-1)).trans()*((*(linRGBclaster[l*y+p])) + m*(-1));
        } 
    }
    C = (C*((double)1/(x*y))); 
    cout << "Ковариационная матрица: " << endl;
    C.Print();
    
    A = C;
    Matr u(1,3), uX1(1,3), uX2(1,3); //с окончанием Х тот же смысл но для анализа FHT
    Matr lam(1,3);

    Sobstv3nach(A,lam);
    A=C;
    lam.index(0,1) = Max(lam.index(0,1), lam.index(0,2));
    cout << lam.index(0,1) <<endl;
    //---------------------------------------------------------------------------------
    SobstvVector(A, u); //Вычисляю главную ось
    MainAxisXaf(x, y, labclaster, uX1, uX2, lam.index(0,1));//главная ось по Хафу
    //---------------------------------------------------------------------------------

    cout << "u: ";
    u.Print(); 
    
    Matr g(1,3), c(1,3), d(1,3), S(3);
    for(int i=0;i<3;i++) {  //задаю параметры для адаптации 
        for(int j=0;j<3;j++) {
            S.index(i,j)=0;
            if(i==j) S.index(i,j) = (u.index(0,0)+u.index(0,1)+u.index(0,2))/(3*u.index(0,i));
        }
    }
    d.index(0,0) = 1./sqrt(3); d.index(0,1) = 1./sqrt(3); d.index(0,2) = 1./sqrt(3); //ось яркости
    c.index(0,0) = 1./2.; c.index(0,1) = 1./2.; c.index(0,2) = 1./2.; //центр куба
    //вычисление серой точки
    g = m + u * ((d * (c + (m*(-1)) ).trans())*(1/((d*u.trans()).index(0,0)))).index(0,0);
    cout << "g: ";
    g.Print();
    //рисую найденные параметры кластера
    glPointSize(5);
    glBegin(GL_POINTS);
    glColor3f(1.0, 1.0, 0);
    glVertex3f(uX1.index(0,0)*255, uX1.index(0,1)*255, uX1.index(0,2)*255);
    glColor3f(1.0, 0, 1.0);
    glVertex3f(uX2.index(0,0)*255, uX2.index(0,1)*255, uX2.index(0,2)*255);
    glEnd();

    glPointSize(5);
    glBegin(GL_POINTS);
    glColor3f(0, 1.0, 0);
    glVertex3f(g.index(0,0)*255, g.index(0,1)*255, g.index(0,2)*255);
    glColor3f(0, 0, 1.0);
    glVertex3f(m.index(0,0)*255, m.index(0,1)*255, m.index(0,2)*255);
    glEnd();
    glPointSize(5);
    glBegin(GL_LINE_STRIP);
    glColor3f(1.0, 0, 1.0);
    glVertex3f(g.index(0,0)*255, g.index(0,1)*255, g.index(0,2)*255);
    glVertex3f(m.index(0,0)*255, m.index(0,1)*255, m.index(0,2)*255);
    glEnd();
    glPointSize(1);


    for (unsigned l = 0; l < x; ++l) //Делаю цветокорекцию      
    { 
        for (unsigned p = 0; p < y; ++p)
        { 
            (*linRGBclaster[l*y+p]) = (*linRGBclaster[l*y+p]) + g*(-1);
            (*linRGBclaster[l*y+p]) =  (*linRGBclaster[l*y+p])*S;
            (*linRGBclaster[l*y+p]) = (*linRGBclaster[l*y+p]) + c;
        } 
    }
    
    cout << "Мaтрица S" << endl;
    S.Print();
        
    for (unsigned l = 0; l < x; ++l) //Меняю пиксели на фото и рисую новый кластер
    { 
        for (unsigned p = 0; p < y; ++p)
        { 
            *linRGBclaster[l*y+p] = (*linRGBclaster[l*y+p]) * 255;
            glBegin(GL_POINTS);
            glPointSize(1);
            glColor3f(0, 0, 1.0);
            glVertex3f(linRGBclaster[l*y+p]->index(0,0), 
            linRGBclaster[l*y+p]->index(0,1), linRGBclaster[l*y+p]->index(0,2));
            glEnd();
            *linRGBclaster[l*y+p] = (*linRGBclaster[l*y+p]) * (1./255);
            
            //преобразование назад в ргб
            for (int i = 0; i < 3; i++)
            {
                if(linRGBclaster[l*y+p]->index(0,i) <= 0.0031308) 
                linRGBclaster[l*y+p]->index(0,i) = linRGBclaster[l*y+p]->index(0,i) * (12.92);
                else 
                linRGBclaster[l*y+p]->index(0,i) = pow((linRGBclaster[l*y+p]->index(0,i)),1./2.4) * (1.055) - 0.055; 
                if(linRGBclaster[l*y+p]->index(0,i) > 1 ) 
                *linRGBclaster[l*y+p] = (*linRGBclaster[l*y+p])*(1./(linRGBclaster[l*y+p]->index(0,i)));
                if(linRGBclaster[l*y+p]->index(0,i) < 0)
                linRGBclaster[l*y+p]->index(0,i) = 0;
            }   
            *linRGBclaster[l*y+p] = (*linRGBclaster[l*y+p]) * 255;

            Color ct(linRGBclaster[l*y+p]->index(0,0),
            linRGBclaster[l*y+p]->index(0,1),
            linRGBclaster[l*y+p]->index(0,2)
            );
            for (int i = 0; i < 3; i++) 
            {
                labclaster[l*y+p]->index(0,i) = Pixel(ct).lab[i];
            }
            MSE += pow(labclaster[l*y+p]->index(0,0),2) + pow(labclaster[l*y+p]->index(0,1),2);
            heroimage.setPixel(l, p, ct);
        } 
    } 
    MSE = (sqrt(MSE))/(x*y);
    cout << "Конечное MSE: " << MSE << endl;
    cout << "Начальное на конечное: " << MSE1/MSE << endl;
    MSE = 0;
    heroimage.saveToFile("Преобразованное.jpg"); //Сохранение преобразованного изображения

    //теперь делаю все то же но для параметров найденных с помощью преобразования Хафа
    for (unsigned l = 0; l < x; ++l)
    { 
        for (unsigned p = 0; p < y; ++p)
        {
            for (int i = 0; i < 3; i++) 
            {
                linRGBclaster[l*y+p]->index(0,i) = claster[l*y+p].linRGB[i];
            }
        } 
    } 

    u = uX2 + uX1*(-1);
    m = uX1;
    cout << "uX: ";
    u.Print();
    g = m + u * ((d * (c + (m*(-1)) ).trans())*(1/((d*u.trans()).index(0,0)))).index(0,0);
    cout << "gX: ";
    g.Print();
    for(int i=0;i<3;i++) {  //задаю параметры для адаптации 
        for(int j=0;j<3;j++) {
            S.index(i,j)=0;
            if(i==j) S.index(i,j) = (u.index(0,0)+u.index(0,1)+u.index(0,2))/(3*u.index(0,i));
        }
    }
    
    for (unsigned l = 0; l < x; ++l) //Делаю цветокорекцию      
    { 
        for (unsigned p = 0; p < y; ++p)
        { 
            (*linRGBclaster[l*y+p]) = (*linRGBclaster[l*y+p]) + g*(-1);
            (*linRGBclaster[l*y+p]) =  (*linRGBclaster[l*y+p])*S;
            (*linRGBclaster[l*y+p]) = (*linRGBclaster[l*y+p]) + c;
        } 
    }
    for (unsigned l = 0; l < x; ++l) //Меняю пиксели на фото и рисую новый кластер
    { 
        for (unsigned p = 0; p < y; ++p)
        { 
            *linRGBclaster[l*y+p] = (*linRGBclaster[l*y+p]) * 255;
            glBegin(GL_POINTS);
            glPointSize(1);
            glColor3f(1.0, 1.0, 1.0);
            glVertex3f(linRGBclaster[l*y+p]->index(0,0), 
            linRGBclaster[l*y+p]->index(0,1), linRGBclaster[l*y+p]->index(0,2));
            glEnd();
            *linRGBclaster[l*y+p] = (*linRGBclaster[l*y+p]) * (1./255);
            
            //преобразование назад в ргб
            for (int i = 0; i < 3; i++)
            {
                if(linRGBclaster[l*y+p]->index(0,i) <= 0.0031308) 
                linRGBclaster[l*y+p]->index(0,i) = linRGBclaster[l*y+p]->index(0,i) * (12.92);
                else 
                linRGBclaster[l*y+p]->index(0,i) = pow((linRGBclaster[l*y+p]->index(0,i)),1./2.4) * (1.055) - 0.055; 
                if(linRGBclaster[l*y+p]->index(0,i) > 1 ) 
                *linRGBclaster[l*y+p] = (*linRGBclaster[l*y+p])*(1./(linRGBclaster[l*y+p]->index(0,i)));
                if(linRGBclaster[l*y+p]->index(0,i) < 0)
                linRGBclaster[l*y+p]->index(0,i) = 0;
            }   
            *linRGBclaster[l*y+p] = (*linRGBclaster[l*y+p]) * 255;

            Color ct(linRGBclaster[l*y+p]->index(0,0),
            linRGBclaster[l*y+p]->index(0,1),
            linRGBclaster[l*y+p]->index(0,2)
            );
            for (int i = 0; i < 3; i++) 
            {
                labclaster[l*y+p]->index(0,i) = Pixel(ct).lab[i];
            }
            MSE += pow(labclaster[l*y+p]->index(0,0),2) + pow(labclaster[l*y+p]->index(0,1),2);
            delete linRGBclaster[l*y+p];
            delete labclaster[l*y+p];
            heroimage.setPixel(l, p, ct);
        } 
    } 
    MSE = (sqrt(MSE))/(x*y);
    cout << "Конечное MSE (Хаф): " << MSE << endl;
    cout << "Начальное на конечное (Хаф): " << MSE1/MSE << endl;
    glPopMatrix();
    glutSwapBuffers();
    heroimage.saveToFile("ПреобразованноеХаф.jpg"); //Сохранение преобразованного изображения

    cout << "Готово" << endl;
    delete[] linRGBclaster;
    delete[] labclaster;
    
}


int main(int argc, char* argv[])
{ 
    setlocale(LC_ALL, "rus");
    if(argc != 2) 
    {
        cout << "Укажите название изображения " << endl; //нет проверки корректности названия
        return 0;
    }
    //все ниже для openGL
    Imia = argv[1];
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(800, 800);
    glutInitWindowPosition(20, 86);
    glutCreateWindow("Preobrazov");
    glClearColor(0, 0, 0, 1.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-400, 200, -500, 100, 500, 10000);
    glMatrixMode(GL_MODELVIEW);
    glutDisplayFunc(display);
    glutMainLoop();
    return 0;
}
    
