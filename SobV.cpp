#include "Matr.hpp"

void QRraz(Matr &a, Matr &Q, Matr &R);

Matr vuchMatrVrash1(Matr &a, int k,int n) {  
    double sk, cosinus, sinus, ModVec;
    int j,i,c;
    k=k+2;
    Matr* T1 = new Matr(n);
    Matr* Tpom = new Matr(n);
    for(i=0;i<n;i++) { //делаем Е
        for(c=0;c<n;c++) {
            T1->index(i,c)=0;
            if(i==c) T1->index(i,c) = 1;
        }
    }
    for(j=k;j<n;j++) {
        ModVec = sqrt(a.index(k-1,k-2)*a.index(k-1,k-2)+a.index(j,k-2)*a.index(j,k-2));
        if(ModVec<EPS) continue;
        cosinus = a.index(k-1,k-2)/ModVec; 
        sinus = (-1)*(a.index(j,k-2)/ModVec);
        for(i=0;i<n;i++) {
             for(c=0;c<n;c++) {
                 T1->index(i,c)=0;
                 if(i==c) T1->index(i,c) = 1;
             }
        }
        
        T1->index(k-1,k-1) = cosinus;
        T1->index(j,j) = cosinus;
        T1->index(k-1,j) = (-1)*sinus;
        T1->index(j,k-1) = sinus;
        *Tpom = a;
        for(i=0;i<n;i++) { //умножаем Т1*a
            Tpom->index(k-1,i) = (a.index(k-1,i))*cosinus + (-1)*(a.index(j,i))*sinus;
            Tpom->index(j,i) =  (a.index(k-1,i))*sinus + (a.index(j,i))*cosinus;
        }
        a = *Tpom;

        *Tpom = a;
        for(i=0;i<n;i++) { //умножаем a*Т1.trans()
            Tpom->index(i,k-1) = (a.index(i,k-1))*cosinus + (-1)*(a.index(i,j))*sinus;
            Tpom->index(i,j) =  (a.index(i,k-1))*sinus + (a.index(i,j))*cosinus;
        }
        a = *Tpom;
    } 
    return (*T1);
}


double CrutimVertim(Matr &a, int n) {
    //разлагаем матрицу А-sI потом собираем и приб sI
    int i,j;                               
    Matr* S = new Matr(n);
    Matr* Pom = new Matr(n);
    Matr* Q = new Matr(n);
    Matr* R = new Matr(n);
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            S->index(i,j)=0;
            if(i==j) S->index(i,j) = 1;
        }
    }
    *S = (*S)*a.index(n-1,n-1);
    *Pom = (a+(*S)*(-1));
    QRraz(*Pom,*Q,*R);
    a = (*R)*(*Q) + (*S);
    return(a.index(n-1,n-1));
}


void QRraz(Matr &a, Matr &Q, Matr &R) {  //только для почти треуг
    double cosinus, sinus, ModVec;
    int i,k,n;
    n = a.getStr();
    Matr* Tpom = new Matr(n);

    for(i=0;i<n;i++) {
        for(k=0;k<n;k++) {
            Q.index(i,k) = 0;
            if(i==k)
                Q.index(i,k) = 1;
        }
    }   
    R = a;
    for(k=0;k<n-1;k++)
    {
        ModVec = sqrt(R.index(k,k)*R.index(k,k)+R.index(k+1,k)*R.index(k+1,k));
        if(ModVec<EPS) continue;
        cosinus = R.index(k,k)/ModVec; 
        sinus = (-1)*(R.index(k+1,k)/ModVec);

        *Tpom = R;
        for(i=0;i<n;i++) {
            Tpom->index(k,i) = (R.index(k,i))*cosinus + (-1)*(R.index(k+1,i))*sinus;
            Tpom->index(k+1,i) =  (R.index(k,i))*sinus + (R.index(k+1,i))*cosinus;
        }
        R = *Tpom;
        *Tpom = Q;
        for(i=0;i<n;i++) { 
            Tpom->index(k,i) = (Q.index(k,i))*cosinus + (-1)*(Q.index(k+1,i))*sinus;
            Tpom->index(k+1,i) =  (Q.index(k,i))*sinus + (Q.index(k+1,i))*cosinus;
        }
        Q = *Tpom;
    } 
    Q = Q.trans();
}
void SobstvVector (Matr &a, Matr &u) 
{
    int n(3);//, arr[3]={0,1,2};//перестановочный массив
    Matr T(n),Qinf(n),Q(n),R(n);
    double temp(0);
    for(int i=0;i<n;i++) { //делаем Е
        for(int j=0;j<n;j++) {
            Qinf.index(i,j)=0;
            if(i==j) Qinf.index(i,j) = 1;
        }
    }
    T = vuchMatrVrash1(a,0,n);
    for(int i = 0; i < 150; i++) 
    {
        QRraz(a,Q,R);
        Qinf = Qinf * Q;
        a = (R)*(Q);
    }
    Qinf = (T.trans())*Qinf;
    u.index(0,0) = Qinf.index(0,0);
    u.index(0,1) = Qinf.index(1,0);
    u.index(0,2) = Qinf.index(2,0); 
}


void Sobstv3nach (Matr &a, Matr &lambda) 
{
    int n(3);
    int i,k,j;
    //приводим матрицу к почти треугольному виду-----------------------------------------------------------
    vuchMatrVrash1(a,0,n);
    //----------------------------------------------------------------------------------------------------------------
    //Строим QR разложение почти треуг матрицы со сдвигами перемножаем RQ и так далее --------------------------------
    int shet(0);
    for(int shet = 0; shet < 50; shet++) {
        CrutimVertim(a,3);
    }
    lambda.index(0,0) = a.index(0,0);
    lambda.index(0,1) = a.index(1,1);
    lambda.index(0,2) = a.index(2,2);
}
//Собственные вектора и значения ------------------------------------------------------------------
//Собственные вектора и значения ------------------------------------------------------------------
//Собственные вектора и значения ------------------------------------------------------------------
