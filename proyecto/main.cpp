#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdlib>
#define MAX 10
using namespace std;
int Menu(int op1)
{
         cout<<endl;
         cout<<"\t\t\t ÉÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ»\n";
         cout<<"\t\t\t º        Aroximacion polinomial         º\n";
         cout<<"\t\t\t ÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼\n";
         cout<<"\t\t\t º     Resolver por:                     º\n";
         cout<<"\t\t\t º     1. Factorizacion  LU              º\n";
         cout<<"\t\t\t º     2. Eliminacion Gaussiana          º\n";
         cout<<"\t\t\t º     3. Salir.                         º\n";
         cout<<"\t\t\t ÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼\n";

         cout<<endl<<"\t Ingrese Opcion  ";
         cout<<endl<<"\t ÄÄÄÄÄÄÄÄÄÄÄÄÄÄ : ";
         cin>>op1;
         return op1;


}
int SubMenu(int op2)
{
         cout<<endl;
         cout<<"\t\t\t ÉÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ»\n";
         cout<<"\t\t\t º        Factorizacion  LU              º\n";
         cout<<"\t\t\t ÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼\n";
         cout<<"\t\t\t º     Resolver por:                     º\n";
         cout<<"\t\t\t º     1. Doolittle                      º\n";
         cout<<"\t\t\t º     2. Crout                          º\n";
         cout<<"\t\t\t º     3. Cholesky                       º\n";
         cout<<"\t\t\t º     4. Salir.                         º\n";
         cout<<"\t\t\t ÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼\n";

         cout<<endl<<"\t Ingrese Opcion  ";
         cout<<endl<<"\t ÄÄÄÄÄÄÄÄÄÄÄÄÄÄ : ";
         cin>>op2;
         return op2;

}
void titulo()
{
     cout<<"\t\t\t ÉÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ»\n";
     cout<<"\t\t\t º                  AJUSTE DE UNA CURVA POLINOMIAL                º\n";
     cout<<"\t\t\t º               POR EL CRITERIO DE MINIMOS CUADRADOS             º\n";
     cout<<"\t\t\t ÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼\n\n";
}
void lectura_datos(double *x,double *y,int N)
{

   cout<<"\n Ingrese los valores del eje-x:\n";                //Entrada de los valores del eje-x
    for (int i=0;i<N;i++)
        cin>>x[i];
    cout<<"\n Ingrese los valores del eje-y:\n";                //Entrada de los valores del eje-y
    for (int i=0;i<N;i++)
        cin>>y[i];
}
void SumatoriasXin(double *X,double *x,int N,int n)
{
    for (int i=0;i<2*n+1;i++)
    {
        X[i]=0;
        for (int j=0;j<N;j++)
            X[i]=X[i]+pow(x[j],i);
    }
}
void SumatoriasYiXin(double *Y,double *x,double *y,int N,int n)
{
    for (int i=0;i<n+1;i++)
    {
        Y[i]=0;
        for (int j=0;j<N;j++)
        Y[i]=Y[i]+pow(x[j],i)*y[j];
    }
}
double **arreglo_matricial_SXin(double **A,double *X,int n )
{
    A=new double*[n+1];
    for (int i=0;i<=n;i++)
    {
        A[i]=new double[n+2];
        for (int j=0;j<=n;j++)
        {
            A[i][j]=X[i+j];
        }

    }
    return A;
}
void arreglo_matricialAumentado_SYiXin(double **A,double *Y,int n )
{
    for (int i=0;i<=n;i++)
        A[i][n+1]=Y[i];
}
double** submatriz(double **matriz, int orden, int i, int j)
{

	double **subm;
	int p, q;				// Indices para la matriz
	int a = 0, b;		// Indices para la submatriz
	subm = new double* [orden - 1];

	for(p = 0; p < orden; p++) {
		if(p==i) continue;				//salta la i-esima fila
			subm[a] = new double[orden - 1];

			b = 0;

		for(q = 0; q< orden; q++) {
				if(q==j) continue;		//salta la j-esima columna
			subm[a][b] = matriz[p][q];
			b++;
		}
		a++; //Incrementa el indice de la fila
	}
	return subm;
}

double determinante(double **matriz, int orden)
{
	if(orden == 1)
		return matriz[0][0]; //Retorna el elemento  de la matriz si es de orden uno

	int i;
	int det = 0;
	for(i = 0; i < orden; i++)
		det += static_cast<int>(pow(-1.0,(int)i)) * matriz[i][0] * determinante(submatriz(matriz, orden, i, 0), orden - 1);
	return det;
}
void mostrar_matriz_A(double **A,int n)
{
     cout<<"\n\n A =\n";
    cout<<setw(5);
    for (int i=0;i<n;i++)            //print the new matrix
    {
        cout<<"|";
        for (int j=0;j<=n;j++){
            cout<<"\t"<<setw(8)<<A[i][j]<<setw(8);
            }cout<<"|";
 			cout<<"\n    ";
    }
}
void mostrar_Coef(double **A,int n)
{
     cout<<"\n\n A =\n";
    cout<<setw(5);
    for (int i=0;i<n;i++)            //print the new matrix
    {
        cout<<"|";
        for (int j=0;j<n;j++){
            cout<<"\t"<<setw(8)<<A[i][j]<<setw(8);
            }cout<<"|";
 			cout<<"\n    ";
    }
}
double **Transpuesta(double **A, int n){
   //1. Se crea la matriz donde se va a guardar la traspuesta.
   double**aux ;
   aux=new double* [MAX];
   for (int i=0;i<n;i++)
   {
       aux[i]=new double [n];
       for(int j=0;j<n;j++)
       {
           aux[i][j]=0;
       }
   }
   //2. Se guarda la traspuesta en aux.
   for (int i=0; i<n; i++)
   {
       for (int j=0; j<n; j++)
       {
           aux[j][i]=A[i][j];
       }
   }
   return aux;
}
bool Determinante_Cero(double **M, int n)
{
            if (determinante(M,n) != 0.000)
            {
                return false;
            }
            else
            {
              return true;
            }

}
bool Es_simetrica(double **M, int n)
{
    double **Tr=Transpuesta(M,n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (M[i][j] != Tr[i][j])
                return false;
    return true;
}
bool Pivotes_Cero(double **M, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (M[i][i] != 0)
                return false;
    return true;
}
void Doolitle(double **a, int n,double **l, double **u){//Factorización de A=LU donde L posee la diagonal de unos
int i,j,k;
double sum;
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			if(i>j)
			{
			  u[i][j]=0;
			}
			else if(i==j)
			{
			  l[i][j]=1;
			}
			else
			{
			  l[i][j]=0;
			}

	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			sum=0;
			if(i<=j)
			{
				for(k=0;k<n;k++)
					if(k!=i)
					sum=sum+l[i][k]*u[k][j];
				        u[i][j]=a[i][j]-sum;
			}
			else
			{
				for(k=0;k<n;k++)
					if(k!=j)
						sum=sum+l[i][k]*u[k][j];
					l[i][j]=(a[i][j]-sum)/u[j][j];
			}
		}
	}
}
void Crout (double **a, int n, double **l, double **u){
    int k=0;
    double sum=0;
    for(int i=0; i<n; i++){// 3. Llena las matrices. Llena una columna de la matriz L y luego una fila de la matriz U.
        sum=0;
            for(int k=0; k<n; k++){
                sum=0;
            for(int j=0; j<i; j++){// 3.1 Multiplica la fila de la matriz L con la columna de la matriz U.
                sum=sum+l[k][j]*u[j][i];
            }
            l[k][i]=a[k][i]-sum;//3.2 Toma el resultado de la multiplicación y le resta a matriz A.
        }//Llena las columnas de la matriz L.
        sum=0;
       for(int k=0; k<n; k++){//3.3  Multiplica la fila de la matriz L con la columna de la matriz U.
            sum=0;
            for(int j=0; j<i; j++){
                sum=sum+l[i][j]*u[j][k];//3.4 Toma el resultado de la multiplicación y le resta a matriz A, luego divide
                // entre una posición dada de la matriz L.
            }
            u[i][k]=(a[i][k]-sum)/l[i][i];
        }//Llena las filas de la matriz U.
    }
}
double **Cholesky (double **A, int n){
   //1. Se crea la matriz L y se inicializa con entradas aleatorias.
   double**L ;
   L=new double* [MAX];
   for (int i=0;i<n;i++)
   {
       L[i]=new double [n];
       for(int j=0;j<n;j++)
       {
           L[i][j]=0;
       }
   }

   //2. Se hace un recorrido primero por columnas y luego filas.
   for (int j=0; j<n; j++){
       for(int i=0; i<n; i++){

           //3. Se plantea el proceso para los elementos en la diagonal principal.
           if (i==j){
               double sum1=0;
               for(int k=0;k<i;k++)
               {
                   sum1+=pow(L[i][k], 2);
               }
               L[i][j]=sqrt(A[i][j]-sum1);
               sum1=0;
           }

           //4. Se plantea el proceso para los elementos en la diagonal principal.
           else if (i>j){
               double sum2=0;
               for (int k=0; k<j; k++){
                   sum2+=(L[i][k]*L[j][k]);
               }
               L[i][j] = A[i][j] - sum2;
               L[i][j]/= L[j][j];
               sum2=0;
           }

           //5. Se llena el triángulo superior de la matriz de ceros.
           else{
               L[i][j]=0;
           }
       }
   }

   //6. La función retorna la matriz L.
   return L;
}
void pivoteo_parcial(double **A, int n)
{
    for (int i=0;i<n;i++)
        for (int k=i+1;k<n;k++)
            if (abs(A[i][i])<abs(A[k][i]))
                for (int j=0;j<=n;j++)
                {
                    double aux=A[i][j];
                    A[i][j]=A[k][j];
                    A[k][j]=aux;
                }
}
void eliminacion_gaussiana(double **A, int n)
{
    for (int i=0;i<n-1;i++)            //loop to perform the gauss elimination
        for (int k=i+1;k<n;k++)
            {
                double m=A[k][i]/A[i][i];
                for (int j=0;j<=n;j++)
                    A[k][j]=A[k][j]-m*A[i][j];    //make the elements below the pivot elements equal to zero or elimnate the variables
            }
}
double **coeficientes(double **A, int n)
{
    double **a;
    a=new double*[n];
    for ( int i = 0 ; i < n ; i ++)
    {
        a[i]=new double[n];
        for (int j = 0 ; j < n ; j ++)
        {
        a[i][j]=A[i][j] ;
        }
    }
    return a;
}
double *terminos_independientes(double **A, int n)
{
    double *b=new double[n];
    for (int i=0; i<n; i++)
    {
        b[i] = A[i][n];
    }
 return b;
}
double *Sustitucion_Progresiva (double **A, double *b, int n)
{
   // 1. Define el vector solución de retorno de la función.
   double *y = new double[n];

   //2. Recorre la matriz por filas de arriba para abajo.
   for (int j=0; j<n; j++){

       //3. Paso base de la sustitución.
       if (j==0){
           y[j] = (b[j]/A[j][j]);
       }

       //4. Paso general de la sustitución.
       else{
           double sum=0;
           for ( int k=0; k<j; k++){
               //5. Realiza la suma de los productos coeficientes-incógnita de las incógnitas ya calculadas.
               sum+= y[k]*A[j][k];
           }

           //6. Resta al termino independiente la suma obtenida y la divide por el coeficiente correspondiente.
           y[j] = (b[j] - sum)/A[j][j];

       }

   }

   return y;
}
double *Sustitucion_Regresiva (double **M, double *Y, int n)
{
   // 1. Define el vector solución de retorno de la función.
   double *x = new double[n];

   //2. Recorre la matriz por filas de abajo para arriba.
   for (int j=n-1; j>=0; j--){

       //2.1. Paso base de la sustitución.
       if (j==n-1){
           x[j] = (Y[j]/M[j][j]);
       }

       //2.2. Paso general de la sustitución.
       else{
           double sum=0;
           //3. Realiza la suma de los productos coeficientes-incógnita de las incógnitas ya calculadas.
           for ( int k=n-1; k>j; k--){
               sum+= x[k]*M[j][k];
           }
           //6. Resta al termino independiente la suma obtenida y la divide por el coeficiente correspondiente.
           x[j] = (Y[j] - sum)/ M[j][j];
       }

   }

   return x;
}
double **Multiplicaion_Matriz(double** A,double** B,int Af,int Ac,int Bf,int Bc)
 {
    double **M;
    M=new double*[Af];
    for(int i=0;i<Af;i++)
    {
        M[i]=new double[Bc];
        for(int j=0;j<Bc;j++)
        {
            M[i][j]=0;
        }
    }
    for(int i=0;i<Af;i++)
    {
        for(int j=0;j<Bc;j++)
        {
            for(int k=0;k<Ac;k++)
            {
                M[i][j]+=A[i][k]*B[k][j];
            }
        }
    }
   return M;
 }
 void mostrar_matriz_L(double **L,int n)
 {
    cout<<"\n\n Matriz triangular inferior ";
	cout<<"\n\n L =\n";
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
            cout<<"\t"<<setw(6)<<L[i][j]<<setw(6);
 			cout<<"\n    ";
	}
 }
 void mostrar_matriz_Lt(double **Lt,int n)
 {
    cout<<"\n\n Matriz triangular inferior transpuesta ";
	cout<<"\n\n Lt =\n";
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
            cout<<"\t"<<setw(6)<<Lt[i][j]<<setw(6);
 			cout<<"\n    ";
	}
 }
 void mostrar_matriz_U(double **U,int n)
 {
    cout<<"\n\n ";
    cout<<"\n\n Matriz triangular superior ";
	cout<<"\n\n U =\n";
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
            cout<<"\t"<<setw(6)<<U[i][j]<<setw(6);
 			cout<<"\n    ";
	}
 }
 void mostrar_solucion_Y(double *Y,int n)
 {
    cout<<"\n\n La solucion intermedia \n";
	cout<<" con Y : \n";

	for(int i=0;i<n;i++)
                cout<<"\n y"<<"["<<i+1<<"]="<<"\t"<<setw(6)<<Y[i];
                cout<<"\n\n";
 }
void mostrar_solucion_X(double *X,int n)
 {
    cout<<"\n\n La solucion final \n";
	cout<<" con X :\n";

	for(int i=0;i<n;i++)
    cout<<"\n x"<<"["<<i+1<<"]="<<"\t"<<setw(6)<<X[i];
	cout<<"\n\n";
 }

double error_ajuste(double *x, double *y, double *a, int N, int n)
{

    double sum =0;
    for (int i = 0; i < N; i++) {
        double yEsperado = 0;
        for (int j = 0; j <= n; j++)
            yEsperado += a[j]* pow(x[i],j);
        double diff = yEsperado - y[i];
        sum += pow( diff, 2);
    }
    sum=sqrt(sum);
   return sum;
}
double error_cuadratico_medio(double *x, double *y, double *a, int N, int n)
{

    double sum =0;
    for (int i = 0; i < N; i++) {
        double yEsperado = 0;
        for (int j = 0; j <= n; j++)
            yEsperado += a[j]* pow(x[i],j);
        double diff = yEsperado - y[i];
        sum += pow( diff, 2);
    }
    sum=sqrt(sum/N);
   return sum;
}
void mostrar_parametros_hallados(double *c, int n, int g)
{
     cout<<"\n los parametros del polinomio de grado "<<g<< " son:\n\n";
    for (int i=0;i<n;i++)
        cout<<"x^"<<i<<"="<<c[i]<<endl;            // Muestra las variables x^0,x^1,x^2,x^3,...
}
void mostrar_polinomio_ajustado(double *c,int n)
{
    cout<<"\n El polinomio que mejor se aproxima al conjunto de datos es:\n\n y=";
    for (int i=0;i<n;i++)
        cout<<" + ("<<c[i]<<")"<<"x^"<<i;
    cout<<"\n";
}
void mostrar_error_ajuste(double *x, double *y, double *c,int N, int n)
{
   cout<<"\n El error del ajuste es:\n";
    cout<<"e= "<<error_ajuste(x,y,c,N,n);
}
 void mostrar_error_cuadratico_medio(double *x, double *y, double *c, int N, int n)
{
   cout<<"\n El error cuadratico medio:\n";
   cout<<"e= "<<error_cuadratico_medio(x,y,c,N,n);
}
int main()
{
 int op1,op2;

 char centinela;
 bool salir1=false;
 bool salir2=false;
 int n,g,N;

    cout.precision(3);
    cout.setf(ios::fixed);

        do
      {
          titulo();
    cout<<"\t Ingrese la cantidad de pares de datos que dispone\n";
    cout<<"\t ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ :\n ";        //Establece el tamaño de del arreglo matricial sobre el cual estara el conjunto de datos
    cin>>N;
    double *x=new double[N];
    double *y=new double[N];
    lectura_datos(x,y,N);
    cout<<"\n\t"<<char(168)<<" Cual es el grado del polinomio sobre el que desea ajustar los datos"<<char(63)<<'\n';
    cout<<"\t Ingrese el grado del polinomio\n";
    cout<<"\t ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ :\n ";
    cin>>n;
    g=n;
    system("cls");
    titulo();                        // n es el grado del polinomio
    double *X=new double[2*n+1];
    SumatoriasXin(X,x,N,n);
    double **A=arreglo_matricial_SXin(A,X,n);
    double *Y=new double[n+1];
    SumatoriasYiXin(Y,x,y,N,n);
    arreglo_matricialAumentado_SYiXin(A,Y,n);
    n=n+1;
    double **k=coeficientes(A,n);
    if(!Determinante_Cero(k,n))
    {
         op1=Menu(op1);


      switch(op1)
            {
             case 1 :
                      do
                    {
                        system("cls");
                        op2=SubMenu(op2);


                        switch(op2)
                    {
                    case 1 :
                        {
                             cout<<"\t\t\t ÉÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ»\n";
                             cout<<"\t\t\t º                  SOLUCION DE SISTEMAS LINEALES AX=B            º\n";
                             cout<<"\t\t\t º                  POR FACTORIZACION LU DE DOOLITTLE             º\n";
                             cout<<"\t\t\t ÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼\n\n";

                     cout<<"\n\n Matriz aumentada:\n";
                     mostrar_matriz_A(A,n);
                     double **a=coeficientes(A,n);
                     double *b=terminos_independientes(A,n);
                     double**L1 ;
                    L1=new double* [MAX];
                    for (int i=0;i<n;i++)
                    {
                    L1[i]=new double [n];
                    }
                    double**U1 ;
                    U1=new double* [MAX];
                    for (int i=0;i<n;i++)
                    {
                    U1[i]=new double [n];
                    }
                     Doolitle(a,n,L1,U1);
                      if(!Pivotes_Cero(U1,n))
                    {
                     double* q=Sustitucion_Progresiva(L1,b,n);
                     double* p=Sustitucion_Regresiva(U1,q,n);
                     mostrar_matriz_L(L1,n);
                     mostrar_matriz_U(U1,n);
                     mostrar_solucion_Y(q,n);
                     mostrar_solucion_X(p,n);
                     mostrar_parametros_hallados(p,n,g);
                     mostrar_polinomio_ajustado(p,n);
                     mostrar_error_ajuste(x,y,p,N,n);
                     mostrar_error_cuadratico_medio(x,y,p,N,n);
                    }
                    else
                        {
                        cout<<"\n Alguno de los elementos obtenidos sobre la diagonal de la matriz triangular superior \n";
                        cout<<"\n son  ceros :\n";

                        cout<<"\n\n No es posible Resolver el sistema por este metodo de factorizacion:\n";
                        cout<<"\nIntentelo de nuevo!!\n";
                        }
                        }



                    break;

                    case 2 :
                        {
                             cout<<"\t\t\t ÉÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ»\n";
                             cout<<"\t\t\t º                  SOLUCION DE SISTEMAS LINEALES AX=B            º\n";
                             cout<<"\t\t\t º                  POR FACTORIZACION LU DE CROUT                 º\n";
                             cout<<"\t\t\t ÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼\n\n";
                    cout<<"\n\n Matriz aumentada:\n";
                    mostrar_matriz_A(A,n);
                    double **a=coeficientes(A,n);
                    double *b=terminos_independientes(A,n);
                    double**L2 ;
                    L2=new double* [MAX];
                    for (int i=0;i<n;i++)
                    {
                    L2[i]=new double [n];
                    }
                    double**U2 ;
                    U2=new double* [MAX];
                    for (int i=0;i<n;i++)
                    {
                    U2[i]=new double [n];
                    }
                    Crout(a,n,L2,U2);
                    if(!Pivotes_Cero(L2,n))
                    {
                    double* q=Sustitucion_Progresiva(L2,b,n);
                    double* p=Sustitucion_Regresiva(U2,q,n);
                    mostrar_matriz_L(L2,n);
                    mostrar_matriz_U(U2,n);
                    mostrar_solucion_Y(q,n);
                    mostrar_solucion_X(p,n);
                    mostrar_parametros_hallados(p,n,g);
                    mostrar_polinomio_ajustado(p,n);
                    mostrar_error_ajuste(x,y,p,N,n);
                    mostrar_error_cuadratico_medio(x,y,p,N,n);
                    }
                    else
                        {
                        cout<<"\n Alguno de los elementos obtenidos sobre la diagonal de la matriz triangular inferior \n";
                        cout<<"\n son  ceros :\n";

                        cout<<"\n\n No es posible Resolver el sistema por este metodo de factorizacion:\n";
                        cout<<"\nIntentelo de nuevo!!\n";
                        }
                        }



                    break;
                    case 3 :
                        {
                             cout<<"\t\t\t ÉÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ»\n";
                             cout<<"\t\t\t º                  SOLUCION DE SISTEMAS LINEALES AX=B            º\n";
                             cout<<"\t\t\t º                  POR FACTORIZACION LU DE CHOLESKI              º\n";
                             cout<<"\t\t\t ÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼\n\n";
                             double **a=coeficientes(A,n);
                             double *b=terminos_independientes(A,n);
                            if(Es_simetrica(a,n))
                                {
                     cout<<"\nla matriz si es simetrica\n";
                     cout<<"\n\n Matriz de coeficientes :\n";
                     mostrar_Coef(a,n);
                     cout<<"\n";
                     cout<<"\n\n Matriz aumentada:\n";
                     mostrar_matriz_A(A,n);
                     double** L=Cholesky(a,n);
                     if(!Pivotes_Cero(L,n))
                    {
                     double** Lt=Transpuesta(L,n);
                     double* q=Sustitucion_Progresiva(L,b,n);
                     double* p=Sustitucion_Regresiva(Lt,q,n);
                     mostrar_matriz_L(L,n);
                     mostrar_matriz_Lt(Lt,n);
                     mostrar_solucion_Y(q,n);
                     mostrar_solucion_X(p,n);
                     mostrar_parametros_hallados(p,n,g);
                     mostrar_polinomio_ajustado(p,n);
                     mostrar_error_ajuste(x,y,p,N,n);
                     mostrar_error_cuadratico_medio(x,y,p,N,n);
                    }
                    else
                    {
                       cout<<"\n Alguno de los elementos obtenidos sobre la diagonal de la matriz triangular inferior \n";
                        cout<<"\n son  ceros :\n";

                        cout<<"\n\n No es posible Resolver el sistema por este metodo de factorizacion:\n";
                        cout<<"\nIntentelo de nuevo!!\n";
                    }
                        }
                        else
                        {
                        cout<<"\nla matriz no es simetrica\n";
                        cout<<"\n\n Matriz de coeficientes :\n";
                        mostrar_Coef(a,n);
                        cout<<"\n\n No es posible Resolver el sistema por este metodo de factorizacion:\n";
                        cout<<"\nIntentelo de nuevo!!\n";
                        }
                                }

                    break;



                    case 4 :
                     salir2 = true;
                    }

                    if(!salir2)
                    {
                    cout<<"\n\n"<<char(168)<<" Desea continuar "<<char(63)<<"\n";
                    cout<<"\t presione c o C para continuar o cualquier otra tecla para salir \n ";
                    cout<<"\t ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ:\n ";
                    cin>>centinela;
                    }

                    system("cls");
                    } while((centinela=='c'||centinela=='C')&&!salir2);

             case 2 :
                 {
                             cout<<"\t\t\t ÉÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ»\n";
                             cout<<"\t\t\t º                  SOLUCION DE SISTEMAS LINEALES AX=B                º\n";
                             cout<<"\t\t\t º             POR ELIMINACION GAUSSIANA CON PIVOTEO PARCIAL          º\n";
                             cout<<"\t\t\t ÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼\n\n";
                        cout<<"                                 SOLUCION DE SISTEMAS LINEALES AX=B                           "<<"\n";
                        cout<<"                           POR ELIMINACION GAUSSIANA CON PIVOTEO PARCIAL                         "<<"\n\n";
                        cout<<"\n Matriz aumentada:\n";
                        mostrar_matriz_A(A,n);
                        pivoteo_parcial(A,n);
                        cout<<"\nPivotasion parcial de las filas de la Matriz es:\n";
                        mostrar_matriz_A(A,n);
                        eliminacion_gaussiana(A,n);
                        cout<<"\n\nMatriz despues de la eliminacion gaussiana:\n";
                        mostrar_matriz_A(A,n);
                        double **a=coeficientes(A,n);
                        double *b=terminos_independientes(A,n);
                        double *c=Sustitucion_Regresiva(a,b,n);
                        mostrar_parametros_hallados(c,n,g);
                        mostrar_polinomio_ajustado(c,n);
                        mostrar_error_ajuste(x,y,c,N,n);
                        mostrar_error_cuadratico_medio(x,y,c,N,n);
                 }



                    break;



             case 3 :
                     salir1 = true;
             }

       if(!salir1)
         {
          cout<<"\n\n"<<char(168)<<" Desea continuar "<<char(63)<<"\n";
          cout<<"\t presione c o C para continuar o cualquier otra tecla para salir \n ";
          cout<<"\t ÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ:\n ";
          cin>>centinela;
         }

         system("cls");
     }
     else
     {
        cout<<"\t La matriz  A de nxn cumple las condiciones para un sistema incompatible: son\n";
        cout<<"\t 1. El determinante de A = "<<determinante(k,n)<<"\n";
        cout<<"\t 2. La matriz A es singular \n";
        cout<<"\t 3. La matriz no se reduce por filas a la matriz identidad In \n";
        break;
     }
      } while((centinela=='c'||centinela=='C')&&!salir1);
}

