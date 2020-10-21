#include <iostream>
#define _USE_MATH_DEFINES
#define _USE_GENERIC_MATH1
//#include <math.h>
#include<cmath>
using namespace std;

double epsilon = 1.0E-15;// значение близкое к нулю
int deep_recursion; // inversion of deep recursion


double func(double x, double param) 
{
	return sin(param*x);
}

double Nfunc(int N, double X[]){
	double Sum = 0;
	for (int i = 0; i < N; i++){ Sum += X[i];}
	//Sum = sin(Sum);
	return Sum;
}

double SimpsonsFunc(double Segments[], double H[], int n, int N, int _param, double* X, double* F) 
{
	deep_recursion--;
	double sum_odd = 0; 
	double sum_even = 0;
	double SumOfIntegral=0;
	int param = _param;
	double a = Segments[2*param];
	double b = Segments[2*param+1];
	double h = H[param];
	//int count = 0; // счетчик для присвоения значений функций
	int count_res = 2 * n*deep_recursion; // счетчик для "external" заполнения значений функций
	for (double i = a, count_i=0; (i <= b) || (count_i<2*n); i += h, count_i++)
	{
		if (param < (N - 1)) { 
			X[param] = i; 
			param++;
			F[count_res] = SimpsonsFunc(Segments, H, n, N, param, X, F); //Заполнение F 2n позиций 
			count_res++;
			param--; 
		} //
		else {
			X[param] = i;
			F[count_res] = Nfunc(N, X);
			count_res++;
		}
	}
//------------------------------------------------------------------Вычисляем внешнюю сумму	
	count_res = 2 * n*deep_recursion;

	SumOfIntegral += F[count_res];
	SumOfIntegral += F[count_res+2*n-1];
	for (int j = count_res + 1; j < (count_res + 2 * n - 1); j++){
		if (j % 2 == 0) { sum_even += F[j]; }
		else { sum_odd += F[j]; }
	}
	SumOfIntegral += 4 * sum_odd + 2 * sum_even;
	SumOfIntegral = SumOfIntegral*h / 3;
	//}

	//param--;
//-----------------------------------------------------------------------------------------
	
	deep_recursion++;

	return SumOfIntegral;
}



int main() {
		
	int N; // размерность интеграла
	int n; // количество точек на отрезке /2
		
	cout << "Enter the dimension of the integral: ";
	cin >> N;

	double *Segments = new double[2 * N]; // Выделение памяти для массива пределов интегрирования
	
	cout << "Enter the number of points on the segment: ";
	cin >> n;

	cout << "Enter pairs of integration limits starting from internal in the amount of " << N << ": ";
	
	for (int i = 0; i < 2*N; i++) { // Заполнение массива отрезков 
		cin >> Segments[i];
	}
	
	double *H = new double[N]; //Вычислили значения шагов
	for (int i = 0; i < N; i++){
		H[i] = (Segments[2 * i + 1] - Segments[2 * i]) / (2 * n - 1);
	}

	double *X = new double[N];
	
	double *F = new double[N*2*n]; 
	
	//while (1) {
	//
	deep_recursion = N;
	int param = 0;
	//
	double result = SimpsonsFunc(Segments, H, n, N, param, X, F);
	cout << "Result: " << result;
	//
	//if (fabs(REAL_result - result) < epsilon) { printf("Решение получено"); break;}
	//
	//func_change_n(&n);				//void func_change_n(int *n) {*n *= 2;}
	//}

	delete[] Segments;
	delete[] H;
	delete[] X;
	delete[] F;
	system("PAUSE");	
	return 0;
}