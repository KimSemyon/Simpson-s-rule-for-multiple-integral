#include <iostream>
#include <thread>
#include <queue>
#include <vector>
#include <condition_variable>
#include <mutex>
#include <cmath>
#include <ctime>
#include <chrono>
#include <atomic>

using namespace std;

double epsilon = 1e-15;// значение близкое к нулю
int deep_recursion; // inversion of deep recursion
atomic<int> count_mainthread = 0; // for pause main thread when main thread wait for other threads
class ThreadPool
{
public:

	ThreadPool(int threads) : shutdown_(false)
	{
		// Create the specified number of threads
		threads_.reserve(threads);
		for (int i = 0; i < threads; ++i)
			threads_.emplace_back(bind(&ThreadPool::threadEntry, this, i));
	}

	~ThreadPool()
	{
		{
			// Unblock any threads and tell them to stop
			unique_lock <mutex> l(lock_);

			shutdown_ = true;
			condVar_.notify_all();
		}

		// Wait for all threads to stop
		cerr << "Joining threads" << endl;
		for (auto& thread : threads_)
			thread.join();
	}

	void doJob(function <void(void)> func)
	{
		// Place a job on the queu and unblock a thread
		unique_lock <mutex> l(lock_);

		jobs_.emplace(move(func));
		condVar_.notify_one();
	}
protected:

	void threadEntry(int i)
	{
		function <void(void)> job;

		while (1)
		{
			{
				unique_lock <mutex> l(lock_);

				while (!shutdown_ && jobs_.empty())
					condVar_.wait(l);

				if (jobs_.empty())
				{
					// No jobs to do and we are shutting down
					cerr << "Thread " << i << " terminates" << endl;
					return;
				}

				//std::cerr << "Thread " << i << " does a job" << std::endl;
				job = move(jobs_.front());
				jobs_.pop();
			}

			// Do the job without holding any locks
			job();
		}

	}

	mutex lock_;
	condition_variable condVar_;
	bool shutdown_;
	queue <function <void(void)>> jobs_;
	vector <thread> threads_;
};

double Nfunc(int N, double X[]){
	double Sum = 0;
	for (int i = 0; i < N; i++){ Sum += X[i];}
	//Sum = sin(Sum);
	return Sum;
	//return sin(X[0]);
}

double NfuncP(int N, double X[], double last_elem){
	double Sum = 0;
	for (int i = 0; i < N - 1; i++){ Sum += X[i]; }
	Sum += last_elem;
	//Sum = sin(Sum);
	return Sum;
	//return sin(last_elem);
}

void Simpf(double a_, double h_, int n, int N, double X[], std::vector <double> &F, int part_from, int part_to){
	for (double step = a_; part_from < part_to; step += h_)
	{
		F[part_from] = NfuncP(N, X, step);
		++part_from;
	}
	++count_mainthread;
}

double SerialSimpsonsFunc(double Segments[], double H[], int n, int N, int _param, double* X, std::vector <double> &F)
{
	deep_recursion--;
	double sum_odd = 0; 
	double sum_even = 0;
	long double SumOfIntegral=0;
	int param = _param;
	double a = Segments[2*param];
	//double b = Segments[2*param+1];
	double h = H[param];
	int count_i = 0;
	int count_res = 2 * n*deep_recursion; // счетчик для "external" заполнения значений функций
	for (double i = a; count_i<2*n; i += h, count_i++)
	{
		if (param < (N - 1)) { 
			X[param] = i; 
			param++;
			F[count_res] = SerialSimpsonsFunc(Segments, H, n, N, param, X, F); //Заполнение F 2n позиций 
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
//-----------------------------------------------------------------------------------------
	
	deep_recursion++;

	return SumOfIntegral;
}

double ParallelSimpsonsFunc(ThreadPool &p, double Segments[], double H[], int n, int N, int _param, double* X, std::vector <double> &F, int nthreads)
{
	deep_recursion--;
	double sum_odd = 0;
	double sum_even = 0;
	double SumOfIntegral = 0;
	int param = _param;
	double a = Segments[2 * param];
	//double b = Segments[2 * param + 1];
	double h = H[param];
	int count_res = 2 * n*deep_recursion; // счетчик для "external" заполнения значений функций
	for (double i = a, count_i = 0; count_i < 2 * n; i += h, count_i++)
	{
		if (param < (N - 1)) {
			X[param] = i;
			param++;
			F[count_res] = ParallelSimpsonsFunc(p, Segments, H, n, N, param, X, F, nthreads); //Заполнение F 2n позиций 
			count_res++;
			param--;
		} //
		else {
			for (int t = 0; t < nthreads; t++){
				double a_ = (t * n / nthreads) * h;
				//double b_ = ((t + 1) * 2 * n / nthreads - 1) * h;
				p.doJob(std::bind(Simpf, a_, h, n, N, X, std::ref(F), t * n / nthreads, (t + 1) * n / nthreads));
			}

			Simpf(n * h, h, n, N, X, F, n, 2 * n); // main thread works

			while (true){// main thread wait for other threads when queue jobs_ is temporarily empty
			if (count_mainthread == nthreads+1) { count_mainthread = 0; break; }
			}
			break;
		}
	}
	
	
	//------------------------------------------------------------------Вычисляем внешнюю сумму	
	count_res = 2 * n*deep_recursion;

	SumOfIntegral += F[count_res];
	SumOfIntegral += F[count_res + 2 * n - 1];
	for (int j = count_res + 1; j < (count_res + 2 * n - 1); j++){
		if (j % 2 == 0) { sum_even += F[j]; }
		else { sum_odd += F[j]; }
	}
	SumOfIntegral += 4 * sum_odd + 2 * sum_even;
	SumOfIntegral = SumOfIntegral*h / 3;
	
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
	

	int nthreads = 1;//std::thread::hardware_concurrency()-1; // default = 2 
	cout << " number threads " << nthreads << endl;
	//cout << "Enter the number of threads: ";
	//cin >> nthreads;
	
	double *H = new double[N]; //Вычислили значения шагов
	for (int i = 0; i < N; i++){
		H[i] = (Segments[2 * i + 1] - Segments[2 * i]) / (2 * n - 1);
	}

	double *X = new double[N];
	
	std::vector <double> F(N*2*n);

	
////////////////////////////////////////////////////////////////////////////////////SERIAL:	
	//while (1) {
	//
	
	/*deep_recursion = N;
	int param = 0;
	
	unsigned int start_time = clock();
	chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
	double result = SerialSimpsonsFunc(Segments, H, n, N, param, X, F);
	chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
	unsigned int end_time = clock();
	cout << "Result: " << result << endl;
	unsigned int search_time = end_time - start_time;
	cout << "time(ALL): " << search_time/ CLOCKS_PER_SEC << " seconds" << endl;
	chrono::duration<double, std::milli> time_span = t2 - t1;
	std::cout << "It took me(ALL) " << time_span.count() << " milliseconds." << endl;
	*/

	//
	//if (fabs(REAL_result - result) < epsilon) { printf("Решение получено"); break;}
	//
	//func_change_n(&n);				//void func_change_n(int *n) {*n *= 2;}
	//}

//////////////////////////////////////////////////////////////////////////////////////PARALLEL:
	
	deep_recursion = N;
	int paramforparallel = 0;
	//
	ThreadPool p(nthreads);
	unsigned int start_timeforparallel = clock();
	
	chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
	
	double resultforparallel = ParallelSimpsonsFunc(p, Segments, H, n, N, paramforparallel, X, F, nthreads);
	
	chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();

	unsigned int end_timeforparallel = clock();
	cout << "Result for parallel: " << resultforparallel << endl;
	unsigned int search_timeforparallel = end_timeforparallel - start_timeforparallel;
	cout << "time(ALL): " << search_timeforparallel / CLOCKS_PER_SEC << " seconds" << endl;;
	
	chrono::duration<double, std::milli> time_span = t2 - t1;
	std::cout << "It took me(ALL) " << time_span.count() << " milliseconds." << endl;

	chrono::duration<double> time_span_seconds = t2 - t1;
	std::cout << "It took me(ALL) " << time_span_seconds.count() << " seconds." << endl;
	

	delete[] Segments;
	delete[] H;
	delete[] X;
	system("PAUSE");	
	return 0;
}