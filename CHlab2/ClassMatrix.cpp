#include <vector>
#include <iomanip>
#include "ClassMartix.h"

ClassMatrix::ClassMatrix(int _N, int _L, int _start, int _end)
{
	N = _N;
	L = _L;
	for (int i = 1; i <= (2 * L - 1); i++) {
		vector<long double> vec;
		ClassVector tmp_vec;
		int distance_from_the_center = L;
		
		int tmp = L - i;
		int k;
		for (k = 0; k < tmp; k++) { // tmp > 0
			vec.push_back(0);
		}
		for (k; k < N + tmp && k < N; k++) {
			vec.push_back(_start + rand() % _end);
		}
		for (k; k < N; k++) {       // tmp < 0
			vec.push_back(0);
		}

		tmp_vec.change_vector(vec);
		vectors.push_back(tmp_vec);
		
	}
	vector<long double> tmp_x;
	vector<long double> tmp_f(N, 0);
	for (int i = 0; i < N; i++) {
		tmp_x.push_back(_start + rand() % _end);
	}
	x.change_vector(tmp_x);
	f.change_vector(calculate_f());
}

ClassMatrix::ClassMatrix(vector<ClassVector> _vectors, ClassVector _x)
{
	L = (_vectors.size() + 1) / 2;
	N = _vectors[0].get_size();
	x = _x;
	vectors = _vectors;
	f = calculate_f();
	
}

void ClassMatrix::print_class_matrix_vectors()
{
	for (int i = 1; i <= N; i++) {
		cout << "i:" << i << " ";
		for (int j = 0; j < vectors.size(); j++) {
			cout << setw(2) << vectors[j].get_elem_vector(i) << " ";
		}
		cout << "|" << setw(2) << x.get_elem_vector(i) << setw(2) << "|"  << setw(2) << f.get_elem_vector(i);
		cout << '\n';
	}
}

void ClassMatrix::print_class_matrix()
{
	for (int i = 1; i <= N; i++) {
		cout << "i:" << i << "| ";
		for (int j = 1; j <= N; j++) {
			cout << setw(8) << get_cell(i, j) << " ";
		}
		cout << '\n';
	}
}

vector<long double> ClassMatrix::calculate_f()
{
	vector<long double> _f(N, 0);
	// L -> (2 * L) - 1;
	// int max_var = (2 * L) - 1;
	// int cnt_var = L;
	int i = 1;
	for (i; i <= N; i++) {
		int j = i - L + 1 > 0 ? i - L + 1: 1;
		int cnt = L - i > 0 ? L - i: 0;
		for (j; j < i; j++) {
			_f[i - 1] += x.get_elem_vector(j) * vectors[cnt++].get_elem_vector(i);
		}
		for (j; j < i + L && j <= N; j++) {
			_f[i - 1] += x.get_elem_vector(j) * vectors[cnt++].get_elem_vector(i);
		}
	}
	return _f;
	
	/*for (i; i <= N - L + 1; i++) {

	}
	for (i; i <= N; i++) {

	}*/
}

void ClassMatrix::solve_ext()
{
	vector<ClassVector> tmp = vectors;
	
	for (int i = 2; i <= N; i++) {
		set_cell(1, i, (get_cell(1, i) * 1.0f) / get_cell(1, 1));
	}
	for (int i = 2; i <= N; i++) {
		// b (low matrix)
		for (int j = 2; j <= i; j++) {
			long double tmp = get_cell(i, j);

			int k = i > L ? i - L + 1 : 1;
			for (k; k <= j - 1; k++) {
				tmp -= get_cell(i, k) * get_cell(k, j);
			}
			set_cell(i, j, tmp);
		}
		//vectors[2].print_vector();
		// c (high matrix)
		int ova = i > N - L ? N : i + L - 1;
		for (int j = i + 1; j <= ova; j++) {
			long double tmp = get_cell(i, j);
			
			int k = j > L ? j - L + 1 : 1;
			for (k; k <= i - 1; k++) {
				tmp -= get_cell(i, k) * get_cell(k, j);
			}
			tmp = tmp * 1.0f / get_cell(i, i);
			set_cell(i, j, tmp);
			
		}
	}
	vector<ClassVector> low, high;
	ClassVector zeroes(N, (long double)0.0);
	ClassVector ones(N, (long double)1.0);
	for (int i = 1; i < L; i++) {
		low.push_back(vectors[i - 1]);
		high.push_back(zeroes);
	}
	low.push_back(vectors[L - 1]);
	high.push_back(ones);
	for (int i = L + 1; i <= (2 * L) - 1; i++) {
		low.push_back(zeroes);
		high.push_back(vectors[i - 1]);
	}

	ClassMatrix Lower_Matrix(low, x);
	ClassMatrix Upper_Matrix(high, x);

	cout << "-=-=-=-=-=-=-=-=- L -=-=-=-=-=-=-=-=-=-=-\n";
	Lower_Matrix.print_class_matrix();
	cout << "-=-=-=-=-=-=-=-=- U -=-=-=-=-=-=-=-=-=-=-\n";
	Upper_Matrix.print_class_matrix();
	cout << "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n";


	ClassVector y(N, (long double)0.0);
	for (int i = 1; i <= N; i++) {
		long double tmp = f.get_elem_vector(i);
		int k = i > L ? i - L + 1 : 1;
		for (k; k <= i - 1; k++) {
			tmp -= Lower_Matrix.get_cell(i, k) * y.get_elem_vector(k);
		}
		y.set_elem_vector(i, tmp * 1.0f / Lower_Matrix.get_cell(i, i));
	}

	ClassVector x_new(N, (long double)0.0);
	for (int i = N; i >= 1; i--) {
		long double tmp = y.get_elem_vector(i);
		int ova = i > N - L ? N : i + L - 1;
		for (int k = i + 1; k <= ova; k++) {
			tmp -= Upper_Matrix.get_cell(i, k) * x_new.get_elem_vector(k);
		}
		x_new.set_elem_vector(i, tmp);
	}

	x.print_vector();
	x_new.print_vector();
	y.print_vector();

	std::cout << "Fallibility to X - Xot = " << setprecision(60) << calculate_falibility(x, x_new) << setprecision(4) << '\n';
	
	
}

long double ClassMatrix::calculate_falibility(ClassVector x, ClassVector xot)
{
	ClassVector fallib(N);
	for (int i = 1; i <= N; i++) {
		fallib.set_elem_vector(i, abs(x.get_elem_vector(i) - xot.get_elem_vector(i)));
	}
	return fallib.get_norm();
}

long double ClassMatrix::get_cell(int i, int j)
{
	if (abs(i - j) >= L) {
		return 0.0;
	} else {  
		// L - 1 -> L - 1 - 1 ... 0 -> 1 -> L - 1
		int num_vec = L + j - i - 1;
		return vectors[num_vec].get_elem_vector(i);
	}
}

void ClassMatrix::set_cell(int i, int j, long double digit)
{
	if (abs(i - j) < L) {
		int num_vec = L + j - i - 1;
		vectors[num_vec].set_elem_vector(i, digit);
	}
}

void ClassMatrix::delim_str(int index, long double delim)
{
	for (int i = 1; i <= N; i++) {
		set_cell(index, i, get_cell(index, i) * 1.0f / delim);
	}
}