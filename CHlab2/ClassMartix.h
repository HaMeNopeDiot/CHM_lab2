#pragma once
#include "ClassVector.h"
#include <ctime>
#include <iomanip>
#include <vector>

class ClassMatrix : ClassVector
{

	int L, N;
protected:
	vector<ClassVector> vectors;
	ClassVector f, x;  

public:
	ClassMatrix(int _N, int _L, int _start = 1, int _end = 10);
	ClassMatrix(vector<ClassVector> _vectors, ClassVector _x);

	void print_class_matrix_vectors();
	void print_class_matrix();
	vector<long double> calculate_f();

	long double get_cell(int i, int j);
	void set_cell(int i, int j, long double digit);
	void solve_ext();
	long double calculate_falibility(ClassVector x, ClassVector xot);
	void delim_str(int index, long double delim);
};