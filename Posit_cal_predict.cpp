#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "myfdlibm.h"

using namespace std;

void load_data(string file_name, vector<vector<Posit64>> &, vector<Posit64> &);
vector<vector<Posit64>> transpose(const vector<vector<Posit64>>);
void d1_normalization(vector<Posit64> &);
void d2_normalization(vector<vector<Posit64>> &);
vector<vector<Posit64>> multiplication(const vector<vector<Posit64>>, const vector<vector<Posit64>>);
vector<Posit64> multiplication(const vector<vector<Posit64>> A, const vector<Posit64> b);
vector<Posit64> gaussian_elimination(vector<vector<Posit64>> &A, vector<Posit64> &b);
Posit64 dot_product(vector<Posit64>, vector<Posit64>);
void caculate(vector<vector<Posit64>> data, vector<Posit64> x, vector<Posit64> true_val);

// g++ -std=c++17 -o Posit_cal_predict Posit_cal_predict.cpp Posit_fabs.cpp Posit_sqrt.cpp

int main(void)
{
	vector<vector<Posit64>> train_A;
	vector<Posit64> train_b;
	load_data("cal_train.txt", train_A, train_b);
	
	vector<vector<Posit64>> test_A;
	vector<Posit64> test_b;
	load_data("cal_test.txt", test_A, test_b);
	
	vector<vector<Posit64>> train_A_t = transpose(train_A);
	vector<vector<Posit64>> test_A_t  = transpose(test_A);
	
	// normalization
	d2_normalization(train_A_t);
	d2_normalization(test_A_t);
	train_A = transpose(train_A_t);
	test_A  = transpose(test_A_t);
	d1_normalization(train_b);
	d1_normalization(test_b);
	
	vector<vector<Posit64>> ATA = multiplication(train_A_t, train_A);
	vector<Posit64> ATB = multiplication(train_A_t, train_b);
	
	vector<Posit64> x = gaussian_elimination(ATA, ATB);
	caculate(test_A, x, test_b);
	
	
	
	return 0;
}

void d1_normalization(vector<Posit64> &src)
{
	int len = src.size();
	Posit64 min = src[0], max = src[0], diff{0.0};
	
	for(int i = 1; i < len; i++)
	{
		if(src[i] > max)
			max = src[i];
		if(src[i] < min)
			min = src[i];
	}
	
	diff = max - min;
	for(int i = 0; i < len; i++)
	{
		src[i] = (src[i] - min) / diff;
	}
}

void d2_normalization(vector<vector<Posit64>> &src)
{
	int len = src.size();
	
	for(int i = 0; i < len; i++)
	{
		d1_normalization(src[i]);
	}
}

void load_data(string file_name, vector<vector<Posit64>> &store_arr, vector<Posit64> &target)
{
	ifstream train_file(file_name, ios::in);
	if(!train_file)
	{
		cout << "open file error." << endl;
		exit(1);
	}
	
	string line;
	getline(train_file, line);
	vector<vector<Posit64>> train_arr;
	while(getline(train_file, line))
	{
		stringstream ss(line);
		string token;
		vector<Posit64> row;
		int col_count = 0;
		
		while(getline(ss, token, ','))
		{
			if(col_count >= 9)
				break;
				
			Posit64 value{0.0};
			if(!token.empty())
			{
				try
				{
					value = Posit64{stod(token)};
				}
				catch(...)
				{
					value = Posit64{0.0};
				}
			}
			
			if(col_count == 8)
				target.push_back(value);
			else
				row.push_back(value);
			col_count++;
		}
		
		store_arr.push_back(row);
	}
}

vector<vector<Posit64>> transpose(const vector<vector<Posit64>> src)
{
	if(src.empty())
		return {};
		
	size_t row = src.size();
	size_t col = src[0].size();
	
	vector<vector<Posit64>> res(col, vector<Posit64>(row));
	
	for(size_t i = 0; i < row; i++)
	{
		for(size_t j = 0; j < col; j++)
		{
			res[j][i] = src[i][j];
		}
	}
	
	return res;
}

vector<vector<Posit64>> multiplication(const vector<vector<Posit64>> a, const vector<vector<Posit64>> b)
{
	int row_a = a.size();
	int col_a = a[0].size();
	int col_b = b[0].size();
	
	vector<vector<Posit64>> res(row_a, vector<Posit64>(col_b, Posit64{0.0}));
	
	for(int i = 0; i < row_a; i++)
	{
		for(int j = 0; j < col_b; j++)
		{
			for(int k = 0; k < col_a; k++)
			{
				res[i][j] += a[i][k] * b[k][j];
			}
		}
	}
	
	return res;
}

vector<Posit64> multiplication(const vector<vector<Posit64>> A, const vector<Posit64> b)
{
    int row = A.size();
    int col = A[0].size();
    vector<Posit64> res(row, 0.0);

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            res[i] += A[i][j] * b[j];
        }
    }
    return res;
}

vector<Posit64> gaussian_elimination(vector<vector<Posit64>> &A, vector<Posit64> &b)
{
	int n = A.size();
	
	for(int i = 0; i < n; i++)
	{
		int max_row = i;
		for(int k = i + 1; k < n; k++)
		{
			if(abs(A[k][i]) > abs(A[max_row][i]))
				max_row = k;
		}
		swap(A[i], A[max_row]);
		swap(b[i], b[max_row]);
		
		for(int k = i + 1; k < n; k++)
		{
			Posit64 factor = A[k][i] / A[i][i];
			
			for(int j = i; j < n; j++)
			{
				A[k][j] -= factor * A[i][j];
			}
			b[k] -= factor * b[i];
		}
	}
	
	vector<Posit64> x(n);
	for(int i = n - 1; i >= 0; i--)
	{
		x[i] = b[i];
		for(int j = i + 1; j < n; j++)
		{
			x[i] -= A[i][j] * x[j];
		}
		x[i] /= A[i][i];
	}
	
	return x;
}

Posit64 dot_product(vector<Posit64> a, vector<Posit64> b)
{
	int len = a.size();
	Posit64 sum{0.0};
	
	for(int i = 0; i < len; i++)
	{
		sum += a[i] * b[i];
	}
	
	return sum;
}


void caculate(vector<vector<Posit64>> data, vector<Posit64> x, vector<Posit64> true_val)
{
	int len = data.size();
	Posit64 squared_error_sum{0.0};
	Posit64 abs_error_sum{0.0};
	Posit64 y_sum{0.0};

	for(int i = 0; i < len; i++)
	{	
		Posit64 prediction = dot_product(data[i], x);
		Posit64 error = prediction - true_val[i];
		Posit64 abs_error = Posit_fabs(error);
		squared_error_sum += error * error;
		abs_error_sum += abs_error;
		y_sum += true_val[i];
	}

	Posit64 rmse = Posit_sqrt(squared_error_sum / len);
	Posit64 mae = abs_error_sum / len;
	Posit64 y_mean = y_sum / len;
	Posit64 rmse_percent = (rmse / y_mean) * Posit64{100.0};

	cout << "RMSE: " << setprecision(20) << rmse << " (" << rmse_percent << "% of mean true value)" << endl;
	cout << "MAE: " << setprecision(20) << mae << endl;
}



