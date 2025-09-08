#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include "myfdlibm.h"
// g++ -std=c++17 -o test_posit_v1 test_posit_v1.cpp Posit_fabs.cpp Posit_sqrt.cpp Posit_exp2.cpp Posit_floor.cpp Posit_log2.cpp
using namespace std;
int LEN = 2500;
#define arr_len 14
char train_arr[arr_len][30] = {"linear10", "linear11", "linear12", "linear13", "linear14", "linear15", "linear16", "linear17", "linear18", "linear19", "linear20", "linear21", "linear23", "linear24"};
char sep_arr[arr_len + 1] = "\t\t\t\t\t\t\t\t\t\t\t\t\t\t";
int total_arr[arr_len] = {11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11};
int target_arr[arr_len] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};

char train_file[30];
char test_file[30];
char sep_word;
int total_col;
int target_col;
void load_data(string file_name, vector<vector<Posit64>> &, vector<Posit64> &);
vector<vector<Posit64>> transpose(const vector<vector<Posit64>>);
void d1_normalization(vector<Posit64> &, Posit64 &, Posit64 &);
void d2_normalization(vector<vector<Posit64>> &);
void add_exp2_features(vector<vector<Posit64>> &data);
vector<vector<Posit64>> multiplication(const vector<vector<Posit64>>, const vector<vector<Posit64>>);
vector<Posit64> multiplication(const vector<vector<Posit64>> A, const vector<Posit64> b);
vector<Posit64> gaussian_elimination(vector<vector<Posit64>> &A, vector<Posit64> &b);
Posit64 dot_product(vector<Posit64>, vector<Posit64>);
void caculate(vector<vector<Posit64>> data, vector<Posit64> x, vector<Posit64> true_val, Posit64&, Posit64&);
void t_frexp(vector<Posit64> &src, Posit64 &max, Posit64 &min);
void move_data(vector<Posit64> &src);

int main(void)
{
    for(int i = 0; i < arr_len; ++i)
    {
        strcpy(train_file, train_arr[i]);
        strcat(train_file, "_train.txt");
        strcpy(test_file, train_arr[i]);
        strcat(test_file, "_test.txt");
        sep_word = sep_arr[i];
        total_col = total_arr[i];
        target_col = target_arr[i];
	// Posit64 test0{0.9}, test1{1.1};
	// __uint32_t test2 = test0.bits(), test3 = test1.bits();
	// cout << "0.9 = " << bitset<32>(test2) << endl;
	// cout << "1.1 = " << bitset<32>(test3) << endl;

        vector<vector<Posit64>> train_A;
        vector<Posit64> train_b;
        load_data(train_file, train_A, train_b);
        
        vector<vector<Posit64>> test_A;
        vector<Posit64> test_b;
        // LEN = 1000;
        load_data(test_file, test_A, test_b);
        
        //cout << "load finish" << endl;

        vector<vector<Posit64>> train_A_t = transpose(train_A);
        vector<vector<Posit64>> test_A_t  = transpose(test_A);

        //cout << "tran finish" << endl;
        move_data(train_b);
        move_data(test_b);

        for(int i = 0; i < train_A_t.size(); i++)
            move_data(train_A_t[i]);
        for(int i = 0; i < test_A_t.size(); i++)
            move_data(test_A_t[i]);

        // normalization
        d2_normalization(train_A_t);
        d2_normalization(test_A_t);
        //cout << "d2 nor finish" << endl;

        train_A = transpose(train_A_t);
        test_A  = transpose(test_A_t);
        //cout << "tran finish" << endl;
        // -------------------

        Posit64 max{0.0}, min{0.0}, test_max, test_min;
        // d1_normalization(train_b, max, min);
        t_frexp(train_b, max, min);
        //d1_normalization(test_b, test_max, test_min);
        //cout << "d1 nor finish" << endl;
        
        vector<vector<Posit64>> ATA = multiplication(train_A_t, train_A);
        vector<Posit64> ATB = multiplication(train_A_t, train_b);
        //cout << "mul finish" << endl;

        vector<Posit64> x = gaussian_elimination(ATA, ATB);
        //cout << "gauss finish" << endl;
        caculate(test_A, x, test_b, max, min);
    }
	return 0;
}

void d1_normalization(vector<Posit64> &src, Posit64 &nor_max, Posit64 &nor_min)
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
	
	nor_max = max;
	nor_min = min;

	diff = max - min;
	if(diff == Posit64{0})
		diff = Posit64{0.000001};
	for(int i = 0; i < len; i++)
	{
		src[i] = Posit64{0.9} + (src[i] - min) / diff * Posit64{0.2};
	}
}

void d2_normalization(vector<vector<Posit64>> &src)
{
	int len = src.size();
	Posit64 max{0}, min{0};
	
	for(int i = 0; i < len; i++)
	{
		// d1_normalization(src[i], max, min);
        t_frexp(src[i], max, min);
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
	
	int cnt = 0;
	string line;
	getline(train_file, line);
	vector<vector<Posit64>> train_arr;
	while(getline(train_file, line)) //cnt < LEN && 
	{
		cnt++;
		stringstream ss(line);
		string token;
		vector<Posit64> row;
		int col_count = 0;
		
		while(getline(ss, token, sep_word))
		{
			if(col_count >= total_col)
				break;
				
			Posit64 value{0.0};
			if(!token.empty())
			{
				try
				{
					value = Posit64{stod(token)}; //  + Posit64{0.000001} 
				}
				catch(...)
				{
					value = Posit64{0.0};
				}
			}
			
			if(col_count == target_col)
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
        // cout << "i = " << i << endl;
		for(int j = 0; j < col_b; j++)
		{
			for(int k = 0; k < col_a; k++)
			{
                //cout << "i = " << i << " j = " << j << " k = " << k << endl;
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
        // cout << "i = " << i << endl;
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


void caculate(vector<vector<Posit64>> data, vector<Posit64> x, vector<Posit64> true_val, Posit64 &max, Posit64 &min)
{
	int len = data.size();
	Posit64 squared_error_sum{0.0};
	Posit64 abs_error_sum{0.0};
	Posit64 y_sum{0.0};
	Posit64 diff = max - min;
	if(diff == 0)
        diff = Posit64{1e-6};

	for(int i = 0; i < len; i++)
	{	
		Posit64 prediction = dot_product(data[i], x);
		// prediction = (prediction - Posit64{0.9}) / Posit64{0.2} * diff + min;
        prediction = prediction * Posit_exp2((max - min) / Posit64{2});

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

	
	cout << train_file << endl;
	cout << "RMSE: " << setprecision(20) << rmse << " (" << rmse_percent << "% of mean true value)" << endl;
	cout << "MAE: " << setprecision(20) << mae << endl << endl;
}

void add_exp2_features(vector<vector<Posit64>> &data)
{
    for (auto &row : data)
    {
        int original_size = row.size();
        for (int i = 0; i < original_size; ++i)
        {
            row.push_back(Posit_exp2(Posit_fabs(row[i]) + 1e-6));
        }
    }
}

void t_frexp(vector<Posit64> &src, Posit64 &max, Posit64 &min)
{
    Posit64 t_max = Posit_log2(src[0]), t_min = Posit_log2(src[0]);

	t_max = Posit_floor(t_max);
	t_min = Posit_floor(t_min);
	vector<Posit64> save_arr;
	save_arr.push_back(t_max);

    for(int i = 1; i < src.size(); i++)
    {
		Posit64 tmp = Posit_log2(src[i]);
		tmp = Posit_floor(tmp);
		save_arr.push_back(tmp);

        if(tmp > t_max)
        {
            t_max = tmp;
        }
        if(tmp < t_min)
            t_min = tmp;
    }

    max = t_max;
    min = t_min;

    for(int i = 0; i < src.size(); i++)
    {
        Posit64 e = save_arr[i];
		Posit64 tmp = (Posit_exp2(e));
		Posit64 v{0.0};
		if(tmp != Posit64{0})
        	v = src[i] / tmp;

	
		if(v > Posit64{1.0} || v < Posit64{-1.0})
		{
			v = v / Posit64{2.0};
			e = e + 1;
		}
        Posit64 e_0 = (t_max - t_min) / 2;

		// if(v == NAR)
		// {
		// 	cout << i << endl;
		// 	cout << src[i] << endl;
		// 	cout << e << endl << endl;
		// }

        src[i] = v * Posit_exp2(e - e_0);
	}
}

void move_data(vector<Posit64> &src)
{
    Posit64 min = src[0];
    for(int i = 1; i < src.size(); ++i)
    {   
        if(src[i] < min)
            min = src[i];
    }
    min = Posit_fabs(min);

    for(int i = 0; i < src.size(); ++i)
        src[i] = src[i] + min + Posit64{1};
}
