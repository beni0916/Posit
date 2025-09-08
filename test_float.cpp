#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <iomanip>
// g++ -std=c++17 -o test_float test_float.cpp
using namespace std;
int LEN = 2500;
// #define arr_len 31
// char train_arr[arr_len][30] = {"AME", "cal", "auto", "concrete", "AHD", "HPD", "HPI", "insurance", "linear", "linear1", "linear2", "linear3", "linear4", "linear5", "linear6", "linear7", "linear8", "linear9", "hetero", "boundary", "collinear", "collinear_extreme", "collinear_high_p", "collinear_high", "collinear_low_p", "collinear_low", "collinear_mid_p", "collinear_mid", "CP", "btsp", "melb"};
// char sep_arr[arr_len + 1] = "\t,\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t";
// int total_arr[arr_len] = {37, 9, 8, 9, 36, 13, 22, 4, 10, 10, 20, 20, 20, 30, 30, 32, 33, 23, 28, 13, 11, 11, 7, 11, 7, 11, 7, 11, 15, 6, 13};
// int target_arr[arr_len] = {36, 8, 0, 8, 35, 0, 21, 3, 9, 9, 19, 19, 19, 29, 29, 31, 32, 22, 27, 12, 10, 10, 6, 10, 6, 10, 6, 10, 14, 5, 1};

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
void d1_normalization(vector<float> &, float &, float&);
void d2_normalization(vector<vector<float>> &);
void load_data(string file_name, vector<vector<float>> &, vector<float> &);
vector<vector<float>> transpose(const vector<vector<float>>);
vector<vector<float>> multiplication(const vector<vector<float>>, const vector<vector<float>>);
vector<float> multiplication(const vector<vector<float>> A, const vector<float> b);
vector<float> gaussian_elimination(vector<vector<float>> &A, vector<float> &b);
float dot_product(vector<float>, vector<float>);
void caculate(vector<vector<float>> data, vector<float> x, vector<float> true_val, float &, float &);

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

        vector<vector<float>> train_A;
        vector<float> train_b;
        load_data(train_file, train_A, train_b);
        
        vector<vector<float>> test_A;
        vector<float> test_b;
        // LEN = 1000;
        load_data(test_file, test_A, test_b);
        
        vector<vector<float>> train_A_t = transpose(train_A);
        vector<vector<float>> test_A_t  = transpose(test_A);
        
        // normalization
        d2_normalization(train_A_t);
        d2_normalization(test_A_t);
        train_A = transpose(train_A_t);
        test_A  = transpose(test_A_t);

        float max = 0, min = 0, test_max, test_min;
        d1_normalization(train_b, max, min);
        //d1_normalization(test_b, test_max, test_min);
        
        vector<vector<float>> ATA = multiplication(train_A_t, train_A);
        vector<float> ATB = multiplication(train_A_t, train_b);
        
        vector<float> x = gaussian_elimination(ATA, ATB);
        caculate(test_A, x, test_b, max, min);
    }
	return 0;
}

void d1_normalization(vector<float> &src, float &nor_max, float &nor_min)
{
	int len = src.size();
	float min = src[0], max = src[0], diff = 0.0f;
	
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
	if(diff == 0)
		diff = 0.000001f;
	for(int i = 0; i < len; i++)
	{
		// src[i] = 0.9 + (src[i] - min) / diff * 0.2;
		src[i] = 0.00390625 + (src[i] - min) / diff * 127.99609375;
	}
}

void d2_normalization(vector<vector<float>> &src)
{
	int len = src.size();
	
	float max = 0, min = 0;

	for(int i = 0; i < len; i++)
	{
		d1_normalization(src[i], max, min);
	}
}

void load_data(string file_name, vector<vector<float>> &store_arr, vector<float> &target)
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
	vector<vector<float>> train_arr;
	while(getline(train_file, line)) //cnt < LEN && 
	{
		cnt++;
		stringstream ss(line);
		string token;
		vector<float> row;
		int col_count = 0;
		
		while(getline(ss, token, sep_word))
		{
			if(col_count >= total_col)
				break;
				
			float value = 0.0f;
			if(!token.empty())
			{
				try
				{
					value = stod(token); // + 0.000001f
				}
				catch(...)
				{
					value = 0.0f;
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

vector<vector<float>> transpose(const vector<vector<float>> src)
{
	if(src.empty())
		return {};
		
	size_t row = src.size();
	size_t col = src[0].size();
	
	vector<vector<float>> res(col, vector<float>(row));
	
	for(size_t i = 0; i < row; i++)
	{
		for(size_t j = 0; j < col; j++)
		{
			res[j][i] = src[i][j];
		}
	}
	
	return res;
}

vector<vector<float>> multiplication(const vector<vector<float>> a, const vector<vector<float>> b)
{
	int row_a = a.size();
	int col_a = a[0].size();
	int col_b = b[0].size();
	
	vector<vector<float>> res(row_a, vector<float>(col_b, 0.0f));
	
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

vector<float> multiplication(const vector<vector<float>> A, const vector<float> b)
{
    int row = A.size();
    int col = A[0].size();
    vector<float> res(row, 0.0);

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            res[i] += A[i][j] * b[j];
        }
    }
    return res;
}

vector<float> gaussian_elimination(vector<vector<float>> &A, vector<float> &b)
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
			float factor = A[k][i] / A[i][i];
			
			for(int j = i; j < n; j++)
			{
				A[k][j] -= factor * A[i][j];
			}
			b[k] -= factor * b[i];
		}
	}
	
	vector<float> x(n);
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

float dot_product(vector<float> a, vector<float> b)
{
	int len = a.size();
	float sum = 0.0f;
	
	for(int i = 0; i < len; i++)
	{
		sum += a[i] * b[i];
	}
	
	return sum;
}


void caculate(vector<vector<float>> data, vector<float> x, vector<float> true_val, float &max, float &min)
{
	int len = data.size();
	float squared_error_sum = 0.0f;
	float abs_error_sum = 0.0f;
	float y_sum = 0.0f;
	float diff = max - min;
	if(diff == 0)
		diff = 1e-6;


	for(int i = 0; i < len; i++)
	{	
		float prediction = dot_product(data[i], x);
		// prediction = (prediction - 0.9) / 0.2 * diff + min;
		prediction = (prediction - 0.00390625) / 127.99609375 * diff + min;

		float error = prediction - true_val[i];
		float abs_error = abs(error);
		squared_error_sum += error * error;
		abs_error_sum += abs_error;
		y_sum += true_val[i];
	}

	float rmse = sqrt(squared_error_sum / len);
	float mae = abs_error_sum / len;
	float y_mean = y_sum / len;
	float rmse_percent = (rmse / y_mean) * 100.0;
    cout << train_file << ":" << endl;
	cout << "RMSE: " << setprecision(20) << rmse << " (" << rmse_percent << "% of mean true value)" << endl;
	cout << "MAE: " << setprecision(20) << mae << endl << endl;
}



