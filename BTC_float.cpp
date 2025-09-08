#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <iomanip>
// g++ -std=c++17 -o BTC_float BTC_float.cpp
using namespace std;
int LEN = 2500;
#define arr_len 31

char train_file[30] = "BTC_train.txt";
char test_file[30] = "BTC_test.txt";
char sep_word = '\t';
int total_col = 6;
int target_col = 4;
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
		src[i] = 0.9 + (src[i] - min) / diff * 0.2;
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
		prediction = (prediction - 0.9) / 0.2 * diff + min;

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



