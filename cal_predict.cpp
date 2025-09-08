#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

int LEN = 100;
void d1_normalization(vector<double> &);
void d2_normalization(vector<vector<double>> &);
void load_data(string file_name, vector<vector<double>> &, vector<double> &);
void add_exp2_features(vector<vector<double>> &data);
vector<vector<double>> transpose(const vector<vector<double>>);
vector<vector<double>> multiplication(const vector<vector<double>>, const vector<vector<double>>);
vector<double> multiplication(const vector<vector<double>> A, const vector<double> b);
vector<double> gaussian_elimination(vector<vector<double>> &A, vector<double> &b);
double dot_product(vector<double>, vector<double>);
void caculate(vector<vector<double>> data, vector<double> x, vector<double> true_val);

// g++ -std=c++17 -o cal_predict cal_predict.cpp

int main(void)
{
	vector<vector<double>> train_A;
	vector<double> train_b;
	load_data("cal_train.txt", train_A, train_b);
	//cout << "size = " << train_A.size() << endl;

	vector<vector<double>> test_A;
	vector<double> test_b;
	load_data("cal_test.txt", test_A, test_b);

	vector<vector<double>> train_A_t = transpose(train_A);
	vector<vector<double>> test_A_t  = transpose(test_A);

	// for(int i = 0; i < test_A[0].size(); i++)
	// 	cout << test_A[0][i] << endl;
	// cout << endl;

	// normalization
	d2_normalization(train_A_t);
	d2_normalization(test_A_t);
	//cout << "After d2" << endl;
	// for(int i = 0; i < test_A_t.size(); i++)
	// 	cout << test_A_t[i][0] << endl;
	// cout << endl;
	train_A = transpose(train_A_t);
	test_A  = transpose(test_A_t);

	// cout << "after tran" << endl;
	// for(int i = 0; i < test_A[0].size(); i++)
	// 	cout << test_A[0][i] << endl;
	// cout << endl;
	
	// for(int i = 0; i < 10; i++)
	// 	cout << test_b[i] << endl;
	// cout << endl;

	// d1_normalization(train_b);
	// d1_normalization(test_b);

	// for(int i = 0; i < 10; i++)
	// 	cout << test_b[i] << endl;
	// cout << endl;

	vector<vector<double>> ATA = multiplication(train_A_t, train_A);
	vector<double> ATB = multiplication(train_A_t, train_b);

	vector<double> x = gaussian_elimination(ATA, ATB);
	// for(int i = 0; i < x.size(); i++)
	// 	cout << x[i] << endl;
	// cout << endl;

	caculate(test_A, x, test_b);
	//cout << train_A[0][0] << endl;
	return 0;
}

void d1_normalization(vector<double> &src)
{
	int len = src.size();
	double min = src[0], max = src[0], diff = 0, sum = 0.000001;
	
	for(int i = 1; i < len; i++)
	{
		if(src[i] > max)
			max = src[i];
		if(src[i] < min)
			min = src[i];
	}
	
	// for(int i = 0; i < len; i++)
	// {
	// 	src[i] = exp(src[i] - max) ;
	// 	sum += src[i];
	// }

	// for(int i = 0; i < len; i++)
	// {
	// 	src[i] = src[i] / sum;
	// }
	
	diff = max - min;
	if(diff == 0)
		diff = 0.000001;
	for(int i = 0; i < len; i++)
	{
		src[i] = (src[i] - min) / diff;
	}
}

void d2_normalization(vector<vector<double>> &src)
{
	int len = src.size();
	
	for(int i = 0; i < len; i++)
	{
		d1_normalization(src[i]);
	}
}

void load_data(string file_name, vector<vector<double>> &store_arr, vector<double> &target)
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
	vector<vector<double>> train_arr;
	while(getline(train_file, line)) //cnt < LEN &&
	{
		cnt++;
		stringstream ss(line);
		string token;
		vector<double> row;
		int col_count = 0;
		
		while(getline(ss, token, ','))
		{
			if(col_count >= 9)
				break;
				
			double value = 0.0;
			if(!token.empty())
			{
				try
				{
					value = stod(token) + 0.000001;
				}
				catch(...)
				{
					value = 0.0;
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

vector<vector<double>> transpose(const vector<vector<double>> src)
{
	if(src.empty())
		return {};
		
	size_t row = src.size();
	size_t col = src[0].size();
	
	vector<vector<double>> res(col, vector<double>(row));
	
	for(size_t i = 0; i < row; i++)
	{
		for(size_t j = 0; j < col; j++)
		{
			res[j][i] = src[i][j];
		}
	}
	
	return res;
}

vector<vector<double>> multiplication(const vector<vector<double>> a, const vector<vector<double>> b)
{
	int row_a = a.size();
	int col_a = a[0].size();
	int col_b = b[0].size();
	
	vector<vector<double>> res(row_a, vector<double>(col_b, 0.0));
	
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

vector<double> multiplication(const vector<vector<double>> A, const vector<double> b)
{
    int row = A.size();
    int col = A[0].size();
    vector<double> res(row, 0.0);

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            res[i] += A[i][j] * b[j];
        }
    }
    return res;
}

vector<double> gaussian_elimination(vector<vector<double>> &A, vector<double> &b)
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
			double factor = A[k][i] / A[i][i];
			
			for(int j = i; j < n; j++)
			{
				A[k][j] -= factor * A[i][j];
			}
			b[k] -= factor * b[i];
		}
	}
	
	vector<double> x(n);
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

double dot_product(vector<double> a, vector<double> b)
{
	int len = a.size();
	double sum = 0.0;
	
	//cout << endl << endl;	
	for(int i = 0; i < len; i++)
	{
		sum += a[i] * b[i];
		//cout << a[i] << " "  << b[i] << " "  <<  sum << endl;
	}
	
	return sum;
}


void caculate(vector<vector<double>> data, vector<double> x, vector<double> true_val)
{
	int len = data.size();
	double squared_error_sum = 0.0;
	double abs_error_sum = 0.0;
	double y_sum = 0.0;

	for(int i = 0; i < len; i++)
	{	
		double prediction = dot_product(data[i], x);
		double error = prediction - true_val[i];
		double abs_error = abs(error);
		squared_error_sum += error * error;
		abs_error_sum += abs_error;
		y_sum += true_val[i];
		
		// if(i == 99)
		// {
		// 	cout << data[99].size() << endl;
		// 	cout << "size " << x.size() << endl;
		// 	cout << x[0] << endl;
		// 	cout << data[0][0] << endl;
		// 	cout << "prediction " << prediction << endl;
		// 	cout << "error " << error << endl;
		// }
	}

	double rmse = sqrt(squared_error_sum / len);
	double mae = abs_error_sum / len;
	double y_mean = y_sum / len;
	double rmse_percent = (rmse / y_mean) * 100.0;

	//cout << "y_mean" << y_mean << endl;
	cout << "RMSE: " << setprecision(20) << rmse << " (" << rmse_percent << "% of mean true value)" << endl;
	cout << "MAE: " << setprecision(20) << mae << endl;
}

void add_exp2_features(vector<vector<double>> &data)
{
    for (auto &row : data)
    {
        int original_size = row.size();
        for (int i = 0; i < original_size; ++i)
        {
            row.push_back(exp2(abs(row[i]) + 1e-6));
        }
    }
}


