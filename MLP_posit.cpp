#include "myfdlibm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <random>
// g++ -std=c++17 -o MLP_posit MLP_posit.cpp Posit_sqrt.cpp
using namespace std;

// char train_file[30] = "cal_train.txt";
// char test_file[30] = "cal_test.txt";
// char sep_word = '\t';
// int total_col = 9;
// int target_col = 8;

char train_file[30] = "CP_train.txt";
char test_file[30] = "CP_test.txt";
char sep_word = '\t';
int total_col = 15;
int target_col = 14;

vector<vector<Posit64>> relu(vector<vector<Posit64>> src);
vector<vector<Posit64>> relu_deriv(vector<vector<Posit64>> res);
Posit64 mse_loss(vector<Posit64> y_true, vector<vector<Posit64>> y_pred);
void load_data(string file_name, vector<vector<Posit64>> &store_arr, vector<Posit64> &target);
vector<vector<Posit64>> transpose(const vector<vector<Posit64>> src);
void get_mean_std(vector<vector<Posit64>> src, vector<Posit64> &mean, vector<Posit64> &std);
void get_mean_std(vector<Posit64> src, Posit64 &mean, Posit64 &std);
void normalize(vector<vector<Posit64>> &src, vector<Posit64> mean, vector<Posit64> std);
void normalize(vector<Posit64> &src, Posit64 mean, Posit64 std);
vector<vector<Posit64>> matrix_mul(const vector<vector<Posit64>> a, const vector<vector<Posit64>> b);
vector<vector<Posit64>> matrix_add(vector<vector<Posit64>> target, vector<Posit64> src);
vector<vector<Posit64>> matrix_add(vector<vector<Posit64>> target, Posit64 src);
vector<vector<Posit64>> matrix_sub(vector<vector<Posit64>> a, vector<Posit64>b);
vector<Posit64> matrix_sub(vector<Posit64> a, vector<vector<Posit64>> b);
vector<vector<Posit64>> matrix_sub(vector<vector<Posit64>> a, vector<vector<Posit64>>b);
vector<vector<Posit64>> matrix_mul(const vector<vector<Posit64>> a, Posit64 b);
vector<vector<Posit64>> matrix_div(const vector<vector<Posit64>> a, Posit64 b);
vector<vector<Posit64>> matrix_sum(const vector<vector<Posit64>> a, int dir);
vector<vector<Posit64>> hadamard(const vector<vector<Posit64>> &a, const vector<vector<Posit64>> &b);

int main(void)
{
    // load_data
    vector<vector<Posit64>> X_train;
    vector<Posit64> y_train;
    load_data(train_file, X_train, y_train);

    vector<vector<Posit64>> X_test;
    vector<Posit64> y_test;
    load_data(test_file, X_test, y_test);

    cout << "load data success" << endl;
    //-------------------------------------------------------------------------------------------------
    // transpose
    vector<vector<Posit64>> X_train_tran;
    vector<vector<Posit64>> X_test_tran;
    X_train_tran = transpose(X_train);
    X_test_tran = transpose(X_test);
    X_train = X_train_tran;
    X_test = X_test_tran;

    cout << "transpose success" << endl;
    //-------------------------------------------------------------------------------------------------
    // get mean and std
    vector<Posit64> X_train_mean, X_train_std;
    Posit64 y_train_mean, y_train_std;
    get_mean_std(X_train, X_train_mean, X_train_std);
    get_mean_std(y_train, y_train_mean, y_train_std);

    vector<Posit64> X_test_mean, X_test_std;
    get_mean_std(X_test, X_test_mean, X_test_std);
    
    cout << "get_mean_std success" << endl;
    //-------------------------------------------------------------------------------------------------
    // normalize
    normalize(X_train, X_train_mean, X_train_std);
    normalize(X_test, X_test_mean, X_test_std);
    normalize(y_train, y_train_mean, y_train_std);

    cout << "normalize success" << endl;
    //-------------------------------------------------------------------------------------------------
    // parameter initial
    srand(42);
    int n_features = X_train.size();
    int n_hidden = 8;
    int n_outputs = 1;
    Posit64 learning_rate = Posit64{0.01};
    Posit64 best_loss = Posit64{1000000000000000};
    int patience = 3;
    int wait = 0;
    int epochs = 10;

    cout << "parameter ini success" << endl;
    //-------------------------------------------------------------------------------------------------
    // weight initial
    mt19937 gen(42);
    normal_distribution<double> dist(0.0, 1.0);

    // 初始化 W1: n_features x n_hidden，乘以 0.01
    vector<vector<Posit64>> W1(n_features, vector<Posit64>(n_hidden));
    for (int i = 0; i < n_features; i++) {
        for (int j = 0; j < n_hidden; j++) {
            W1[i][j] = Posit64{dist(gen)} * Posit64{0.01};
        }
    }

    // 初始化 b1: 1 x n_hidden，全部 0
    vector<Posit64> b1(n_hidden, Posit64{0.0});

    // 初始化 W2: n_hidden x n_outputs，乘以 0.01
    vector<vector<Posit64>> W2(n_hidden, vector<Posit64>(n_outputs));
    for (int i = 0; i < n_hidden; i++) {
        for (int j = 0; j < n_outputs; j++) {
            W2[i][j] = Posit64{dist(gen)} * Posit64{0.01};
        }
    }

    // 初始化 b2: 1 x n_outputs，全部 0
    vector<Posit64> b2(n_outputs, Posit64{0.0});

    cout << "weight ini success" << endl;
    //-------------------------------------------------------------------------------------------------
    // train
    X_train = transpose(X_train);

    for(int i = 0; i < epochs; ++i)
    {
        // z1 = X_train @ W1 + b1
        vector<vector<Posit64>> z1 = matrix_mul(X_train, W1);
        // cout << "111" << endl;
        z1 = matrix_add(z1, b1);

        // a1 = relu(z1)
        vector<vector<Posit64>> a1 = relu(z1);

        // z2 = a1 @ W2 + b2
        vector<vector<Posit64>> z2 = matrix_mul(a1, W2);
        z2 = matrix_add(z2, b2);
    
        // y_pred = z2
        vector<vector<Posit64>> y_pred = z2;
        
        // cout << "before before loss" << endl;
        // loss = mse_loss(y_train, y_pred)
        Posit64 loss = mse_loss(y_train, y_pred);

        // dL_dy = 2 * (y_pred - y_train) / X_train.shape[0]
        vector<vector<Posit64>> dL_dy = matrix_sub(y_pred, y_train);
        dL_dy = matrix_mul(dL_dy, 2);
        dL_dy = matrix_div(dL_dy, X_train.size());

        // dL_dW2 = a1.T @ dL_dy
        vector<vector<Posit64>> dL_dW2 = matrix_mul(transpose(a1), dL_dy);
        
        // dL_db2 = np.sum(dL_dy, axis=0, keepdims=True)
        vector<vector<Posit64>> dL_db2 = matrix_sum(dL_dy, 0);
        
        // da1 = dL_dy @ W2.T
        vector<vector<Posit64>> da1 = matrix_mul(dL_dy, transpose(W2));

        // dz1 = da1 * relu_deriv(z1)
        vector<vector<Posit64>> dz1 = hadamard(da1, relu_deriv(z1));

        // dL_dW1 = X_train.T @ dz1
        vector<vector<Posit64>> dL_dW1 = matrix_mul(transpose(X_train), dz1);

        // dL_db1 = np.sum(dz1, axis=0, keepdims=True)
        vector<vector<Posit64>> dL_db1 = matrix_sum(dz1, 0);

        // W2 -= learning_rate * dL_dW2
        W2 = matrix_sub(W2, matrix_mul(dL_dW2, learning_rate));

        // b2 -= learning_rate * dL_db2
        b2 = matrix_sub(b2, matrix_mul(dL_db2, learning_rate));

        // W1 -= learning_rate * dL_dW1
        W1 = matrix_sub(W1, matrix_mul(dL_dW1, learning_rate));

        // b1 -= learning_rate * dL_db1
        b1 = matrix_sub(b1, matrix_mul(dL_db1, learning_rate));

        // cout << "before loss" << endl;

        if(loss < best_loss)
        {
            best_loss = loss;
            wait = 0;
        }
        else 
        {
            wait++;
            if (wait >= patience) 
            {
                learning_rate *= Posit64{0.5};   // 減半
                cout << "Reduce lr to " << learning_rate << endl;
                wait = 0;
            }
        }
        
        cout << "i = " << i << endl;
    }
    cout << "train success" << endl;

    // y_pred_test = relu(X_test @ W1 + b1) @ W2 + b2
    vector<vector<Posit64>> y_pred_test = matrix_mul(transpose(X_test), W1);
    y_pred_test = matrix_add(y_pred_test, b1);
    y_pred_test = relu(y_pred_test);
    y_pred_test = matrix_mul(y_pred_test, W2);
    y_pred_test = matrix_add(y_pred_test, b2);

    // y_pred_test = y_pred_test * y_std + y_mean
    y_pred_test = matrix_mul(y_pred_test, y_train_std);
    y_pred_test = matrix_add(y_pred_test, y_train_mean);

    // y_test_orig = y_test
    vector<Posit64> y_test_orig = y_test;

    // test_loss = mse_loss(y_test_orig, y_pred_test)
    Posit64 test_loss = mse_loss(y_test_orig, y_pred_test);
    cout << "Test Loss = " << test_loss << endl;

    // rmse = np.sqrt(test_loss)
    Posit64 rmse = Posit_sqrt(test_loss);

    // rmse_percent = (rmse / np.mean(y_test)) * 100
    Posit64 y_test_mean = 0;
    for(int i = 0; i < y_test.size(); ++i)
        y_test_mean = y_test_mean + y_test[i];
    y_test_mean = y_test_mean / y_test.size();
    Posit64 rmse_percent = (rmse / y_test_mean) * Posit64{100};
    cout << "rmse_percent = " << rmse_percent << "%" << endl;


    cout << "finish" << endl;
    return 0;
}

vector<vector<Posit64>> relu(vector<vector<Posit64>> src)
{
    vector<vector<Posit64>> res = src;  

    for(size_t i = 0; i < src.size(); ++i) 
    {
        for(size_t j = 0; j < src[0].size(); ++j) 
        {
            res[i][j] = (src[i][j] > Posit64{0.0}) ? src[i][j] : Posit64{0.0};
        }
    }
    
    return res;
}

vector<vector<Posit64>> relu_deriv(vector<vector<Posit64>> src)
{
    vector<vector<Posit64>> res = src;  

    for(size_t i = 0; i < src.size(); ++i) 
    {
        for(size_t j = 0; j < src[0].size(); ++j) 
        {
            res[i][j] = (src[i][j] > Posit64{0.0}) ? Posit64{1.0} : Posit64{0.0};
        }
    }
    
    return res;
}

Posit64 mse_loss(vector<Posit64> y_true, vector<vector<Posit64>> y_pred)
{
    Posit64 sum = Posit64{0.0};
    
    for(int i = 0; i < y_true.size(); ++i)
    {
        Posit64 diff = y_true[i] - y_pred[i][0];
        sum = sum + diff * diff;
    }

    return sum / y_true.size();
}

void load_data(string file_name, vector<vector<Posit64>> &store_arr, vector<Posit64> &target)
{
	ifstream train_file(file_name, ios::in);
	if(!train_file)
	{
		perror("open file failed");
        cout << "檔名: " << file_name << endl;
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
				
			Posit64 value = Posit64{0.0};
			if(!token.empty())
			{
				try
				{
					value = stod(token) ; // 0.000001
				}
				catch(...)
				{
					value = Posit64{0.0};
				}
			}
			
			if(col_count == target_col)
			{
                target.push_back(value);
            }
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

void get_mean_std(vector<vector<Posit64>> src, vector<Posit64> &mean, vector<Posit64> &std)
{
    // cout << "start get mean std" << endl;
    for(int i = 0; i < src.size(); ++i)
    {
        // cout << "i = " << i << " start" << endl;
        mean.push_back(0.0);
        for(int j = 0; j < src[i].size(); ++j)
            mean[i] = mean[i] + src[i][j];
        mean[i] = mean[i] / src[i].size();

        std.push_back(0.0);
        for(int j = 0; j < src[i].size(); ++j)
            std[i] = std[i] + (src[i][j] - mean[i]) * (src[i][j] - mean[i]);
        std[i] = std[i] / src[i].size();
        std[i] = Posit_sqrt(std[i]);
        // cout << "i = " << i << " end" << endl;
    }
}

void get_mean_std(vector<Posit64> src, Posit64 &mean, Posit64 &std)
{
    mean = 0;
    for(int i = 0; i < src.size(); ++i)
        mean = mean + src[i];
    mean = mean / src.size();

    std = 0;
    for(int i = 0; i < src.size(); ++i)
        std = std + (src[i] - mean) * (src[i] - mean);
    std = std / src.size();
    std = Posit_sqrt(std);
}

void normalize(vector<vector<Posit64>> &src, vector<Posit64> mean, vector<Posit64> std)
{
    for(int i = 0; i < src.size(); ++i)
    {
        for(int j = 0; j < src[i].size(); ++j)
            src[i][j] = (src[i][j] - mean[i]) / std[i];
    }
}

void normalize(vector<Posit64> &src, Posit64 mean, Posit64 std)
{
    for(int i = 0; i < src.size(); ++i)
        src[i] = (src[i] - mean) / std;
}

vector<vector<Posit64>> matrix_mul(const vector<vector<Posit64>> a, const vector<vector<Posit64>> b)
{
	int row_a = a.size();
	int col_a = a[0].size();
	int col_b = b[0].size();
	
	vector<vector<Posit64>> res(row_a, vector<Posit64>(col_b, 0.0));
	
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

vector<vector<Posit64>> matrix_add(vector<vector<Posit64>> target, vector<Posit64> src)
{
    vector<vector<Posit64>> res = target;

    for(int i = 0; i < target.size(); ++i)
    {
        for(int j = 0; j < target[i].size(); ++j)
        {
            res[i][j] = target[i][j] + src[j];
        }
    }

    return res;
}

vector<vector<Posit64>> matrix_add(vector<vector<Posit64>> target, Posit64 src)
{
    vector<vector<Posit64>> res = target;

    for(int i = 0; i < target.size(); ++i)
    {
        for(int j = 0; j < target[i].size(); ++j)
        {
            res[i][j] = target[i][j] + src;
        }
    }

    return res;
}

vector<vector<Posit64>> matrix_sub(vector<vector<Posit64>> a, vector<Posit64>b)
{
    vector<vector<Posit64>> res = a;

    for(int i = 0; i < a.size(); ++i)
        res[i][0] = a[i][0] - b[i];
    return res;
}

vector<Posit64> matrix_sub(vector<Posit64> a, vector<vector<Posit64>> b)
{
    vector<Posit64> res = a;

    for(int i = 0; i < a.size(); ++i)
        res[i] = a[i] - b[0][i];
    
    return res;
}

vector<vector<Posit64>> matrix_sub(vector<vector<Posit64>> a, vector<vector<Posit64>>b)
{
    vector<vector<Posit64>> res = a;

    for(int i = 0; i < a.size(); ++i)
        for(int j = 0; j < a[i].size(); ++j)
            res[i][j] = a[i][j] - b[i][j];
    return res;
}

vector<vector<Posit64>> matrix_mul(const vector<vector<Posit64>> a, Posit64 b)
{
    vector<vector<Posit64>> res = a;
    for(int i = 0; i < a.size(); ++i)
        for(int j = 0; j < a[i].size(); ++j)
            res[i][j] = a[i][j] * b;
    return res;
}

vector<vector<Posit64>> matrix_div(const vector<vector<Posit64>> a, Posit64 b)
{
    vector<vector<Posit64>> res = a;
    for(int i = 0; i < a.size(); ++i)
        for(int j = 0; j < a[i].size(); ++j)
            res[i][j] = a[i][j] / b;
    return res;
}

vector<vector<Posit64>> matrix_sum(const vector<vector<Posit64>> a, int dir)
{
    vector<vector<Posit64>> res(1, vector<Posit64>(a[0].size(), 0.0));

    if(dir == 0)
    {
        for(int i = 0; i < a[0].size(); ++i)
        {
            for(int j = 0; j < a.size(); ++j)
            {
                res[0][i] = res[0][i] + a[j][i];
            }
        }
    }

    return res;
}

vector<vector<Posit64>> hadamard(const vector<vector<Posit64>> &a, const vector<vector<Posit64>> &b) 
{
    if (a.size() == 0 || a[0].size() == 0) return {};
    int rows = a.size();
    int cols = a[0].size();
    // 可加檢查 b shape 是否一樣
    if ((int)b.size() != rows || (int)b[0].size() != cols) {
        cerr << "hadamard shape mismatch\n";
        exit(1);
    }
    vector<vector<Posit64>> res(rows, vector<Posit64>(cols, 0.0));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            res[i][j] = a[i][j] * b[i][j];
        }
    }
    return res;
}
