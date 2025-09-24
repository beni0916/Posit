#include "myfdlibm.h"
#include <string>
#include <bits/stdc++.h>
#include <fstream>
#include <filesystem>
#include <vector>
#include <gmp.h>
#include <mpfr.h>
#include <bits/stdc++.h>

namespace fs = std::filesystem;
using namespace std;

double RMSE(const std::vector<double> &vec);

std::string difference(std::string num1, std::string num2);

std::string add(std::string num1, std::string num2);

std::string formatPositiveFloatString(const std::string& inputStr, bool isNegative);

std::string toString(Posit64 input);

string toString(double input);

string toString(mpfr_t input);