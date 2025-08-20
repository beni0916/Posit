#include "Posit/myfdlibm.h"
#include <string>
#include <fstream>
#include <filesystem>
#include <vector>
#include <mpfr.h>
#include <bits/stdc++.h>

namespace fs = std::filesystem;
using namespace std;

void analyzeFloatExponentDistribution(const std::string& analysisName, const double* data, int size, std::ofstream& outputFile);

double calculateMeanRelativeError(const std::vector<double>& errors, const std::vector<double>& mpfrValues);

std::string getFilenameWithoutExtension(const std::string& filename);

std::string getFileExtension(const std::string& filename);

double calculateMedian(const std::vector<double>& vec);

double calculateStandardDeviation(const std::vector<double>& vec);

double RMSE(const std::vector<double> &vec);

std::string difference(std::string& num1, std::string& num2);

std::string formatPositiveFloatString(const std::string& inputStr);

std::string toString(Posit64 input);

string toString(double input);

string toString(mpfr_t input);