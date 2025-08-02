#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"
#include "Posit/myfdlibm.h"
#include <gmp.h>
#include <mpfr.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <filesystem>
#include <string>
#include <iomanip>
#include <bits/stdc++.h>
#include <fstream>

namespace fs = std::filesystem;
using namespace std;

std::string getFilenameWithoutExtension(const std::string& filename);

std::string getFileExtension(const std::string& filename);

double RMSE(const std::vector<double> &vec);

std::string difference(std::string& num1, std::string& num2);

std::string formatPositiveFloatString(const std::string& inputStr);

std::string toString(Posit64 input);

string toString(double input);

string toString(mpfr_t input);