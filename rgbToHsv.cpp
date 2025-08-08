#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"
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
#include "testTool.h"

namespace fs = std::filesystem;
using namespace std;
using Posit32 = sw::universal::posit<32, 2>;
using Posit16_1 = sw::universal::posit<16, 1>;
using Posit16_2 = sw::universal::posit<16, 2>;

// RGB to HSV 轉換函數 (浮點數版本)
void rgbToHsvFloatStb(int width, int height, const unsigned char* rgbData, float* hsvData) {
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int index = (i * width + j) * 3;
            float r = static_cast<float>(rgbData[index + 0]) / 255.0f;
            float g = static_cast<float>(rgbData[index + 1]) / 255.0f;
            float b = static_cast<float>(rgbData[index + 2]) / 255.0f;

            float h, s, v;
            float maxValue = std::max(r, std::max(g, b));
            float minValue = std::min(r, std::min(g, b));
            float delta = maxValue - minValue;

            v = maxValue;

            if (delta == 0) {
                h = 0;
                s = 0;
            } else {
                s = delta / maxValue;
                if (r == maxValue) {
                    h = (g - b) / delta;
                } else if (g == maxValue) {
                    h = 2 + (b - r) / delta;
                } else {
                    h = 4 + (r - g) / delta;
                }
                h *= 60;
                if (h < 0) {
                    h += 360;
                }
            }

            int hsvIndex = (i * width + j) * 3;
            hsvData[hsvIndex + 0] = h;
            hsvData[hsvIndex + 1] = s;
            hsvData[hsvIndex + 2] = v;
        }
    }
}

// RGB to HSV 轉換函數 (Posit64 版本)
void rgbToHsvPosit64Stb(int width, int height, const unsigned char* rgbData, Posit64* hsvData) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            int index = (i * width + j) * 3;
            Posit64 r = Posit64(rgbData[index + 0]) / Posit64(255.0);
            Posit64 g = Posit64(rgbData[index + 1]) / Posit64(255.0);
            Posit64 b = Posit64(rgbData[index + 2]) / Posit64(255.0);

            Posit64 h, s, v;
            Posit64 maxValue = std::max(r, std::max(g, b));
            Posit64 minValue = std::min(r, std::min(g, b));
            Posit64 delta = maxValue - minValue;

            v = maxValue;

            if (delta == Posit64(0)) {
                h = Posit64(0);
                s = Posit64(0);
            } else {
                s = delta / maxValue;
                if (r == maxValue) {
                    h = (g - b) / delta;
                } else if (g == maxValue) {
                    h = Posit64(2) + (b - r) / delta;
                } else {
                    h = Posit64(4) + (r - g) / delta;
                }
                h *= Posit64(60);
                if (h < Posit64(0)) {
                    h += Posit64(360);
                }
            }

            int hsvIndex = (i * width + j) * 3;
            hsvData[hsvIndex + 0] = h;
            hsvData[hsvIndex + 1] = s;
            hsvData[hsvIndex + 2] = v;
        }
    }
}

// RGB to HSV 轉換函數 (Posit32 版本)
void rgbToHsvPosit32Stb(int width, int height, const unsigned char* rgbData, Posit32* hsvData) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            int index = (i * width + j) * 3;
            Posit32 r = Posit32(rgbData[index + 0]) / Posit32(255.0);
            Posit32 g = Posit32(rgbData[index + 1]) / Posit32(255.0);
            Posit32 b = Posit32(rgbData[index + 2]) / Posit32(255.0);

            Posit32 h, s, v;
            Posit32 maxValue = std::max(r, std::max(g, b));
            Posit32 minValue = std::min(r, std::min(g, b));
            Posit32 delta = maxValue - minValue;

            v = maxValue;

            if (delta == Posit32(0)) {
                h = Posit32(0);
                s = Posit32(0);
            } else {
                s = delta / maxValue;
                if (r == maxValue) {
                    h = (g - b) / delta;
                } else if (g == maxValue) {
                    h = Posit32(2) + (b - r) / delta;
                } else {
                    h = Posit32(4) + (r - g) / delta;
                }
                h *= Posit32(60);
                if (h < Posit32(0)) {
                    h += Posit32(360);
                }
            }

            int hsvIndex = (i * width + j) * 3;
            hsvData[hsvIndex + 0] = h;
            hsvData[hsvIndex + 1] = s;
            hsvData[hsvIndex + 2] = v;
        }
    }
}

// RGB to HSV 轉換函數 (Posit16_1 版本)
void rgbToHsvPosit16_1Stb(int width, int height, const unsigned char* rgbData, Posit16_1* hsvData) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            int index = (i * width + j) * 3;
            Posit16_1 r = Posit16_1(rgbData[index + 0]) / Posit16_1(255.0);
            Posit16_1 g = Posit16_1(rgbData[index + 1]) / Posit16_1(255.0);
            Posit16_1 b = Posit16_1(rgbData[index + 2]) / Posit16_1(255.0);

            Posit16_1 h, s, v;
            Posit16_1 maxValue = std::max(r, std::max(g, b));
            Posit16_1 minValue = std::min(r, std::min(g, b));
            Posit16_1 delta = maxValue - minValue;

            v = maxValue;

            if (delta == Posit16_1(0)) {
                h = Posit16_1(0);
                s = Posit16_1(0);
            } else {
                s = delta / maxValue;
                if (r == maxValue) {
                    h = (g - b) / delta;
                } else if (g == maxValue) {
                    h = Posit16_1(2) + (b - r) / delta;
                } else {
                    h = Posit16_1(4) + (r - g) / delta;
                }
                h *= Posit16_1(60);
                if (h < Posit16_1(0)) {
                    h += Posit16_1(360);
                }
            }

            int hsvIndex = (i * width + j) * 3;
            hsvData[hsvIndex + 0] = h;
            hsvData[hsvIndex + 1] = s;
            hsvData[hsvIndex + 2] = v;
        }
    }
}

// RGB to HSV 轉換函數 (Posit16_2 版本)
void rgbToHsvPosit16_2Stb(int width, int height, const unsigned char* rgbData, Posit16_2* hsvData) {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            int index = (i * width + j) * 3;
            Posit16_2 r = Posit16_2(rgbData[index + 0]) / Posit16_2(255.0);
            Posit16_2 g = Posit16_2(rgbData[index + 1]) / Posit16_2(255.0);
            Posit16_2 b = Posit16_2(rgbData[index + 2]) / Posit16_2(255.0);

            Posit16_2 h, s, v;
            Posit16_2 maxValue = std::max(r, std::max(g, b));
            Posit16_2 minValue = std::min(r, std::min(g, b));
            Posit16_2 delta = maxValue - minValue;

            v = maxValue;

            if (delta == Posit16_2(0)) {
                h = Posit16_2(0);
                s = Posit16_2(0);
            } else {
                s = delta / maxValue;
                if (r == maxValue) {
                    h = (g - b) / delta;
                } else if (g == maxValue) {
                    h = Posit16_2(2) + (b - r) / delta;
                } else {
                    h = Posit16_2(4) + (r - g) / delta;
                }
                h *= Posit16_2(60);
                if (h < Posit16_2(0)) {
                    h += Posit16_2(360);
                }
            }

            int hsvIndex = (i * width + j) * 3;
            hsvData[hsvIndex + 0] = h;
            hsvData[hsvIndex + 1] = s;
            hsvData[hsvIndex + 2] = v;
        }
    }
}


// MPFR 版本的 RGB to HSV 轉換函數 (高精度)
void rgbToHsvMpfrStb(int width, int height, const unsigned char* rgbData, mpfr_t* hsvData) {
    mpfr_prec_t prec = 256; // 设置 MPFR 精度 
    mpfr_t r, g, b, maxVal, minVal, delta, temp1, temp2, temp60;

    // 初始化 MPFR 变量
    mpfr_init2(r, prec);
    mpfr_init2(g, prec);
    mpfr_init2(b, prec);
    mpfr_init2(maxVal, prec);
    mpfr_init2(minVal, prec);
    mpfr_init2(delta, prec);
    mpfr_init2(temp1, prec);
    mpfr_init2(temp2, prec);
    mpfr_init2(temp60, prec);
    mpfr_set_d(temp60, 60.0, MPFR_RNDN);

    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int index = (i * width + j) * 3;

            // 将 RGB 值转换为 MPFR 浮点数
            mpfr_set_d(r, static_cast<double>(rgbData[index + 0]) / 255.0, MPFR_RNDN);
            mpfr_set_d(g, static_cast<double>(rgbData[index + 1]) / 255.0, MPFR_RNDN);
            mpfr_set_d(b, static_cast<double>(rgbData[index + 2]) / 255.0, MPFR_RNDN);

            // 找到最大值、最小值和色差
            mpfr_max(maxVal, r, g, MPFR_RNDN);
            mpfr_max(maxVal, maxVal, b, MPFR_RNDN);
            mpfr_min(minVal, r, g, MPFR_RNDN);
            mpfr_min(minVal, minVal, b, MPFR_RNDN);
            mpfr_sub(delta, maxVal, minVal, MPFR_RNDN);

            // 计算 V
            mpfr_t& v = hsvData[(i * width + j) * 3 + 2];
            mpfr_init2(v, prec);
            mpfr_set(v, maxVal, MPFR_RNDN);

            // 如果 delta 是 0，则 S 和 H 都是 0
            if (mpfr_cmp_d(delta, 0.0) == 0) {
                mpfr_t& h = hsvData[(i * width + j) * 3 + 0];
                mpfr_init2(h, prec);
                mpfr_set_d(h, 0.0, MPFR_RNDN);

                mpfr_t& s = hsvData[(i * width + j) * 3 + 1];
                mpfr_init2(s, prec);
                mpfr_set_d(s, 0.0, MPFR_RNDN);

            } else {
                // 计算 S
                mpfr_t& s = hsvData[(i * width + j) * 3 + 1];
                mpfr_init2(s, prec);
                mpfr_div(s, delta, maxVal, MPFR_RNDN);

                // 计算 H
                mpfr_t& h = hsvData[(i * width + j) * 3 + 0];
                mpfr_init2(h, prec);
                if (mpfr_equal_p(maxVal, r)) {  // maxVal == r
                    mpfr_sub(temp1, g, b, MPFR_RNDN);
                    mpfr_div(h, temp1, delta, MPFR_RNDN);
                } else if (mpfr_equal_p(maxVal, g)) { // maxVal == g
                    mpfr_set_d(temp1, 2.0, MPFR_RNDN);
                    mpfr_sub(temp2, b, r, MPFR_RNDN);
                    mpfr_div(temp2, temp2, delta, MPFR_RNDN);
                    mpfr_add(h, temp1, temp2, MPFR_RNDN);
                } else { // maxVal == b
                    mpfr_set_d(temp1, 4.0, MPFR_RNDN);
                    mpfr_sub(temp2, r, g, MPFR_RNDN);
                    mpfr_div(temp2, temp2, delta, MPFR_RNDN);
                    mpfr_add(h, temp1, temp2, MPFR_RNDN);
                }
                mpfr_mul(h, h, temp60, MPFR_RNDN);
                if (mpfr_cmp_d(h, 0.0) < 0) {
                    mpfr_add_d(h, h, 360.0, MPFR_RNDN);
                }
            }
        }
    }

    mpfr_clear(r);
    mpfr_clear(g);
    mpfr_clear(b);
    mpfr_clear(maxVal);
    mpfr_clear(minVal);
    mpfr_clear(delta);
    mpfr_clear(temp1);
    mpfr_clear(temp2);
    mpfr_clear(temp60);
}

// HSV to RGB 轉換函數 (Float 版本)
void hsvToRgbFloatStb(int width, int height, const float* hsvData, unsigned char* rgbData) {
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int hsvIndex = (i * width + j) * 3;
            float h = hsvData[hsvIndex + 0];
            float s = hsvData[hsvIndex + 1];
            float v = hsvData[hsvIndex + 2];

            float r = 0, g = 0, b = 0;
            if (s == 0) {
                r = g = b = v; // achromatic
            } else {
                float hh = h / 60.0f;
                int ii = static_cast<int>(std::floor(hh));
                float ff = hh - static_cast<float>(ii);
                float p = v * (1.0f - s);
                float q = v * (1.0f - s * ff);
                float t = v * (1.0f - s * (1.0f - ff));

                switch (ii) {
                    case 0:
                        r = v;
                        g = t;
                        b = p;
                        break;
                    case 1:
                        r = q;
                        g = v;
                        b = p;
                        break;
                    case 2:
                        r = p;
                        g = v;
                        b = t;
                        break;
                    case 3:
                        r = p;
                        g = q;
                        b = v;
                        break;
                    case 4:
                        r = t;
                        g = p;
                        b = v;
                        break;
                    default: // case 5
                        r = v;
                        g = p;
                        b = q;
                        break;
                }
            }

            rgbData[(i * width + j) * 3 + 0] = static_cast<unsigned char>(std::min(255.0f, std::max(0.0f, r * 255.0f)));
            rgbData[(i * width + j) * 3 + 1] = static_cast<unsigned char>(std::min(255.0f, std::max(0.0f, g * 255.0f)));
            rgbData[(i * width + j) * 3 + 2] = static_cast<unsigned char>(std::min(255.0f, std::max(0.0f, b * 255.0f)));
        }
    }
}


// HSV to RGB 轉換函數 (Posit64 版本)
void hsvToRgbPosit64Stb(int width, int height, const Posit64* hsvData, unsigned char* rgbData) {
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int hsvIndex = (i * width + j) * 3;
            Posit64 h = hsvData[hsvIndex + 0];
            Posit64 s = hsvData[hsvIndex + 1];
            Posit64 v = hsvData[hsvIndex + 2];

            Posit64 r = 0, g = 0, b = 0;
            if (s == Posit64(0)) {
                r = g = b = v; // achromatic
            } else {
                Posit64 hh = h / Posit64(60);
                int ii = static_cast<int>(Posit_floor(hh));
                Posit64 ff = hh - Posit64(ii);
                Posit64 p = v * (Posit64(1) - s);
                Posit64 q = v * (Posit64(1) - s * ff);
                Posit64 t = v * (Posit64(1) - s * (Posit64(1) - ff));

                switch (ii) {
                    case 0:
                        r = v;
                        g = t;
                        b = p;
                        break;
                    case 1:
                        r = q;
                        g = v;
                        b = p;
                        break;
                    case 2:
                        r = p;
                        g = v;
                        b = t;
                        break;
                    case 3:
                        r = p;
                        g = q;
                        b = v;
                        break;
                    case 4:
                        r = t;
                        g = p;
                        b = v;
                        break;
                    default:
                        r = v;
                        g = p;
                        b = q;
                        break;
                }
            }

            rgbData[(i * width + j) * 3 + 0] = static_cast<unsigned char>((int)Posit_floor(r * Posit64(255.0)));
            rgbData[(i * width + j) * 3 + 1] = static_cast<unsigned char>((int)Posit_floor(g * Posit64(255.0)));
            rgbData[(i * width + j) * 3 + 2] = static_cast<unsigned char>((int)Posit_floor(b * Posit64(255.0)));
        }
    }
}

// HSV to RGB 轉換函數 (Posit32 版本)
void hsvToRgbPosit32Stb(int width, int height, const Posit32* hsvData, unsigned char* rgbData) {
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int hsvIndex = (i * width + j) * 3;
            Posit32 h = hsvData[hsvIndex + 0];
            Posit32 s = hsvData[hsvIndex + 1];
            Posit32 v = hsvData[hsvIndex + 2];

            Posit32 r = 0, g = 0, b = 0;
            if (s == Posit32(0)) {
                r = g = b = v; // achromatic
            } else {
                Posit32 hh = h / Posit32(60);
                int ii = static_cast<int>(Posit_floor(hh));
                Posit32 ff = hh - Posit32(ii);
                Posit32 p = v * (Posit32(1) - s);
                Posit32 q = v * (Posit32(1) - s * ff);
                Posit32 t = v * (Posit32(1) - s * (Posit32(1) - ff));

                switch (ii) {
                    case 0:
                        r = v;
                        g = t;
                        b = p;
                        break;
                    case 1:
                        r = q;
                        g = v;
                        b = p;
                        break;
                    case 2:
                        r = p;
                        g = v;
                        b = t;
                        break;
                    case 3:
                        r = p;
                        g = q;
                        b = v;
                        break;
                    case 4:
                        r = t;
                        g = p;
                        b = v;
                        break;
                    default:
                        r = v;
                        g = p;
                        b = q;
                        break;
                }
            }

            rgbData[(i * width + j) * 3 + 0] = static_cast<unsigned char>((int)Posit_floor(r * Posit32(255.0)));
            rgbData[(i * width + j) * 3 + 1] = static_cast<unsigned char>((int)Posit_floor(g * Posit32(255.0)));
            rgbData[(i * width + j) * 3 + 2] = static_cast<unsigned char>((int)Posit_floor(b * Posit32(255.0)));
        }
    }
}

// HSV to RGB 轉換函數 (Posit16_1 版本)
void hsvToRgbPosit16_1Stb(int width, int height, const Posit16_1* hsvData, unsigned char* rgbData) {
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int hsvIndex = (i * width + j) * 3;
            Posit16_1 h = hsvData[hsvIndex + 0];
            Posit16_1 s = hsvData[hsvIndex + 1];
            Posit16_1 v = hsvData[hsvIndex + 2];

            Posit16_1 r = 0, g = 0, b = 0;
            if (s == Posit16_1(0)) {
                r = g = b = v; // achromatic
            } else {
                Posit16_1 hh = h / Posit16_1(60);
                int ii = static_cast<int>(Posit_floor(hh));
                Posit16_1 ff = hh - Posit16_1(ii);
                Posit16_1 p = v * (Posit16_1(1) - s);
                Posit16_1 q = v * (Posit16_1(1) - s * ff);
                Posit16_1 t = v * (Posit16_1(1) - s * (Posit16_1(1) - ff));

                switch (ii) {
                    case 0:
                        r = v;
                        g = t;
                        b = p;
                        break;
                    case 1:
                        r = q;
                        g = v;
                        b = p;
                        break;
                    case 2:
                        r = p;
                        g = v;
                        b = t;
                        break;
                    case 3:
                        r = p;
                        g = q;
                        b = v;
                        break;
                    case 4:
                        r = t;
                        g = p;
                        b = v;
                        break;
                    default:
                        r = v;
                        g = p;
                        b = q;
                        break;
                }
            }

            rgbData[(i * width + j) * 3 + 0] = static_cast<unsigned char>((int)Posit_floor(r * Posit16_1(255.0)));
            rgbData[(i * width + j) * 3 + 1] = static_cast<unsigned char>((int)Posit_floor(g * Posit16_1(255.0)));
            rgbData[(i * width + j) * 3 + 2] = static_cast<unsigned char>((int)Posit_floor(b * Posit16_1(255.0)));
        }
    }
}

// HSV to RGB 轉換函數 (Posit16_2 版本)
void hsvToRgbPosit16_2Stb(int width, int height, const Posit16_2* hsvData, unsigned char* rgbData) {
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int hsvIndex = (i * width + j) * 3;
            Posit16_2 h = hsvData[hsvIndex + 0];
            Posit16_2 s = hsvData[hsvIndex + 1];
            Posit16_2 v = hsvData[hsvIndex + 2];

            Posit16_2 r = 0, g = 0, b = 0;
            if (s == Posit16_2(0)) {
                r = g = b = v; // achromatic
            } else {
                Posit16_2 hh = h / Posit16_2(60);
                int ii = static_cast<int>(Posit_floor(hh));
                Posit16_2 ff = hh - Posit16_2(ii);
                Posit16_2 p = v * (Posit16_2(1) - s);
                Posit16_2 q = v * (Posit16_2(1) - s * ff);
                Posit16_2 t = v * (Posit16_2(1) - s * (Posit16_2(1) - ff));

                switch (ii) {
                    case 0:
                        r = v;
                        g = t;
                        b = p;
                        break;
                    case 1:
                        r = q;
                        g = v;
                        b = p;
                        break;
                    case 2:
                        r = p;
                        g = v;
                        b = t;
                        break;
                    case 3:
                        r = p;
                        g = q;
                        b = v;
                        break;
                    case 4:
                        r = t;
                        g = p;
                        b = v;
                        break;
                    default:
                        r = v;
                        g = p;
                        b = q;
                        break;
                }
            }

            rgbData[(i * width + j) * 3 + 0] = static_cast<unsigned char>((int)Posit_floor(r * Posit16_2(255.0)));
            rgbData[(i * width + j) * 3 + 1] = static_cast<unsigned char>((int)Posit_floor(g * Posit16_2(255.0)));
            rgbData[(i * width + j) * 3 + 2] = static_cast<unsigned char>((int)Posit_floor(b * Posit16_2(255.0)));
        }
    }
}

int main() {
    std::string folderPath = "testImg";

    // 確保 output 資料夾存在
    if (!fs::exists("output")) {
        fs::create_directory("output");
    }

    // 開啟 RMSE 結果輸出檔案
    std::ofstream rmseOutputFile("output/rmse_results.txt", std::ios::app);
    if (!rmseOutputFile.is_open()) {
        std::cerr << "無法開啟 RMSE 輸出檔案: output/rmse_results.txt" << std::endl;
        return 1; // 開啟檔案失敗，終止程式
    }
    
    // 開啟 Exponent 分佈結果輸出檔案
    std::ofstream exponentOutputFile("output/exponent_distribution.txt", std::ios::app);
    if (!exponentOutputFile.is_open()) {
        std::cerr << "無法開啟指數分佈輸出檔案: output/exponent_distribution.txt" << std::endl;
        return 1; // 開啟檔案失敗，終止程式
    }
    

    for (const auto& entry : fs::directory_iterator(folderPath)) {
        if (fs::is_regular_file(entry)) {
            std::string filename = entry.path().string();
            std::string extension = getFileExtension(filename);

            if (extension == "jpg" || extension == "jpeg" || extension == "png" || extension == "bmp" || extension == "tif" || extension == "tiff") {
                int width, height, channels;
                unsigned char* rgbData = stbi_load(filename.c_str(), &width, &height, &channels, 3); // 強制讀取為 3 通道 (RGB)
                if (!rgbData) {
                    std::cerr << "無法讀取影像: " << filename << std::endl;
                    continue; // 讀取失敗，繼續處理下一個檔案
                }
                
                int pixelCount = width * height * 3;

                std::cout << "處理影像: " << filename << ", 寬度: " << width << ", 高度: " << height << std::endl;
                rmseOutputFile << "處理影像: " << filename << ", 寬度: " << width << ", 高度: " << height << std::endl;
                exponentOutputFile << "處理影像: " << filename << ", 寬度: " << width << ", 高度: " << height << std::endl;

                // ***** 獲取原始浮點數結果來進行分析 *****
                double* originalFloatData = new double[pixelCount];
                for(int i = 0; i < pixelCount; ++i) {
                     originalFloatData[i] = static_cast<double>(rgbData[i]) / 255.0;
                }
                analyzeFloatExponentDistribution("原始影像", originalFloatData, pixelCount, exponentOutputFile);

                // ***** RGB to HSV 轉換 *****
                // 浮點數處理
                float* hsvDataFloat = new float[pixelCount];
                rgbToHsvFloatStb(width, height, rgbData, hsvDataFloat);
                double* hsvFloatDoubleData = new double[pixelCount];
                for(int i = 0; i < pixelCount; ++i) {
                    hsvFloatDoubleData[i] = (double)hsvDataFloat[i];
                }
                analyzeFloatExponentDistribution("IEEE_HSV", hsvFloatDoubleData, pixelCount, exponentOutputFile);

                // Posit64 處理 
                Posit64* hsvDataPosit64 = new Posit64[pixelCount];
                rgbToHsvPosit64Stb(width, height, rgbData, hsvDataPosit64);
                double* hsvPosit64DoubleData = new double[pixelCount];
                for(int i = 0; i < pixelCount; ++i) {
                    hsvPosit64DoubleData[i] = (double)hsvDataPosit64[i];
                }
                analyzeFloatExponentDistribution("Posit64_HSV", hsvPosit64DoubleData, pixelCount, exponentOutputFile);

                // Posit32 處理
                Posit32* hsvDataPosit32 = new Posit32[pixelCount];
                rgbToHsvPosit32Stb(width, height, rgbData, hsvDataPosit32);
                double* hsvPosit32DoubleData = new double[pixelCount];
                for(int i = 0; i < pixelCount; ++i) {
                    hsvPosit32DoubleData[i] = (double)hsvDataPosit32[i];
                }
                analyzeFloatExponentDistribution("Posit32_HSV", hsvPosit32DoubleData, pixelCount, exponentOutputFile);

                // Posit16_1 處理
                Posit16_1* hsvDataPosit16_1 = new Posit16_1[pixelCount];
                rgbToHsvPosit16_1Stb(width, height, rgbData, hsvDataPosit16_1);
                double* hsvPosit16_1DoubleData = new double[pixelCount];
                for(int i = 0; i < pixelCount; ++i) {
                    hsvPosit16_1DoubleData[i] = (double)hsvDataPosit16_1[i];
                }
                analyzeFloatExponentDistribution("Posit16_1_HSV", hsvPosit16_1DoubleData, pixelCount, exponentOutputFile);

                // Posit16_2 處理
                Posit16_2* hsvDataPosit16_2 = new Posit16_2[pixelCount];
                rgbToHsvPosit16_2Stb(width, height, rgbData, hsvDataPosit16_2);
                double* hsvPosit16_2DoubleData = new double[pixelCount];
                for(int i = 0; i < pixelCount; ++i) {
                    hsvPosit16_2DoubleData[i] = (double)hsvDataPosit16_2[i];
                }
                analyzeFloatExponentDistribution("Posit16_2_HSV", hsvPosit16_2DoubleData, pixelCount, exponentOutputFile);

                // MPFR 處理
                mpfr_t* hsvDataMpfr = new mpfr_t[pixelCount];
                for(int i = 0; i < pixelCount; ++i) {
                    mpfr_init2(hsvDataMpfr[i], 256); // 初始化每個 mpfr_t 元素
                }
                rgbToHsvMpfrStb(width, height, rgbData, hsvDataMpfr);
                double* hsvMpfrDoubleData = new double[pixelCount];
                for(int i = 0; i < pixelCount; ++i) {
                    hsvMpfrDoubleData[i] = mpfr_get_d(hsvDataMpfr[i], MPFR_RNDN);
                }
                analyzeFloatExponentDistribution("MPFR_HSV", hsvMpfrDoubleData, pixelCount, exponentOutputFile);

                // 將數值轉str並進行誤差運算 (在浮點數層級)
                vector<double> ieeeRMSEVals, pos64RMSEVals, pos32RMSEVals, pos16_1RMSEVals, pos16_2RMSEVals;
                vector<double> MPFRVals;
                // 新增用於計算整數 RMSE 的向量
                vector<double> ieeeIntRMSEVals, pos64IntRMSEVals, pos32IntRMSEVals, pos16_1IntRMSEVals, pos16_2IntRMSEVals;
                
                unsigned char* rgbDataOutMpfrRef = new unsigned char[width * height * 3];
                std::string filenameOnly = fs::path(filename).filename().string();
                std::string filenameWithoutExt = getFilenameWithoutExtension(filenameOnly);
                std::string outputTextFilename = "output/" + filenameWithoutExt + "_hsv_fp_values.txt";
                std::ofstream outputFile(outputTextFilename);
                if (!outputFile.is_open()) {
                    std::cerr << "無法開啟浮點數輸出檔案: " << outputTextFilename << std::endl;
                }

                for (int i = 0; i < pixelCount; i++) {
                    // 將 MPFR 的結果轉換為 double 作為參考值
                    double mpfr_double_val = mpfr_get_d(hsvDataMpfr[i], MPFR_RNDN);
                    MPFRVals.push_back(mpfr_double_val);

                    // 將數值轉str
                    std::string ieeeStr = toString(hsvDataFloat[i]);
                    std::string posit64Str = toString(hsvDataPosit64[i]);
                    std::string posit32Str = toString(hsvDataPosit32[i]);
                    std::string posit16_1Str = toString(hsvDataPosit16_1[i]);
                    std::string posit16_2Str = toString(hsvDataPosit16_2[i]);
                    std::string mpfrStr = toString(hsvDataMpfr[i]);

                    // 將結果寫入檔案
                    outputFile << "Pixel " << i << ":\n";
                    outputFile << "  MPFR: " << mpfrStr << "\n";
                    outputFile << "  IEEE: " << ieeeStr << "\n";
                    outputFile << "  Posit64: " << posit64Str << "\n";
                    outputFile << "  Posit32: " << posit32Str << "\n";
                    outputFile << "  Posit16_1: " << posit16_1Str << "\n";
                    outputFile << "  Posit16_2: " << posit16_2Str << "\n";
                    
                    // 計算浮點數層級的誤差
                    double ieee754ResultDiff = stod(difference(mpfrStr, ieeeStr));
                    double posit64ResultDiff = stod(difference(mpfrStr, posit64Str));
                    double posit32ResultDiff = stod(difference(mpfrStr, posit32Str));
                    double posit16_1ResultDiff = stod(difference(mpfrStr, posit16_1Str));
                    double posit16_2ResultDiff = stod(difference(mpfrStr, posit16_2Str));
                    
                    outputFile << "  IEEE d: " << ieee754ResultDiff << "\n";
                    outputFile << "  Posit64 d: " << posit64ResultDiff << "\n";
                    outputFile << "  Posit32 d: " << posit32ResultDiff << "\n";
                    outputFile << "  Posit16_1 d: " << posit16_1ResultDiff << "\n";
                    outputFile << "  Posit16_2 d: " << posit16_2ResultDiff << "\n";
                    outputFile << "--------------------\n";

                    ieeeRMSEVals.push_back(ieee754ResultDiff);
                    pos64RMSEVals.push_back(posit64ResultDiff);
                    pos32RMSEVals.push_back(posit32ResultDiff);
                    pos16_1RMSEVals.push_back(posit16_1ResultDiff);
                    pos16_2RMSEVals.push_back(posit16_2ResultDiff);
                    
                    // ***** HSV to RGB 轉換後的整數層級誤差計算 *****
                    // 注意：這裡需要先將 HSV 轉換回 RGB，然後再比較整數值
                    // 為了簡化，我們將直接使用 HSV 浮點數值來計算 RMSE，
                    // 如果需要 RGB 整數層級的 RMSE，則需要額外呼叫 hsvToRgb 函式並進行比較。
                    // 這裡先假設您只需要 HSV 浮點數層級的詳細輸出。
                }
                outputFile.close(); // 關閉浮點數結果檔案
                std::cout << "浮點數結果已寫入: " << outputTextFilename << std::endl;


                std::cout << fixed  << setprecision(50) << "--- HSV 浮點數層級 RMSE (vs MPFR) ---" << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "--- HSV 浮點數層級 RMSE (vs MPFR) ---" << std::endl;
                std::cout << fixed  << setprecision(50) << "IEEE fpRMSE:" << RMSE(ieeeRMSEVals) << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "IEEE fpRMSE:" << RMSE(ieeeRMSEVals) << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT64 fpRMSE:" << RMSE(pos64RMSEVals) << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT64 fpRMSE:" << RMSE(pos64RMSEVals) << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT32 fpRMSE:" << RMSE(pos32RMSEVals) << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT32 fpRMSE:" << RMSE(pos32RMSEVals) << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT16_1 fpRMSE:" << RMSE(pos16_1RMSEVals) << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT16_1 fpRMSE:" << RMSE(pos16_1RMSEVals) << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT16_2 fpRMSE:" << RMSE(pos16_2RMSEVals) << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT16_2 fpRMSE:" << RMSE(pos16_2RMSEVals) << std::endl;
                
                // 計算 MRE
                double ieeeFpMRE = calculateMeanRelativeError(ieeeRMSEVals, MPFRVals);
                double posit64FpMRE = calculateMeanRelativeError(pos64RMSEVals, MPFRVals);
                double posit32FpMRE = calculateMeanRelativeError(pos32RMSEVals, MPFRVals);
                double posit16_1FpMRE = calculateMeanRelativeError(pos16_1RMSEVals, MPFRVals);
                double posit16_2FpMRE = calculateMeanRelativeError(pos16_2RMSEVals, MPFRVals);

                std::cout << fixed  << setprecision(50) << "--- HSV 浮點數層級 MRE (vs MPFR) ---" << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "--- HSV 浮點數層級 MRE (vs MPFR) ---" << std::endl;
                std::cout << fixed  << setprecision(50) << "IEEE FpMRE:" << ieeeFpMRE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "IEEE FpMRE:" << ieeeFpMRE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT64 FpMRE:" << posit64FpMRE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT64 FpMRE:" << posit64FpMRE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT32 FpMRE:" << posit32FpMRE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT32 FpMRE:" << posit32FpMRE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT16_1 FpMRE:" << posit16_1FpMRE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT16_1 FpMRE:" << posit16_1FpMRE << std::endl;
                std::cout << fixed  << setprecision(50) << "POSIT16_2 FpMRE:" << posit16_2FpMRE << std::endl;
                rmseOutputFile << fixed  << setprecision(50) << "POSIT16_2 FpMRE:" << posit16_2FpMRE << std::endl;
                rmseOutputFile << "----------------------------------------\n";


                // ***** HSV to RGB 轉換並儲存影像 *****
                unsigned char* rgbDataOutFloat = new unsigned char[pixelCount];
                hsvToRgbFloatStb(width, height, hsvDataFloat, rgbDataOutFloat);

                unsigned char* rgbDataOutPosit64 = new unsigned char[pixelCount];
                hsvToRgbPosit64Stb(width, height, hsvDataPosit64, rgbDataOutPosit64);

                unsigned char* rgbDataOutPosit32 = new unsigned char[pixelCount];
                hsvToRgbPosit32Stb(width, height, hsvDataPosit32, rgbDataOutPosit32);

                unsigned char* rgbDataOutPosit16_1 = new unsigned char[pixelCount];
                hsvToRgbPosit16_1Stb(width, height, hsvDataPosit16_1, rgbDataOutPosit16_1);

                unsigned char* rgbDataOutPosit16_2 = new unsigned char[pixelCount];
                hsvToRgbPosit16_2Stb(width, height, hsvDataPosit16_2, rgbDataOutPosit16_2);


                std::string outputFilenameFloat = "output/" + filenameWithoutExt + "_hsv_float." + extension;
                stbi_write_png(outputFilenameFloat.c_str(), width, height, 3, rgbDataOutFloat, width * 3);

                std::string outputFilenamePosit64 = "output/" + filenameWithoutExt + "_hsv_Posit64." + extension;
                stbi_write_png(outputFilenamePosit64.c_str(), width, height, 3, rgbDataOutPosit64, width * 3);

                std::string outputFilenamePosit32 = "output/" + filenameWithoutExt + "_hsv_Posit32." + extension;
                stbi_write_png(outputFilenamePosit32.c_str(), width, height, 3, rgbDataOutPosit32, width * 3);

                std::string outputFilenamePosit16_1 = "output/" + filenameWithoutExt + "_hsv_Posit16_1." + extension;
                stbi_write_png(outputFilenamePosit16_1.c_str(), width, height, 3, rgbDataOutPosit16_1, width * 3);

                std::string outputFilenamePosit16_2 = "output/" + filenameWithoutExt + "_hsv_Posit16_2." + extension;
                stbi_write_png(outputFilenamePosit16_2.c_str(), width, height, 3, rgbDataOutPosit16_2, width * 3);

                // ***** 記憶體釋放 *****
                stbi_image_free(rgbData);
                delete[] originalFloatData;
                delete[] hsvDataFloat;
                delete[] hsvFloatDoubleData;
                delete[] hsvDataPosit64;
                delete[] hsvPosit64DoubleData;
                delete[] hsvDataPosit32;
                delete[] hsvPosit32DoubleData;
                delete[] hsvDataPosit16_1;
                delete[] hsvPosit16_1DoubleData;
                delete[] hsvDataPosit16_2;
                delete[] hsvPosit16_2DoubleData;
                for (int i = 0; i < pixelCount; ++i) {
                    mpfr_clear(hsvDataMpfr[i]);
                }
                delete[] hsvDataMpfr;
                delete[] hsvMpfrDoubleData;
                delete[] rgbDataOutFloat;
                delete[] rgbDataOutPosit64;
                delete[] rgbDataOutPosit32;
                delete[] rgbDataOutPosit16_1;
                delete[] rgbDataOutPosit16_2;

                std::cout << "浮點數 HSV 影像已儲存為: " << outputFilenameFloat << std::endl;
                std::cout << "Posit64 HSV 影像已儲存為: " << outputFilenamePosit64 << std::endl;
                std::cout << "Posit32 HSV 影像已儲存為: " << outputFilenamePosit32 << std::endl;
                std::cout << "Posit16_1 HSV 影像已儲存為: " << outputFilenamePosit16_1 << std::endl;
                std::cout << "Posit16_2 HSV 影像已儲存為: " << outputFilenamePosit16_2 << std::endl;
            }
        }
    }

    // 關閉 RMSE 檔案 
    rmseOutputFile.close();
    // 關閉 exponent 檔案
    exponentOutputFile.close();

    return 0;
}
