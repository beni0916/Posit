CXX = g++
CXXFLAGS = -std=c++17
LIBFLAGS = -lgmp -lmpfr

RES_COSPI = Result_cosPi
RES_SINPI = Result_sinPi
RES_TANPI = Result_tanPi

RES_COS = Result_cos
RES_SIN = Result_sin
RES_TAN = Result_tan

RES_ASINPI = Result_asinPi
RES_ACOSPI = Result_acosPi
RES_ATANPI = Result_atanPi

RES_EXP = Result_exp
RES_EXPM1 = Result_expM1
RES_EXP2 = Result_exp2
RES_EXP2M1 = Result_exp2M1
RES_EXP10 = Result_exp10
RES_EXP10M1 = Result_exp10M1

RES_RSQRT = Result_rSqrt

SRC_RES_COSPI = Result_cosPi.cpp
SRC_RES_SINPI = Result_sinPi.cpp
SRC_RES_TANPI = Result_tanPi.cpp

SRC_RES_COS = Result_cos.cpp
SRC_RES_SIN = Result_sin.cpp
SRC_RES_TAN = Result_tan.cpp

SRC_RES_ASINPI = Result_asinPi.cpp
SRC_RES_ACOSPI = Result_acosPi.cpp
SRC_RES_ATANPI = Result_atanPi.cpp

SRC_RES_EXP = Result_exp.cpp
SRC_RES_EXPM1 = Result_expM1.cpp
SRC_RES_EXP2 = Result_exp2.cpp
SRC_RES_EXP2M1 = Result_exp2M1.cpp
SRC_RES_EXP10 = Result_exp10.cpp
SRC_RES_EXP10M1 = Result_exp10M1.cpp

SRC_RES_RSQRT = Result_rSqrt.cpp

SRC_SINPI = Posit_sinPi.cpp
SRC_COSPI = Posit_cosPi.cpp
SRC_TANPI = Posit_tanPi.cpp

SRC_SIN = Posit_sin.cpp
SRC_COS = Posit_cos.cpp
SRC_TAN = Posit_tan.cpp

SRC_ASINPI = Posit_asinPi.cpp
SRC_ACOSPI = Posit_acosPi.cpp
SRC_ATANPI = Posit_atanPi.cpp

SRC_KER_SIN = __kernel_sin.cpp
SRC_KER_COS = __kernel_cos.cpp
SRC_KER_TAN = __kernel_tan.cpp

SRC_EXP = Posit_exp.cpp
SRC_EXPM1 = Posit_expMinus1.cpp
SRC_EXP2 = Posit_exp2.cpp
SRC_EXP2M1 = Posit_exp2Minus1.cpp
SRC_EXP10 = Posit_exp10.cpp
SRC_EXP10M1 = Posit_exp10Minus1.cpp

SRC_FABS = Posit_fabs.cpp
SRC_REMP = Posit_rempio2.cpp
SRC_REMH = Posit_remhalf.cpp
SRC_SQRT = Posit_sqrt.cpp
SRC_RSQRT = Posit_rSqrt.cpp
SRC_FLOOR = Posit_floor.cpp

all: $(RES_COSPI) $(RES_SINPI) $(RES_TANPI) $(RES_ASINPI) $(RES_ACOSPI) $(RES_ATANPI) $(RES_EXP) $(RES_EXPM1) $(RES_EXP2) $(RES_EXP2M1) $(RES_EXP10) $(RES_EXP10M1) $(RES_RSQRT)

$(RES_COSPI):
	$(CXX) $(CXXFLAGS) -o $(RES_COSPI) $(SRC_RES_COSPI) $(SRC_COSPI) $(SRC_KER_SIN) $(SRC_KER_COS) $(SRC_FABS) $(SRC_REMH) $(LIBFLAGS)

$(RES_SINPI):
	$(CXX) $(CXXFLAGS) -o $(RES_SINPI) $(SRC_RES_SINPI) $(SRC_SINPI) $(SRC_KER_SIN) $(SRC_KER_COS) $(SRC_FABS) $(SRC_REMH) $(LIBFLAGS)
	
$(RES_TANPI):
	$(CXX) $(CXXFLAGS) -o $(RES_TANPI) $(SRC_RES_TANPI) $(SRC_TANPI) $(SRC_KER_TAN) $(SRC_FABS) $(SRC_REMH) $(LIBFLAGS)

$(RES_ASINPI):
	$(CXX) $(CXXFLAGS) -o $(RES_ASINPI) $(SRC_RES_ASINPI) $(SRC_ASINPI) $(SRC_FABS) $(SRC_SQRT) $(LIBFLAGS)

$(RES_ACOSPI):
	$(CXX) $(CXXFLAGS) -o $(RES_ACOSPI) $(SRC_RES_ACOSPI) $(SRC_ACOSPI) $(SRC_FABS) $(SRC_SQRT) $(LIBFLAGS)

$(RES_ATANPI):
	$(CXX) $(CXXFLAGS) -o $(RES_ATANPI) $(SRC_RES_ATANPI) $(SRC_ATANPI) $(SRC_FABS) $(SRC_SQRT) $(LIBFLAGS)

$(RES_EXP):
	$(CXX) $(CXXFLAGS) -o $(RES_EXP) $(SRC_RES_EXP) $(SRC_EXP) $(LIBFLAGS)

$(RES_EXPM1):
	$(CXX) $(CXXFLAGS) -o $(RES_EXPM1) $(SRC_RES_EXPM1) $(SRC_EXPM1) $(LIBFLAGS)
	
$(RES_EXP2):
	$(CXX) $(CXXFLAGS) -o $(RES_EXP2) $(SRC_RES_EXP2) $(SRC_EXP2) $(SRC_FLOOR) $(SRC_FABS) $(LIBFLAGS)

$(RES_EXP2M1):
	$(CXX) $(CXXFLAGS) -o $(RES_EXP2M1) $(SRC_RES_EXP2M1) $(SRC_EXP2M1) $(SRC_FLOOR) $(SRC_FABS) $(LIBFLAGS)
	
$(RES_EXP10):
	$(CXX) $(CXXFLAGS) -o $(RES_EXP10) $(SRC_RES_EXP10) $(SRC_EXP10) $(SRC_EXP) $(LIBFLAGS)

$(RES_EXP10M1):
	$(CXX) $(CXXFLAGS) -o $(RES_EXP10M1) $(SRC_RES_EXP10M1) $(SRC_EXP10M1) $(SRC_EXP) $(LIBFLAGS)

$(RES_RSQRT):
	$(CXX) $(CXXFLAGS) -o $(RES_RSQRT) $(SRC_RES_RSQRT) $(SRC_RSQRT) $(LIBFLAGS)
	
clean:
	rm -f $(RES_COSPI) $(RES_SINPI) $(RES_TANPI) $(RES_ASINPI) $(RES_ACOSPI) $(RES_ATANPI) $(RES_EXP) $(RES_EXPM1) $(RES_EXP2) $(RES_EXP2M1) $(RES_EXP10) $(RES_EXP10M1) $(RES_RSQRT)

