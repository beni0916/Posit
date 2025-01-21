CXX = g++
CXXFLAGS = -std=c++17

KER_SIN = test0
KER_COS = test1
KER_TAN = test2
COSPI = test5
SINPI = test11
REMH = test12

SRC_T0 = test0.cpp
SRC_T1 = test1.cpp
SRC_T2 = test2.cpp
SRC_T5 = test5.cpp
SRC_T11 = test11.cpp
SRC_T12 = test12.cpp

SRC_SINPI = Posit_sinPi.cpp
SRC_COS = Posit_cos.cpp
SRC_COSPI = Posit_cosPi.cpp
SRC_KER_SIN = __kernel_sin.cpp
SRC_KER_COS = __kernel_cos.cpp
SRC_KER_TAN = __kernel_tan.cpp
SRC_FABS = Posit_fabs.cpp
SRC_REMP = Posit_rempio2.cpp
SRC_REMH = Posit_remhalf.cpp

all: $(KER_SIN) $(KER_COS) $(KER_TAN) $(COSPI) $(SINPI) $(REMH)

$(KER_SIN):
	$(CXX) $(CXXFLAGS) -o $(KER_SIN) $(SRC_T0) $(SRC_KER_SIN)
	
$(KER_COS):
	$(CXX) $(CXXFLAGS) -o $(KER_COS) $(SRC_T1) $(SRC_KER_COS)
	
$(KER_TAN):
	$(CXX) $(CXXFLAGS) -o $(KER_TAN) $(SRC_T2) $(SRC_KER_TAN) $(SRC_FABS)

$(COSPI):
	$(CXX) $(CXXFLAGS) -o $(COSPI) $(SRC_T5) $(SRC_COSPI) $(SRC_COS) $(SRC_KER_SIN) $(SRC_KER_COS) $(SRC_FABS) $(SRC_REMH) $(SRC_REMP)

$(SINPI): 
	$(CXX) $(CXXFLAGS) -o $(SINPI) $(SRC_T11) $(SRC_SINPI) $(SRC_KER_SIN) $(SRC_KER_COS) $(SRC_FABS) $(SRC_REMH)
	
$(REMH):
	$(CXX) $(CXXFLAGS) -o $(REMH) $(SRC_T12) $(SRC_REMH) $(SRC_REMP)
	
clean:
	rm -f $(KER_SIN) $(KER_COS) $(KER_TAN) $(COSPI) $(SINPI) $(REMH)

.PHONY: all clean
