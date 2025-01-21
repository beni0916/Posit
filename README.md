test0: __kernel_sin.cpp
test1: __kernel_cos.cpp
test2: __kernel_tan.cpp
test5: Posit_cosPi.cpp
test11: Posit_sinPi.cpp
test12: Posit_remhalf.cpp

compile: make test#
ex: make test0

compile all: make all
clean all: make clean
