echo --- f2c -C++ *.f ---
f2c -C++ *.f
echo ---  cc -o ne2001 *.c -lf2c -lm ---
cc -o ne2001 *.c -lf2c -lm

# cat *.f | f2c > ne21.c
# gcc -c -Wall -Werror -fPIC ne21.c
# gcc -shared -o libne21.so ne21.o
#c++ -L./ -O3 -shared -fPIC -std=c++11 $(python3 -m pybind11 --includes)  -Wall -o example$(python3-config --extension-suffix) main.c -lne21 -lf2c -lm
