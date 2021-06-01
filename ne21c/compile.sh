c++ -L./ -O3 -shared -fPIC -std=c++11 $(python3 -m pybind11 --includes)  -Wall -o ne21c$(python3-config --extension-suffix) main.cpp -lf2c -lm
