echo --- f2c -C++ *.f ---
f2c -C++ *.f
echo ---  cc -o ne2001 *.c -lf2c -lm ---
cc -o ne2001 *.c -lf2c -lm
