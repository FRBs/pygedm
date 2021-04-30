echo --- f2c *.f ---
f2c *.f
echo ---  cc -o ne2001 *.c -lf2c -lm ---
cc -o ne2001 *.c -lf2c -lm
