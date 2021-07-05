# Download and attempt to compile libf2c.a
wget http://www.netlib.org/f2c/libf2c.zip
unzip libf2c.zip
mv makefile.u Makefile
sed -i 's/CFLAGS = -O/CFLAGS = -O -fPIC/' Makefile
make 
cp libf2c.a /usr/lib64/