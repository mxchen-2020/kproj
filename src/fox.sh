tar zxvf fox.tgz 
cd fox
./configure FC=ifort 
make
value1=`./FoX-config --fcflags`
value2=`./FoX-config --libs`
cd ..
sed -i      " 36i #================Fox_xml=================================="     Makefile
sed -i      " 37i  finclude =   $value1 "                                         Makefile
sed -i      " 38i  foxlib   =   $value2 "                                         Makefile 
sed -i      " 39i #========================================================="     Makefile
