# minimap
cd uLTRA

cd minimap2
make clean
make
cd ..

# slamem
cd slaMEM
make clean
make
cd ..

# mummer
rm -rf mummer-4.0.0rc1
tar -xvf mummer-4.0.0rc1.tar
cd mummer-4.0.0rc1
./configure
make

cd ..

#desalt
installed with conda