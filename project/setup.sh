git submodule update --init --recursive

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

# graphmap
cd graphmap2
make modules  
make  

cd ..

#desalt
cd deSALT/src/deBGA-master/
make   ## built deBGA for RdBG-index
cd ..
make   ## built deSALT for alignment
