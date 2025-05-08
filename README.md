# PMatch

1. Datasets can be downloaded at: https://drive.google.com/drive/folders/1pwlrT5_nOIRown6lLKSyS74t3rf5PR2_?usp=sharing or https://snap.stanford.edu/data/index.html

2. install g++ "brew install g++"

3. install GMP 
   1) download gmp-6.1.2 and unzip "xz -d gmp-6.1.2.tar.xz", "tar xvf gmp-6.1.2.tar "
   2) "cd gmp-6.1.2"
   2) "./configure --prefix=${HOME}/gmp/6.1.2"
   3) "make"
   4) "make check"
   5) "make install"

(or "brew install gmp")

4. Use "Preprocessing.sh" to generate queries and to build the path tables and NL-Index

5. Use "overall-test.sh" to evaluate PMatch

