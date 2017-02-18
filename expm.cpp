#include <boost/random.hpp>
#include <boost/assert.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>
#include <boost/numeric/ublas/expm.hpp>
using namespace boost::numeric::ublas;
using namespace std;
int main(void)
{
        matrix<double> gen(3,3);
        gen(0,0) = 1  ; gen(0,1) = 1; gen(0,2) = 0;
        gen(1,0) = 0; gen(1,1) = 0   ; gen(1,2) = 2;
        gen(2,0) = 0  ; gen(2,1) = 0   ; gen(2,2) = -1;

        cout<< expm_pad(gen) <<"\n\n";
        return 0;
}

