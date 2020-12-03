//
// Created by emery on 19/11/2020.
//
#include <iostream>
using namespace std;
template <int N>
void fonctionTest(const unsigned char (&test)[N]){
    cout<<N<<"\n";
    for (int i = 0; i <N; ++i) {
        cout<<test[i]<<"\n";
    }
}

int main() {
    unsigned char test[3];
    test[0]= 'b';
    test[1]='c';
    test[2]='a';
    fonctionTest(test);


}