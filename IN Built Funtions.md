## min_element/ max_element
```
int sm = *min_element(a.begin() + l, a.begin() + r + 1); // r is non inclusive so, added 1
int mx = *max_element(a.begin() + l, a.begin() + r + 1);
```
## Count set bits
```
    __builtin_popcount(int number);
```
## GCD
```
__gcd();
```
## Distnace function : distance between two iterator position
```
distance();
```
## catalan number

## Tuples in C++
```
#include <iostream>
#include <tuple>
using namespace std;

int main() {
  
    // Declaring tuple
    tuple<char, int, float> geek;

    // Assigning values to tuple using make_tuple()
    geek = make_tuple('a', 10, 15.5);

    // Printing initial tuple values using get()
    cout << "The initial values of tuple are: ";
    cout << get<0>(geek) << " " << get<1>(geek);
    cout << " " << get<2>(geek) << endl;

    // Use of get() to change values of tuple
    get<0>(geek) = 'b';
    get<2>(geek) = 20.5;

    // Printing modified tuple values
    cout << "The modified values of tuple are: ";
    cout << get<0>(geek) << " " << get<1>(geek);
    cout << " " << get<2>(geek) << endl;

    return 0;
}
```
