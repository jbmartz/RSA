#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <gmp.h>

int gcd(int a, int b)
{
    if (a == 0)
        return b;
    return gcd(b % a, a);
}

int calc_inverse(int mod, int num) {
    num %= mod;

    int t = 0;
    int newt = 1;
    int r = mod;
    int newr = num;
    int quo;

    while (newr != 0) {
        quo = r / newr;
        int temp = newr;
        newr = r - quo * newr;
        r = temp;
        temp = newt;
        newt = t - quo * newt;
        t = temp;
    }

    if (r > 1) {
        std::cout << "Not invertible\n";
        return -1;
    }
    if (t < 0) {
        t += mod;
    }

    return t;
}
int modular_exponentiation(int b, int exp, int mod)
{
    int res = 1;
    while (exp > 0)
    {
        if (exp % 2 == 1)
            res = (res * b) % mod;
        exp = exp >> 1;
        b = (b * b) % mod;
    }
    return res;
}
int choose_e(int n, int phi_n) {
    int e;

    do {
        e = phi_n + rand() % ( n + 1 ) - phi_n;
    }
    while(gcd(phi_n, e) != 1);
    return e;
}

int totient(int p, int q)
{
   return (p - 1) * (q - 1);
}



bool is_prime(int n, int repeat) {

    int a;
    for(int i = 1; i <= repeat; i++) {
        while((a =  2 + rand() % (( n + 1 ) - 2)) % 2 == 0);
        if(gcd(n, a) != 1)
            return false;


//        if((a^(n - 1) % n) == 1)
//            return false;

        int exponent = pow(a, n - 1);
        if(exponent % n == 1)
            return false;
    }

    return true;
}
std::pair<int, int> choose_primes() {
    std::pair<int, int> pq;

    int max = 200, min = 3;

        do {
            while ((pq.first = (rand() % ((max - min) + 1) + min)) % 2 ==
                   0);       // std::cout << pq.first << std::endl;
        } while (!is_prime(pq.first, 100));

        do {
            while ((pq.second = (rand() % ((max - min) + 1) + min)) % 2 == 0);
        } while (!is_prime(pq.second, 100));


    return pq;

}

int encrypt(std::pair<int, int> pub_key, int plaintext) {
    return modular_exponentiation(plaintext, pub_key.first, pub_key.second);

}

int decrypt(std::pair<int, int> priv_key, int ciphertext) {
    return modular_exponentiation(ciphertext, priv_key.first, priv_key.second);
}

std::string convert_to_base26(std::string str) {

    std::string new_str;
    for (int i = 0; i < str.length() - 2; i += 2) {
        int sum = 0;
        int k = 2;

        for (int j = i; j < i + 3; j++) {
            int n = (((int) str[j] - 'a') * pow(26, k));

            sum += n;
            --k;
        }




        new_str.push_back((sum % 26) + 'a');
    }
    std::cout << new_str << std::endl;


    return new_str;


}

std::string convert_from_base26(std::string str) {

    std::string new_str;
    int a, b, c;
    for(int i = 0; i < str.length(); i++) {
        int k = (int) str[i] / 26;

    }

    return "";
}



int main() {

    std::srand(std::time(nullptr));
    std::pair<int, int> pq = choose_primes();
    int n = pq.first * pq.second;
    int phi_n = totient(pq.first, pq.second);
    int e = choose_e(n, phi_n);
    int d = calc_inverse(phi_n, e);

    std::pair<int, int> pub_key(e, n);
    std::pair<int, int> priv_key(d, n);

    std::cout << "\np: " << pq.first << "\nq: " << pq.second << std::endl;
    std::cout << "e: " << e << std::endl;
    std::cout << "d: " << d << std::endl;
    std::cout << "n: " << n << std::endl;
    std::cout << "phi_n: " << phi_n << std::endl;
    std::string m = "this message is to be encrypted";

  //std::string bet = "bet";

    m = convert_to_base26(m);

    std::cout << "encrypted text: ";
    std::string ciphertext = "";
    int cipher[m.size()];
    for(int i = 0; i < m.size(); i++) {
        int a = encrypt(pub_key, (int) m[i]);
        cipher[i] = a;

       std::cout << (char)((a % 26) + 'a');

    }
    std::cout << "\ndecrypted text: ";

    for(int i = 0; i < m.size(); i++) {
        int b = decrypt(priv_key, cipher[i]);
        std::cout << (char) b;

    }


    std::cout << std::endl;

  // std::cout << (int) ((799 / pow(26, 2)) / 26);
 // convert_from_base26(4750);



    return 0;
}