#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <gmpxx.h>
#include <random>

#define KEY_SIZE 2048
std::random_device rd;

mpz_class gcd_mpz(mpz_class &a, mpz_class &b) {
    if (a == 0)
        return b;
    return gcd(b % a, a);
}

void calc_inverse(mpz_class &inverse, mpz_class &mod, mpz_class &num) {
    num %= mod;

    mpz_class t = 0;
    mpz_class newt = 1;
    mpz_class r = mod;
    mpz_class newr = num;
    mpz_class quo;

    while (newr != 0) {
        quo = r / newr;
        mpz_class temp = newr;
        newr = r - quo * newr;
        r = temp;
        temp = newt;
        newt = t - quo * newt;
        t = temp;
    }

    if (r > 1) {
        std::cout << "Not invertible\n";
        inverse = -1;
        return;
    }
    if (t < 0) {
        t += mod;
    }

    inverse = t;
}

mpz_class modular_exponentiation_mpz(mpz_class &b, mpz_class exp, mpz_class &mod)
{
    mpz_class res = 1;
    while (exp > 0)
    {
        if (exp % 2 == 1)
            res = (res * b) % mod;
        exp = exp >> 1;
        b = (b * b) % mod;
    }
    return res;
}
void choose_e(mpz_class &e, mpz_class &n, mpz_class &phi_n) {
    gmp_randclass r(gmp_randinit_mt);
    r.seed(rd());
    do {
        //change to less than phi_n
      //  e = r.get_z_bits(1024);

        e = r.get_z_range(phi_n);
    }
    while(gcd_mpz(phi_n, e) != 1);
}

void totient(mpz_class &phi_n, mpz_class &p, mpz_class &q)
{
    phi_n = (p - 1) * (q - 1);
}



bool is_prime(mpz_class &n, int repeat) {

    gmp_randclass r(gmp_randinit_mt);



    mpz_class a, gcd, one, exponent;

    gmp_randstate_t rand_state;


    for (int i = 1; i <= repeat; i++) {

       a = r.get_z_range(n - 2);



        if(gcd_mpz(n, a) == 1)
            return false;

        if(modular_exponentiation_mpz(a, n - 1, n) == 1)
            return false;


    }
    return true;

}

bool pre_filter_is_prime(mpz_class &val) {
    for(int i = 2; i < 100; i++) {
        if((val % i) == 0) {
            return false;
        }
    }

    return true;
}

//if x is congruent to 0, 2, 3, or 4 mod 6 not prime

void choose_prime_mpz(mpz_class &val) {
    gmp_randclass r(gmp_randinit_mt);
    r.seed(rd());


    do {
        val = r.get_z_bits(KEY_SIZE);
    } //while(!pre_filter_is_prime(val) && !is_prime(val, 100000));


    while (mpz_probab_prime_p(val.get_mpz_t(), 100) == 0);


}


mpz_class encrypt(std::pair<mpz_class, mpz_class> &pub_key, mpz_class &plaintext) {
    return modular_exponentiation_mpz(plaintext, pub_key.first, pub_key.second);

}

mpz_class decrypt(std::pair<mpz_class, mpz_class> &priv_key, mpz_class &ciphertext) {
    return modular_exponentiation_mpz(ciphertext, priv_key.first, priv_key.second);
}


int main() {

    mpz_class p_val, q_val, n_val, e_val, d_val, phi_n;
    choose_prime_mpz(p_val);
    choose_prime_mpz(q_val);
    n_val = p_val * q_val;
    totient(phi_n, p_val, q_val);
    choose_e(e_val, n_val, phi_n);
    calc_inverse(d_val, phi_n, e_val);

    std::cout << "p: " << p_val.get_str() << std::endl << std::endl;
    std::cout << "q: " << q_val.get_str() << std::endl << std::endl;
    std::cout << "n: " << n_val.get_str() << std::endl << std::endl;
    std::cout << "e: " << e_val.get_str() << std::endl << std::endl;
    std::cout << "d: " << d_val.get_str() << std::endl << std::endl;




    std::pair<mpz_class, mpz_class> pub_key(e_val, n_val);
    std::pair<mpz_class, mpz_class> priv_key(d_val, n_val);


    std::string m = "Massive numbers RSA!";
    std::cout << "Encrypted text: ";
    std::string ciphertext = "";
    mpz_class cipher[m.size()];
    for(int i = 0; i < m.size(); i++) {
        mpz_class convert = (int) m[i];
        mpz_class a = encrypt(pub_key, convert);
        cipher[i] = a;


        a = a % 26;
        std::cout << (char) (a.get_ui() + 'a');


    }


    std::cout << "\nEncrypted text: ";

    for(int i = 0; i < m.size(); i++) {
        mpz_class b = decrypt(priv_key, cipher[i]);
        std::cout << (char) b.get_ui();


    }

    std::cout << std::endl;

    return 0;
}
