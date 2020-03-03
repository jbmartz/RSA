#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <math.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <gmpxx.h>
#include <random>
#include <limits>
#include <iomanip>
#include <fstream>

#define PROBAB_PRIME_100 "99.999999999999999999999999999999211139094778988194588271434717213770326793564890976995229721069335938" //probability of prime is (1/2)^k. I do k = 100 repetitions
#define KEY_SIZE 1024 //size of p and q - will give 2048 bit rsa key
#define BLOCK_SIZE 256 //block size to pad to
#define MAX_DATA_SIZE 214 //bytes in block that can hold data
#define REPS 100 //number of reps for solovay strassen


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

mpz_class modular_exponentiation_mpz(mpz_class &b, mpz_class exp, mpz_class &mod) {
    mpz_class res = 1;
    while (exp > 0) {
        if ((exp % 2) == 1)
            res = (res * b) % mod;
        exp = exp >> 1;
        b = (b * b) % mod;
    }
    return res;
}

void choose_e(gmp_randclass &r, mpz_class &e, mpz_class &n, mpz_class &phi_n) {
    do {
        e = r.get_z_range(phi_n);
    } while (gcd_mpz(phi_n, e) != 1);
}

void totient(mpz_class &phi_n, mpz_class &p, mpz_class &q) {
    phi_n = (p - 1) * (q - 1);
}


mpz_class jacobi_symbol(mpz_class a, mpz_class b) {
    int res = 1;

    if (!a)
        return 0;

    if (a < 0) {
        a *= -1;
        if (b % 4 == 3)
            res *= -1;
    }

    if (a == 1)
        return res;

    while (a) {
        if (a < 0) {
            a *= -1;
            if (b % 4 == 3)
                res *= -1;
        }

        while (a % 2 == 0) {
            a = a / 2;
            if (b % 8 == 3 || b % 8 == 5)
                res *= -1;

        }

        std::swap(a, b);

        if (a % 4 == 3 && b % 4 == 3)
            res *= -1;
        a = a % b;

        if (a > b / 2)
            a = a - b;

    }

    if (b == 1)
        return res;

    return 0;
}

bool solovay_strassen_primality_test(gmp_randclass &r, mpz_class n, int repetitions) {
    if (n < 2)
        return false;
    if (n != 2 && n % 2 == 0)
        return false;

    while (repetitions-- > 0) {
        mpz_class a = r.get_z_range(n - 2) + 2;
        mpz_class jac_sym = (n + jacobi_symbol(a, n)) % n;

        if (!jac_sym || modular_exponentiation_mpz(a, (n - 1) / 2, n) != jac_sym)
            return false;
    }

    return true;
}

bool fermats_primaility_test(gmp_randclass &r, mpz_class &n, int repeat) {

    mpz_class a;

    for (int i = 1; i <= repeat; i++) {
        a = r.get_z_range(n - 2) + 2;


        //if gcd != 1 and modular_exponentation != 1
        if (gcd_mpz(n, a) != 1)
            return false;

        if (modular_exponentiation_mpz(a, n - 1, n) != 1)
            return false;

    }
    return true;

}


void choose_prime_mpz(gmp_randclass &r, mpz_class &val, int repetitions) {

    int prime_count = 0;
    do {
        val = r.get_z_bits(KEY_SIZE);
        ++prime_count;
    } while (!solovay_strassen_primality_test(r, val, repetitions));

    std::cout << prime_count << " random numbers tried before a finding prime with " << repetitions
              << " values used. The probability of primality is " << PROBAB_PRIME_100 << std::endl;

}


mpz_class encrypt(std::pair<mpz_class, mpz_class> &pub_key, mpz_class &plaintext) {
    return modular_exponentiation_mpz(plaintext, pub_key.first, pub_key.second);

}

mpz_class decrypt(std::pair<mpz_class, mpz_class> &priv_key, mpz_class &ciphertext) {
    return modular_exponentiation_mpz(ciphertext, priv_key.first, priv_key.second);
}

bool write_file(std::string filename, std::pair<mpz_class, mpz_class> &key) {
    std::ofstream myfile(filename);
    if (myfile.is_open()) {
        myfile << key.first << "\n" << key.second;
        myfile.flush();
        myfile.close();
        return true;
    } else {
        std::cerr << "Error writing file" << std::endl;
        return false;
    }
}


bool write_ciphertext_to_file(std::string filename, std::string ciphertext) {
    std::ofstream myfile(filename);
    if (myfile.is_open()) {
        myfile << ciphertext;
        myfile.flush();
        myfile.close();
        return true;
    } else {
        std::cerr << "Error writing file." << std::endl;
        return false;
    }
}

std::string read_file(std::string &file_string) {
    std::ifstream t(file_string);
    if (t.is_open()) {
        std::stringstream buffer;
        buffer << t.rdbuf();
        return buffer.str();
    } else {
        std::cerr << "File " << file_string << " does not exist." << std::endl;
        return "";
    }

}


//remove padding to reform string
void join_blocks(const char *plaintext_padded, char *plaintext, int n_blocks) {
    for (int i = 0; i < n_blocks; i++) {
        memcpy(&plaintext[i * MAX_DATA_SIZE], &plaintext_padded[i * BLOCK_SIZE], MAX_DATA_SIZE);
    }

}

//split blocks into 214 byte blocks and pad as necessary
void split_blocks(char *plaintext_padded, std::string string_split, int n_blocks) {

    for (int i = 0; i < n_blocks; i++) {
        std::string substr = string_split.substr(i * MAX_DATA_SIZE, MAX_DATA_SIZE);
        memcpy(&plaintext_padded[i * BLOCK_SIZE], substr.data(), MAX_DATA_SIZE);
    }
}

//encrypt padded blocks
void encrypt_blocks(char *ciphertext, char *plaintext_padded, int n_blocks, std::pair<mpz_class, mpz_class> &pub_key) {

    mpz_t encode_list[n_blocks];
    mpz_class encode_wrapper[n_blocks];

    for (int i = 0; i < n_blocks; i++) {
        mpz_init(encode_list[i]);
        mpz_import(encode_list[i], 1, -1, BLOCK_SIZE, -1, 0, &plaintext_padded[i * BLOCK_SIZE]);
        mpz_class encode_wrap(encode_list[i]);
        encode_wrapper[i] = encrypt(pub_key, encode_wrap);
        mpz_export(&ciphertext[i * BLOCK_SIZE], NULL, -1, BLOCK_SIZE, -1, 0, encode_wrapper[i].get_mpz_t());

    }


}


//decrypt padded blocks
void decrypt_blocks(char *decrypted_text, char *ciphertext, int n_blocks, std::pair<mpz_class, mpz_class> &priv_key) {


    mpz_t decode_list[n_blocks];
    for (int i = 0; i < n_blocks; i++) {
        mpz_init(decode_list[i]);
        mpz_import(decode_list[i], 1, -1, BLOCK_SIZE, -1, 0, &ciphertext[i * BLOCK_SIZE]);

        mpz_class decode_wrap(decode_list[i]);

        mpz_export(&decrypted_text[i * BLOCK_SIZE], NULL, -1, MAX_DATA_SIZE, -1, 0,
                   decrypt(priv_key, decode_wrap).get_mpz_t());
    }


}


//generate RSA key
void gen_key(std::string &key_file_name) {

    std::random_device rd;
    gmp_randclass r(gmp_randinit_mt);
    r.seed(rd());
    mpz_class p_val, q_val, n_val, e_val, d_val, phi_n;
    choose_prime_mpz(r, p_val, REPS);
    choose_prime_mpz(r, q_val, REPS);
    n_val = p_val * q_val;
    totient(phi_n, p_val, q_val);
    //choose_e(r, e_val, n_val, phi_n); for random e
    e_val = 65537; //preferred e
    calc_inverse(d_val, phi_n, e_val);

    std::pair<mpz_class, mpz_class> pub_key(e_val, n_val);
    std::pair<mpz_class, mpz_class> priv_key(d_val, n_val);


    write_file(key_file_name + ".public_key", pub_key);
    write_file(key_file_name + ".private_key", priv_key);


}

//Import key from file
std::pair<mpz_class, mpz_class> construct_key_from_file(std::string &key_str) {
    size_t pos = key_str.find("\n");
    mpz_class first(key_str.substr(0, pos));
    mpz_class second(key_str.substr(pos + 1, std::string::npos));
    return std::pair<mpz_class, mpz_class>(first, second);
}

//encrypt a file with a key
void encrypt_file(std::string &file_to_encrypt, std::string &public_key) {


    std::string file_contents = read_file(file_to_encrypt);
    std::string public_key_str = read_file(public_key);

    if (file_contents.empty() || public_key_str.empty())
        return;

    std::pair<mpz_class, mpz_class> pub_key = construct_key_from_file(public_key_str);


    int n_blocks = ((file_contents.length() * sizeof(char)) + MAX_DATA_SIZE - 1) / MAX_DATA_SIZE;


    char *plaintext_p = (char *) malloc(BLOCK_SIZE * n_blocks);
    memset(plaintext_p, 0, BLOCK_SIZE * n_blocks);
    char *ciphertext = (char *) malloc(BLOCK_SIZE * n_blocks);
    memset(ciphertext, 0, BLOCK_SIZE * n_blocks);
    char *pl = (char *) malloc(MAX_DATA_SIZE * n_blocks + 1);
    memset(pl, 0, MAX_DATA_SIZE * n_blocks + 1);

    split_blocks(plaintext_p, file_contents, n_blocks);
    encrypt_blocks(ciphertext, plaintext_p, n_blocks, pub_key);

    mpz_t ciphertext_as_number;
    mpz_init(ciphertext_as_number);

    mpz_import(ciphertext_as_number, 1, -1, BLOCK_SIZE * n_blocks, -1, 0, ciphertext);

    char *tmp = mpz_get_str(NULL, 10, ciphertext_as_number);
    std::string str_to_write = tmp;
    write_ciphertext_to_file(file_to_encrypt + ".encrypted", tmp);
    std::cout << "Ciphertext written to file: " << file_to_encrypt + ".encrypted" << std::endl;


    free(plaintext_p);
    free(ciphertext);
    free(pl);


}

//decrypt a file with a key
std::string decrypt_file(std::string &file_to_decrypt, std::string &private_key) {
    std::string file_contents = read_file(file_to_decrypt);
    std::string private_key_str = read_file(private_key);
    if (file_contents.empty() || private_key_str.empty())
        return "";

    std::pair<mpz_class, mpz_class> priv_key = construct_key_from_file(private_key_str);

    int n_blocks = ((file_contents.length() * sizeof(char)) + MAX_DATA_SIZE - 1) / MAX_DATA_SIZE;

    char *ciphertext = (char *) malloc(BLOCK_SIZE * n_blocks);
    memset(ciphertext, 0, BLOCK_SIZE * n_blocks);
    char *decrypted_text = (char *) malloc(BLOCK_SIZE * n_blocks);
    memset(decrypted_text, 0, BLOCK_SIZE * n_blocks);
    char *pl = (char *) malloc(MAX_DATA_SIZE * n_blocks + 1);
    memset(pl, 0, MAX_DATA_SIZE * n_blocks + 1);


    mpz_t ciphertext_as_number;
    mpz_init(ciphertext_as_number);
    mpz_set_str(ciphertext_as_number, file_contents.c_str(), 10);
    mpz_export(ciphertext, NULL, -1, BLOCK_SIZE * n_blocks, -1, 0, ciphertext_as_number);


    decrypt_blocks(decrypted_text, ciphertext, n_blocks, priv_key);
    join_blocks(decrypted_text, pl, n_blocks);


    std::string plaintext(pl);
    free(ciphertext);
    free(decrypted_text);
    free(pl);

    return plaintext;
}

void invalid_usage() {
    std::cout << "Usage:" << std::endl;
    std::cout << "To generate public key: RSA -g key_name" << std::endl;
    std::cout << "for encryption: RSA -e file_to_encrypt public_key_file" << std::endl;
    std::cout << "for decryption: RSA -d file_to_encrypt private_key_file" << std::endl;
}

int main(int argc, char *argv[]) {

    if (argc < 3) {
        invalid_usage();
        return 0;
    }


    std::string flag(argv[1]);
    std::string file_to_encrypt_decrypt(argv[2]);


    if (flag == "-g" && argc == 3) {
        gen_key(file_to_encrypt_decrypt);
    } else if (flag == "-e" && argc == 4) {
        std::string key_file(argv[3]);
        encrypt_file(file_to_encrypt_decrypt, key_file);
    } else if (flag == "-d" && argc == 4) {
        std::string key_file(argv[3]);
        std::string decrypted_string = decrypt_file(file_to_encrypt_decrypt, key_file);
        if (!decrypted_string.empty())
            std::cout << "Decrypted string: " << decrypted_string << std::endl;
    } else {
        invalid_usage();
    }


    return 0;
}
