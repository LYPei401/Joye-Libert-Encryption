
#include <gmp.h>
#include <gmpxx.h>
#include <chrono>
using namespace std;
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cstring>

unsigned long long get_seed(unsigned int i)
{
    auto now = std::chrono::high_resolution_clock::now();
    auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count();
    return millis ^ i;
}

// 这个函数能得到r值
void find_r(mpz_t N, mpz_t r)
{
    gmp_randstate_t state2;
    gmp_randinit_default(state2);
    unsigned long long seed = get_seed(1);
    gmp_randseed_ui(state2, seed);
    mpz_urandomm(r, state2, N);
    gmp_randclear(state2);
}

// 这个函数能得到x值
void find_x(mpz_t N, mpz_t x, int k)
{
    mpz_t a;
    mpz_t two_k_power;
    mpz_init(a);
    mpz_init(two_k_power);
    gmp_randstate_t state2;
    gmp_randinit_default(state2);
    int seed_i = 0;
    while (true)
    {
        seed_i++;
        unsigned long long seed = get_seed(seed_i) ^ (seed_i);
        gmp_randseed_ui(state2, seed);
        mpz_urandomm(a, state2, N);
        mpz_ui_pow_ui(two_k_power, 2, k);
        mpz_powm(x, a, two_k_power, N);
        // gmp_printf("a = %Zd\n", a);
        if (mpz_cmp(x, N) < 0)
        {
            break;
        }
    }

    mpz_clear(a);
    mpz_clear(two_k_power);
    gmp_randclear(state2);
}

// 加密，给出 N x y k m的值计算出 c的值
void calculate_c(const mpz_t y, const mpz_t m, const mpz_t h, mpz_t N, mpz_t c)
{
    mpz_t y_pow_m;
    mpz_t h_pow_r;
    mpz_t r;

    // 初始化临时变量
    mpz_init(y_pow_m);
    mpz_init(h_pow_r);
    mpz_init(r);
    // 计算
    find_r(N, r);
    mpz_powm(y_pow_m, y, m, N);
    mpz_powm(h_pow_r, h, r, N);

    // 计算 c = (y^m) * (x^(2^k)) mod N
    mpz_mul(c, y_pow_m, h_pow_r);

    mpz_mod(c, c, N);

    // 释放临时变量
    mpz_clear(y_pow_m);
    mpz_clear(h_pow_r);
    mpz_clear(r);
}

void read_pq(mpz_t N, mpz_t y, mpz_t h, int *k)
{

    FILE *input_file = fopen("pq_values.txt", "r");
    char line[1024]; // Assuming no line exceeds 1024 characters; adjust if necessary

    if (input_file)
    {
        while (fgets(line, sizeof(line), input_file))
        {
            if (strncmp(line, "N=", 2) == 0)
            {
                mpz_set_str(N, line + 2, 10); // Convert the string (after "N=") to mpz_t
            }
            else if (strncmp(line, "y=", 2) == 0)
            {
                mpz_set_str(y, line + 2, 10); // Convert the string (after "y=") to mpz_t
            }
            else if (strncmp(line, "k=", 2) == 0)
            {
                *k = atoi(line + 2); // Convert the string (after "k=") to integer
            }
            else if (strncmp(line, "h=", 2) == 0)
            {
                mpz_set_str(h, line + 2, 10); // Convert the string (after "h=") to integer
            }
        }
        fclose(input_file);
    }
}

int main()
{
    // 记录开始时间
    auto start = std::chrono::high_resolution_clock::now();

    mpz_t N, y, c, m, h;
    mpz_init(N);
    mpz_init(y);
    mpz_init(c);
    mpz_init(m);
    mpz_init(h);
    int k;
    read_pq(N, y, h, &k);
    gmp_printf("N = %Zd\n", N);
    gmp_printf("y = %Zd\n", y);
    gmp_printf("h = %Zd\n", h);
    printf("k = %d\n", k);
    // Encrypt
    std::ifstream inFile("m.txt");
    std::ofstream outFile("C.txt");

    if (!inFile)
    {
        std::cerr << "Unable to open m.txt for reading." << std::endl;
        return 1;
    }

    std::string numberStr;
    while (std::getline(inFile, numberStr))
    {
        mpz_set_str(m, numberStr.c_str(), 10);
        calculate_c(y, m, h, N, c);
        outFile << mpz_get_str(NULL, 10, c) << std::endl; // 使用get_str()方法将mpz_class转换为字符串
    }

    inFile.close();
    outFile.close();

    // Cleanup
    mpz_clear(h);
    mpz_clear(y);
    mpz_clear(N);
    mpz_clear(c);
    mpz_clear(m);
    std::cout << "finish Encryption!" << std::endl;
    // 记录结束时间
    auto stop = std::chrono::high_resolution_clock::now();

    // 计算差值
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;

    return 0;
}
