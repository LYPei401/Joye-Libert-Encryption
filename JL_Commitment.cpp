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

// 这个函数能从range中取值
void find_pow_range(const mpz_t n, mpz_t a, int k)
{
    gmp_randstate_t state2;
    gmp_randinit_default(state2);
    unsigned long long seed = get_seed(1);
    gmp_randseed_ui(state2, seed);
    mpz_t two_pow_k;
    mpz_init(two_pow_k);
    mpz_ui_pow_ui(two_pow_k, 2, k);
    mpz_mul(two_pow_k, n, two_pow_k);

    mpz_urandomm(a, state2, two_pow_k);
    gmp_randclear(state2);
    mpz_clear(two_pow_k);
}

// c = y^(m*2^k)*h^(r*2^k) mod N  calculate_c(y, m, h, r, N, c, k);
void calculate_c(const mpz_t y, const mpz_t m, const mpz_t h, const mpz_t r, const mpz_t N, mpz_t c, int k)
{
    mpz_t y_pow_m;
    mpz_t h_pow_r;
    mpz_t two_pow_k;
    mpz_t m_power;
    mpz_t r_power;
    // 初始化临时变量
    mpz_init(y_pow_m);
    mpz_init(h_pow_r);
    mpz_init(two_pow_k);
    mpz_init(r_power);
    mpz_init(m_power);
    // 计算
    mpz_ui_pow_ui(two_pow_k, 2, k);
    mpz_mul(m_power, m, two_pow_k);
    mpz_mul(r_power, r, two_pow_k);
    mpz_powm(y_pow_m, y, m_power, N);
    mpz_powm(h_pow_r, h, r_power, N);

    // 计算 c = (y^m) * (x^(2^k)) mod N
    mpz_mul(c, y_pow_m, h_pow_r);
    mpz_mod(c, c, N);

    // 释放临时变量
    mpz_clear(y_pow_m);
    mpz_clear(h_pow_r);
    mpz_clear(two_pow_k);
    mpz_clear(m_power);
    mpz_clear(r_power);
}

int Jacobi(mpz_t a, mpz_t n)
{
    mpz_t n_copy, a_copy;
    mpz_init(n_copy);
    mpz_init(a_copy);

    // 拷贝a和n，以防修改原始值
    mpz_set(n_copy, n);
    mpz_set(a_copy, a);

    if (mpz_sgn(n_copy) <= 0 || mpz_even_p(n_copy))
    {
        // Jacobi符号定义仅对n为奇素数时有效
        mpz_clear(n_copy);
        mpz_clear(a_copy);
        return 0;
    }

    if (mpz_sgn(a_copy) < 0)
    {
        // 对Jacobi符号中的负数的时候我们暂不讨论
        return 0;
    }

    int s = 1;
    mpz_t two, r, temp;
    mpz_init(two);
    mpz_init(r);
    mpz_init(temp);

    mpz_t remainder_a, remainder_n;
    mpz_init(remainder_a);
    mpz_init(remainder_n);
    mpz_set_ui(two, 2);

    while (mpz_sgn(a_copy) != 0)
    {
        while (mpz_divisible_p(a_copy, two))
        {
            mpz_divexact(a_copy, a_copy, two);
            mpz_mod_ui(r, n_copy, 8);
            if (mpz_cmp_ui(r, 3) == 0 || mpz_cmp_ui(r, 5) == 0)
            {
                s = -s;
            }
        }

        // 交换a和n
        mpz_set(temp, a_copy);
        mpz_set(a_copy, n_copy);
        mpz_set(n_copy, temp);

        mpz_mod_ui(remainder_a, a, 4);
        mpz_mod_ui(remainder_n, n, 4);

        // 处理符号变换
        mpz_mod_ui(remainder_a, a_copy, 4);
        mpz_mod_ui(remainder_n, n_copy, 4);
        if (mpz_cmp_ui(remainder_a, 3) == 0 && mpz_cmp_ui(remainder_n, 3) == 0)
        {
            s = -s;
        }

        mpz_mod(a_copy, a_copy, n_copy);
    }
    mpz_clear(remainder_a);
    mpz_clear(remainder_n);

    mpz_clear(a_copy);
    mpz_clear(two);
    mpz_clear(r);
    mpz_clear(temp);

    if (mpz_cmp_ui(n_copy, 1) == 0)
    {
        mpz_clear(n_copy);
        return s;
    }
    else
    {
        mpz_clear(n_copy);
        return 0;
    }
}

void read_pq(mpz_t N, mpz_t h, mpz_t y, int *k)
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

    mpz_t N, y, c, m, h, r, B;
    mpz_init(N);
    mpz_init(y);
    mpz_init(c);
    mpz_init(m);
    mpz_init(h);
    mpz_init(r);
    mpz_init(B);
    int k;
    int s = 76;
    int t = 84;
    // prepare
    read_pq(N, y, h, &k);
    gmp_printf("N = %Zd\n", N);
    gmp_printf("y = %Zd\n", y);
    gmp_printf("h = %Zd\n", h);
    printf("k = %d\n", k);
    find_pow_range(N, r, 0);
    // mpz_ui_pow_ui(B, 2, k);
    mpz_set_ui(B, 1024);
    find_pow_range(B, m, 0);
    //mpz_set_ui(m, 0);

    calculate_c(y, m, h, r, N, c, k);

    // Prover send 1
    mpz_t v, w, d;
    mpz_init(v);
    mpz_init(w);
    mpz_init(d);
    find_pow_range(y, v, s);
    //   find_pow_range(B, v, s + t);
    mpz_set_ui(v, 0);

    find_pow_range(N, w, s + t);
    calculate_c(y, v, h, w, N, d, k);

    // V send 1
    mpz_t e;
    mpz_t o;
    mpz_init(e);
    mpz_init(o);
    mpz_set_ui(o, 1);
    find_pow_range(o, e, t);

    // P send 2
    mpz_t z_m, z_r;
    mpz_init(z_m);
    mpz_init(z_r);
    mpz_mul(z_m, e, m);
    /*
    gmp_printf("e = %Zd\n", e);
    gmp_printf("m = %Zd\n", m);
    gmp_printf("v = %Zd\n", v);
    gmp_printf("z_m = %Zd\n", z_m);
    */
    mpz_mul(z_r, e, r);
    mpz_add(z_m, z_m, v);
    mpz_add(z_r, z_r, w);

    // V accept
    int j_c = Jacobi(c, N);
    int j_d = Jacobi(d, N);
    mpz_t accept_left, accept_right;
    mpz_init(accept_left);
    mpz_init(accept_right);
    mpz_powm(accept_left, c, e, N);
    mpz_mul(accept_left, accept_left, d);
    mpz_mod(accept_left, accept_left, N);
    calculate_c(y, z_m, h, z_r, N, accept_right, k);

    mpz_t two_pow_st;
    mpz_init(two_pow_st);
    mpz_ui_pow_ui(two_pow_st, 2, s + t);
    mpz_mul(two_pow_st, two_pow_st, B);

    if (mpz_cmp(z_m, two_pow_st) < 0 || mpz_cmp(z_m, two_pow_st) == 0)
    {
        printf("1, z_m in the range B\n");
    }
    if (j_c == 1)
    {
        printf("2, J(c) = 1\n");
    }
    if (j_d == 1)
    {
        printf("3, J(d) = 1\n");
    }
    if (mpz_cmp(accept_left, accept_right) == 0)
    {
        printf("4, equation mod N stands\n ");
    }

    // clean
    mpz_clear(v);
    mpz_clear(w);
    mpz_clear(d);
    mpz_clear(e);
    mpz_clear(o);
    mpz_clear(z_m);
    mpz_clear(z_r);
    mpz_clear(accept_left);
    mpz_clear(accept_right);

    mpz_clear(two_pow_st);

    // Cleanup
    mpz_clear(h);
    mpz_clear(y);
    mpz_clear(N);
    mpz_clear(c);
    mpz_clear(m);
    mpz_clear(r);
    mpz_clear(B);

    std::cout << "finish!" << std::endl;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;

    return 0;
}