
#include <iostream>
#include <gmp.h>
#include <chrono>
using namespace std;

#include <iostream>
#include <vector>
#include <cmath>
#include <list>

unsigned long long get_seed(unsigned int i)
{
    auto now = std::chrono::high_resolution_clock::now();
    auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count();
    return millis ^ i;
}

int Jacobi(mpz_t a, mpz_t n)
{
    mpz_t n_copy, a_copy;
    mpz_init(n_copy);
    mpz_init(a_copy);

    mpz_set(n_copy, n);
    mpz_set(a_copy, a);

    if (mpz_sgn(n_copy) <= 0 || mpz_even_p(n_copy))
    {
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

        mpz_set(temp, a_copy);
        mpz_set(a_copy, n_copy);
        mpz_set(n_copy, temp);

        mpz_mod_ui(remainder_a, a, 4);
        mpz_mod_ui(remainder_n, n, 4);

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

void find_y(mpz_t p, mpz_t q, mpz_t N, mpz_t y)
{
    mpz_t a;
    mpz_init(a);
    gmp_randstate_t state2;
    gmp_randinit_default(state2);
    int seed_i = 0;
    while (true)
    {
        seed_i++;
        unsigned long long seed = get_seed(seed_i) ^ (seed_i);
        gmp_randseed_ui(state2, seed);
        mpz_urandomm(a, state2, N);
        if (Jacobi(a, p) == -1 && Jacobi(a, q) == -1)
        {
            mpz_set(y, a);
            break;
        }
    }

    mpz_clear(a);
    gmp_randclear(state2);
}

void find_h(mpz_t p, mpz_t q, mpz_t N, mpz_t h, int k)
{
    mpz_t x;
    mpz_t two_pow_k;
    mpz_init(x);
    mpz_init(two_pow_k);

    find_y(p, q, N, x);
    mpz_ui_pow_ui(two_pow_k, 2, k);
    mpz_powm(h, x, two_pow_k, N);
}

void generate_P_pri(mpz_t prime, int bits, unsigned long seed)
{
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, seed);

    mpz_t min_value, max_value;
    mpz_init(min_value);
    mpz_init(max_value);

    // 设置最小值为 2^(bits-1)
    mpz_ui_pow_ui(min_value, 2, bits - 1);

    // 设置最大值为 2^bits - 1
    mpz_ui_pow_ui(max_value, 2, bits);
    mpz_sub_ui(max_value, max_value, 1);

    do
    {
        mpz_urandomm(prime, state, max_value);
        if (mpz_cmp(prime, min_value) < 0)
        {
            mpz_add(prime, prime, min_value);
        }
        mpz_nextprime(prime, prime);
    } while (mpz_sizeinbase(prime, 2) != bits);

    mpz_clear(min_value);
    mpz_clear(max_value);
    gmp_randclear(state);
}

void findPrime_Q(int k, int bit, mpz_t p)
{
    // Calculate 2^k
    mpz_t two_power_k;

    mpz_t P_copy;
    mpz_t quotient, r;
    mpz_init(quotient);
    mpz_init(r);
    mpz_init(P_copy);
    mpz_init(two_power_k);
    mpz_ui_pow_ui(two_power_k, 2, k);

    mpz_init(p);
    int seed_i = 0;
    int is_prime = 0;

    int count = 0;
    do
    {
        seed_i = (seed_i + 1);
        // unsigned long seed = time(NULL) ^ (seed_i + 1);
        unsigned long long seed = get_seed(seed_i) ^ (seed_i);
        generate_P_pri(P_copy, bit, seed);
        mpz_set(p, P_copy);
        mpz_sub_ui(P_copy, P_copy, 1);
        count++;
        mpz_tdiv_qr(quotient, r, P_copy, two_power_k);
        if (mpz_sgn(r) != 0)
        {
            gmp_printf("This way to generate q is not suitable.\n");
        }
        // 这里产生一个bit的素数，然后再检查是否符合条件
        is_prime = mpz_probab_prime_p(quotient, 50);
    } while (is_prime == 0);
    mpz_clear(two_power_k);
    mpz_clear(P_copy);
    mpz_clear(quotient);
    mpz_clear(r);
    printf("count of genQ = %d\n", count);
}

void findPrime_P(int k, int bit, mpz_t p)
{
    // Calculate 2^k
    mpz_t two_power_k;
    mpz_t P_pri;
    mpz_init(P_pri);
    mpz_init(two_power_k);
    mpz_ui_pow_ui(two_power_k, 2, k);
    mpz_init(p);
    int seed_i = 0;
    int is_prime = 0;
    do
    {
        seed_i = (seed_i + 1);
        unsigned long long seed = get_seed(seed_i) ^ (seed_i);
        generate_P_pri(P_pri, bit - k, seed);
        // 这里产生的P_pri的比特应该等于 大bit-k
        mpz_mul(p, P_pri, two_power_k);
        if (mpz_tstbit(two_power_k, 0) == 0)
        {
            mpz_add_ui(p, p, 1);
        }
        is_prime = mpz_probab_prime_p(p, 50);
    } while (is_prime == 0);
    mpz_clear(two_power_k);
    mpz_clear(P_pri);
}

void n_th_pow_res_towpower(mpz_t a, mpz_t p, int two_pow, int i, mpz_t mod)
{

    mpz_t n_copy, a_copy, p_copy;
    mpz_init(n_copy);
    mpz_init(a_copy);
    mpz_init(p_copy);

    mpz_ui_pow_ui(n_copy, 2, two_pow);
    mpz_set(a_copy, a);
    mpz_set(p_copy, p);
    mpz_t exponent, gcd, p_minus_one;
    mpz_init(gcd);
    mpz_init(exponent);
    mpz_init(p_minus_one);

    mpz_sub_ui(p_minus_one, p_copy, 1);
    mpz_gcd(gcd, n_copy, p_minus_one);
    mpz_divexact(exponent, p_minus_one, gcd);
    mpz_mul_ui(exponent, exponent, i);
    mpz_powm(mod, a_copy, exponent, p);
    mpz_clear(a_copy);
    mpz_clear(n_copy);
    mpz_clear(p_copy);
    mpz_clear(exponent);
    mpz_clear(gcd);
    mpz_clear(p_minus_one);
}

std::vector<std::string> find_T_i(int l, int n_i, mpz_t p, mpz_t y)
{
    size_t length = static_cast<size_t>(std::pow(2, l));
    std::vector<std::string> vec(length);

    for (int i = 0; i < length; ++i)
    {
        mpz_t v;
        mpz_init(v);
        n_th_pow_res_towpower(y, p, l * n_i, i, v); // 假设该函数适用于你的需求
        char *str = mpz_get_str(NULL, 10, v);
        vec[i] = str;
        free(str);
        mpz_clear(v);
    }

    return vec;
}
std::vector<std::vector<std::string>> find_T(int l, int n, mpz_t p, mpz_t y)
{
    std::vector<std::vector<std::string>> all_vectors;

    for (int i = 1; i < 1 + n; ++i)
    {
        all_vectors.push_back(find_T_i(l, i, p, y));
    }

    return all_vectors;
}

int main()
{

    mpz_t p, q, N, y, h;
    mpz_init(p);
    mpz_init(q);
    mpz_init(N);
    mpz_init(y);
    mpz_init(h);
    int k, bit;

    //  Generate p
    printf("generate p: ");
    printf("k = ? ");
    scanf("%d", &k);
    printf("bit = ? ");
    scanf("%d", &bit);
    findPrime_P(k, bit, p);
    gmp_printf("p = %Zd\n", p);

    // Generate q
    printf("generate q: ");
    printf("bit = ? ");
    scanf("%d", &bit);
    findPrime_Q(1, bit, q);
    gmp_printf("q = %Zd\n", q);
    mpz_mul(N, p, q);

    find_y(p, q, N, y);
    find_h(p, q, N, h, k);
    // gmp_printf("N = %Zd\n", N);
    // gmp_printf("y = %Zd\n", y);
    // gmp_printf("h = %Zd\n", h);

    FILE *output_file = fopen("F_pq_values.txt", "w");
    fprintf(output_file, "k=%d\n", k);
    gmp_fprintf(output_file, "p=%Zd\n", p);
    gmp_fprintf(output_file, "q=%Zd\n", q);
    gmp_fprintf(output_file, "N=%Zd\n", N);
    gmp_fprintf(output_file, "y=%Zd\n", y);
    gmp_fprintf(output_file, "h=%Zd\n", h);

    int l;
    printf("l = ? ");
    scanf("%d", &l);
    int n = k / l;
    auto T = find_T(l, n, p, y);
    int vec_count = 1;

    for (const auto &vec : T)
    {
        fprintf(output_file, "T=%d\n", vec_count);
        for (const auto &val : vec)
        {
            fprintf(output_file, "%s\n", val.c_str()); // 每个元素后换行
        }
        vec_count++;
    }

    fclose(output_file);
    // Cleanup
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(N);
    mpz_clear(y);
    mpz_clear(h);
    printf("finish!");

    return 0;
}
