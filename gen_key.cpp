
#include <iostream>
#include <gmp.h>
#include <chrono>
using namespace std;

#include <iostream>
#include <gmp.h>
#include <chrono>
using namespace std;

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

    int count = 0;
    do
    {
        seed_i = (seed_i + 1);
        // unsigned long seed = time(NULL) ^ (seed_i + 1);
        unsigned long long seed = get_seed(seed_i) ^ (seed_i);
        generate_P_pri(P_pri, bit - k, seed);
        // 这里产生的P_pri的比特应该等于 大bit-k
        mpz_mul(p, P_pri, two_power_k);
        count++;
        if (mpz_tstbit(two_power_k, 0) == 0)
        {
            mpz_add_ui(p, p, 1);
        }

        is_prime = mpz_probab_prime_p(p, 50);
    } while (is_prime == 0);

    mpz_clear(two_power_k);
    mpz_clear(P_pri);

    printf("count of genP = %d\n", count);
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
    /*
    findPrime_P(k, bit, p);
    gmp_printf("p = %Zd\n", p);
    printf("generate q: ");
    printf("bit = ? ");
    scanf("%d", &bit);
    findPrime_Q(1, bit, q);
    findPrime_P(1, bit, q);
    gmp_printf("q = %Zd\n", q);
    */
    // 时间
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 10; ++i)
    {
        findPrime_P(k, bit, p);
    }
    //  记录结束时间
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Time of genP P: " << duration.count() << " microseconds" << std::endl;

    printf("generate q: ");

    // 时间
    auto start1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 10; ++i)
    {
        findPrime_Q(1, bit, q);
    }

    // 记录结束时间
    auto stop1 = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1);
    std::cout << "Time of genQ Q: " << duration1.count() << " microseconds" << std::endl;

    // 时间
    auto start2 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 10; ++i)
    {
        findPrime_P(1, bit, q);
    }

    // 记录结束时间
    auto stop2 = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);
    std::cout << "Time of genP Q: " << duration2.count() << " microseconds" << std::endl;

    mpz_mul(N, p, q);

    find_y(p, q, N, y);
    find_h(p, q, N, h, k);
    // gmp_printf("N = %Zd\n", N);
    // gmp_printf("y = %Zd\n", y);
    // gmp_printf("h = %Zd\n", h);

    // Save p and q to a txt file
    FILE *output_file = fopen("pq_values.txt", "w");
    fprintf(output_file, "k=%d\n", k);
    gmp_fprintf(output_file, "p=%Zd\n", p);
    gmp_fprintf(output_file, "q=%Zd\n", q);
    gmp_fprintf(output_file, "N=%Zd\n", N);
    gmp_fprintf(output_file, "y=%Zd\n", y);
    gmp_fprintf(output_file, "h=%Zd\n", h);

    fclose(output_file);

    // Cleanup
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(N);
    mpz_clear(y);
    mpz_clear(h);

    return 0;
}
