using namespace std;
#include <gmp.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <chrono>

#include <vector>
#include <list>
#include <cstdio>
#include <algorithm>

void read_pq(mpz_t p, mpz_t y, int *k)
{
    FILE *input_file = fopen("F_pq_values.txt", "r");
    char line[1024];

    if (input_file)
    {
        while (fgets(line, sizeof(line), input_file))
        {
            if (strncmp(line, "p=", 2) == 0)
            {
                mpz_set_str(p, line + 2, 10);
            }
            else if (strncmp(line, "y=", 2) == 0)
            {
                mpz_set_str(y, line + 2, 10);
            }
            else if (strncmp(line, "k=", 2) == 0)
            {
                *k = atoi(line + 2);
            }
        }
        fclose(input_file);
    }
}

// (a/p)^(i)_(2^two_pow)
void n_th_pow_res_towpower(mpz_t a, mpz_t p, int two_pow, int i, mpz_t mod)
{

    mpz_t n_copy, a_copy, p_copy;
    mpz_init(n_copy);
    mpz_init(a_copy);
    mpz_init(p_copy);

    // 拷贝，以防修改原始值
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

std::vector<int> de(std::vector<std::vector<string>> T, mpz_t C, mpz_t p, int l, int n)
{
    mpz_t z, inv, numStr;
    std::vector<int> M(n);
    mpz_init(z);
    mpz_init(inv);
    mpz_init(numStr);
    n_th_pow_res_towpower(C, p, l, 1, z);
    char *str = mpz_get_str(NULL, 10, z);
    auto it = std::find(T[0].begin(), T[0].end(), str);
    if (it != T[0].end())
    {
        int m = std::distance(T[0].begin(), it);
        M[0] = m;
    }
    else
    {
        // std::cout << "the alg is not suitable for step one. \n";
        return M;
    }

    for (int i = 1; i < n; i++)
    {
        n_th_pow_res_towpower(C, p, (i + 1) * l, 1, z);
        mpz_set_ui(inv, 1);
        for (int j = 0; j < i; j++)
        {
            std::string Str = T[i - j][M[j]];
            mpz_set_str(numStr, Str.c_str(), 10);
            mpz_mul(inv, inv, numStr);
        }
        mpz_invert(inv, inv, p);
        mpz_mul(z, z, inv);
        mpz_mod(z, z, p);
        char *str = mpz_get_str(NULL, 10, z);
        auto it = std::find(T[0].begin(), T[0].end(), str);
        if (it != T[0].end())
        {
            int m = std::distance(T[0].begin(), it);
            M[i] = m;
        }
        else
        {
            // std::cout << "the alg is not suitable. \n";
            return M;
        }
    }

    free(str);
    mpz_clear(z);
    mpz_clear(inv);
    mpz_clear(numStr);
    return M;
}

void out(mpz_t m_out, std::vector<int> M, int l, int n)
{

    mpz_t m, two_pow;
    mpz_init_set_ui(m, 0);
    mpz_init(two_pow);

    for (int i = 0; i < n; ++i)
    {
        mpz_ui_pow_ui(two_pow, 2, i * l);
        mpz_mul_ui(two_pow, two_pow, M[i]); //  M[i] * 2^(i*l)
        mpz_add(m, m, two_pow);             // 累加到 m
    }

    // gmp_printf("The sum m is: %Zd\n", m);
    mpz_set(m_out, m);
    mpz_clear(m);
    mpz_clear(two_pow);
}

std::vector<std::vector<std::string>> read_T()
{
    std::vector<std::vector<std::string>> all_vectors;
    FILE *input_file = fopen("F_pq_values.txt", "r");
    char line[1024];
    std::vector<std::string> current_vector;
    bool read_next_line = false;
    if (input_file)
    {
        while (fgets(line, sizeof(line), input_file))
        {
            size_t len = strlen(line);
            if (len > 0 && line[len - 1] == '\n')
            {
                line[len - 1] = '\0';
            }

            if (strncmp(line, "T=", 2) == 0)
            {
                if (!current_vector.empty())
                {
                    all_vectors.push_back(current_vector);
                    current_vector.clear();
                }
                read_next_line = true;
            }
            else if (read_next_line && strlen(line) > 0)
            {
                current_vector.push_back(std::string(line));
            }
        }
        if (!current_vector.empty())
        {
            all_vectors.push_back(current_vector);
        }
        fclose(input_file);
    }
    return all_vectors;
}

int main()
{
    // 记录开始时间
    auto start = std::chrono::high_resolution_clock::now();

    mpz_t p, y, c, m_out;
    mpz_init(p);
    mpz_init(y);
    mpz_init(c);
    mpz_init(m_out);
    int k;
    read_pq(p, y, &k);
    // gmp_printf("y = %Zd\n", y);
    // gmp_printf("p = %Zd\n", p);
    printf("k = %d\n", k);

    std::vector<std::vector<string>> T = read_T();

    int n = T.size();
    int l = k / n;
    printf("n = %d\n", n);
    printf("l = %d\n", l);

    // Decrypt
    std::ifstream inFile("C.txt");
    std::ofstream outFile("m_out.txt");

    if (!inFile)
    {
        std::cerr << "Unable to open m.txt for reading." << std::endl;
        return 1;
    }

    std::string numberStr;
    while (std::getline(inFile, numberStr))
    {
        mpz_set_str(c, numberStr.c_str(), 10);
        auto M = de(T, c, p, l, n);
        out(m_out, M, l, n);
        outFile << mpz_get_str(NULL, 10, m_out) << std::endl; // 使用get_str()方法将mpz_class转换为字符串
    }
    inFile.close();
    outFile.close();

    std::cout << "finish Decryption!" << std::endl;
    // Cleanup
    mpz_clear(y);
    mpz_clear(c);
    mpz_clear(p);
    mpz_clear(m_out);
    // 记录结束时间
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;

    return 0;
}
