using namespace std;
#include <gmp.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <chrono>

void Decryption(mpz_t m, const mpz_t c, const mpz_t y, const mpz_t p, int k)
{
    mpz_t B, C, z, D, p_minus_1, two_power, q_D, two_power_k_minus_j;
    mpz_init(D);
    mpz_init(q_D);

    mpz_init(p_minus_1);
    mpz_init(two_power);
    mpz_init(two_power_k_minus_j);
    mpz_init_set_ui(m, 0);
    mpz_init_set_ui(B, 1);
    mpz_init(C);

    mpz_sub_ui(p_minus_1, p, 1);
    mpz_ui_pow_ui(two_power, 2, k);
    mpz_cdiv_q(q_D, p_minus_1, two_power);

    mpz_powm(C, c, q_D, p);
    mpz_powm(D, y, q_D, p);

    mpz_invert(D, D, p);

    mpz_t Y;
    mpz_t Y_exp;
    mpz_init(Y);
    mpz_init(Y_exp);
    mpz_mul_ui(Y_exp, q_D, 24);
    mpz_powm(Y, y, Y_exp, p);
    for (int j = 1; j <= k - 1; j++)
    {
        mpz_init(z);
        mpz_ui_pow_ui(two_power_k_minus_j, 2, k - j);
        mpz_powm(z, C, two_power_k_minus_j, p);
        if (mpz_cmp_ui(z, 1) != 0)
        {
            mpz_add(m, m, B);
            mpz_mul(C, C, D);
            mpz_mod(C, C, p);
        }

        mpz_mul_ui(B, B, 2);
        mpz_mul(D, D, D);
        mpz_mod(D, D, p);
    }

    if (mpz_cmp_ui(C, 1) != 0)
    {
        mpz_add(m, m, B);
    }
    mpz_clear(z);
    mpz_clear(p_minus_1);
    mpz_clear(q_D);
    mpz_clear(two_power);
    mpz_clear(two_power_k_minus_j);
    mpz_clear(B);
    mpz_clear(D);
    mpz_clear(C);
}

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
        Decryption(m_out, c, y, p, k);
        outFile << mpz_get_str(NULL, 10, m_out) << std::endl;
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
