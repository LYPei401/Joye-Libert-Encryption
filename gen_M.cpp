#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>

int main()
{
    std::ofstream outFile("m.txt");

    // 使用当前时间作为随机数种子
    std::srand(static_cast<unsigned int>(std::time(NULL)));
    int k = 50;
    int number = 2 << (k - 1);
    printf("num = %d\n", number);
    for (int i = 0; i < 10000; ++i)
    {
        int randomNumber = std::rand() % number; //
        outFile << randomNumber << std::endl;
    }

    outFile.close();
    std::cout << "Numbers have been written to m.txt" << std::endl;

    return 0;
}