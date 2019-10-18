#include <iostream>
#include <string>
#include "mime_io.hpp"

int main(int argc,
         char** argv)
{
    constexpr int required_arg_c{7};
    if (argc not_eq required_arg_c)
    {
        std::cerr << "Provide 4 files and 2 int:" << std::endl;
        std::cerr << "./sam2file reffile samfile1 samfile2 outfile dimension qualityTrehsold" << std::endl;
        return -1;
    }

    const std::string ref_file{argv[1]};
    const std::string sam_file_A{argv[2]};
    const std::string sam_file_B{argv[3]};
    const std::string out_file{argv[4]};
    const int dimension = std::stoi(argv[5]);
    const int qualityTrehsold = std::stoi(argv[6]);

    std::cout << "Reffile " << ref_file << std::endl;
    std::cout << "Sam1 " << sam_file_A << std::endl;
    std::cout << "Sam2 " << sam_file_B << std::endl;
    std::cout << "outfile " << out_file << std::endl;

    mime_io::analyse_positions(ref_file,
                               sam_file_A,
                               sam_file_B,
                               out_file,
                               dimension,
                               qualityTrehsold);
}
