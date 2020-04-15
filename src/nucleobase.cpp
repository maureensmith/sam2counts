#include "nucleobase.hpp"

#include <stdexcept>
#include <utility>


namespace nucleotid
{

const std::array<char, 6> nucleobase::bases{'A', 'C', 'G', 'T', 'N','X'};

nucleobase::nucleobase(char base)
{
    //bit shift = divide by 2
    //base = (base - 'A') >> 1;
    //const auto id = (base - 'A') / 2;
//    if (std::max(std::min(10, id), 0) != id)
    {
 //       throw std::invalid_argument(std::string{base} + std::string{" is not a valid base"});
    }
    //TODO: sicherer is switchen
    //mapping of the nucleotide bases to IDs (eg A-A/2 = 0, C-A/2 = 1 etc)
    //A=0, C=1, G=2, T=3, N=4, X=5
    switch(base)
    {
        case 'A': this->base = 0; break;
        case 'C': this->base = 1; break;
        case 'G': this->base = 2; break;
        case 'T': this->base = 3; break;
        case 'N': this->base = 4; break;
        default: this->base = 5; break;
    }

//    static uint8_t bla[] = {0, 1, 5, 2, 5, 5, 4, 5, 5, 3, 3};
//    this->base = bla[static_cast<unsigned> (base)];
}

nucleobase::nucleobase(const int base_id) 
{
    //TODO changed to 5, because N is considered here
    if (std::max(std::min(5, base_id), 1) != base_id)
    {
        throw std::invalid_argument(std::to_string(base) + std::string{" is not a valid base id"});
    }

    base = base_id - 1;
}

}

