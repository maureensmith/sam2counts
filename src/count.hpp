#ifndef COUNT_HPP_
#define COUNT_HPP_

#include "reference.hpp"
#include "ref_map.hpp"
#include "utils.hpp"
#include "nucleobase.hpp"

#include <array>
#include <string>
#include <vector>
#include <utility>

namespace count
{
class counter_1 final
{
  public:
    using count_type = unsigned;

    // Delete default constructor
    counter_1() = delete;
    ~counter_1() = default;


    counter_1(const ref::reference& ref)
    : nucleobase_count(nucleotide::numberOfValidSymbols(false)), data(ref.size(), std::vector<count_type>(nucleobase_count))
    {
    }

    counter_1(const ref::reference& ref, bool a)
            : ambig{a}, nucleobase_count(nucleotide::numberOfValidSymbols(a)), data(ref.size(), std::vector<count_type>(nucleobase_count))
    {
    }

    void count(const ref::ref_map& read);

    void count(const ref::ref_map& read, const unsigned times);

    void count(const std::size_t index,
               const char base)
    {
        ++data[index][nucleotide::nucleobase{base}.to_id()];
    }

    void write_to_file(const std::string& out_file);


  private:
    //containing counts for all 4 nucleotides for each position of the reference (counted from 0, for input/output -1)
    //std::vector<std::array<count_type, nucleobase_count>> data;
    const bool ambig{false};
    const unsigned nucleobase_count;
    std::vector<std::vector<count_type>> data;
};


class counter_2 final
{
  public:
    using count_type = unsigned;

    counter_2() = delete;
    ~counter_2() = default;

    counter_2(const ref::reference& ref)
    : size{ref.size()},  data(((size * size) - size) / 2, {0})
    {
    }

    void count(const ref::ref_map& read);
    void count(const ref::ref_map& read, const unsigned times);
    void count(const std::size_t pos1, const std::size_t pos2,const char base1, const char base2, const int times)
    {
        const auto pos1_idx = (pos1 - 1);
        const auto i = pos1_idx * (size-1) - utils::choose(pos1_idx, 2);
        std::size_t pairwise_pos_index = i + pos2 - pos1_idx - 2;
        std::size_t pairwise_nucl_index = nucleobase_count * nucleotide::nucleobase{base1}.to_id() + nucleotide::nucleobase{base2}.to_id();
        data[pairwise_pos_index][pairwise_nucl_index] += times;
    }

    void write_to_file(const std::string& out_file);

  private:
    static constexpr unsigned nucleobase_count{nucleotide::numberOfBasicSymbols};
    const std::size_t size;
    std::vector<std::array<count_type, nucleobase_count * nucleobase_count>> data;
};

class counter_3 final
{
  public:
    using count_type = unsigned;

    counter_3() = delete;
    ~counter_3() = default;

    counter_3(const ref::reference& ref)
        : size{ref.size()},  data(utils::choose(size,3), {0})
    {
    }

    void count(const ref::ref_map& read);

    void count(const ref::ref_map& read, const unsigned times);

    void write_to_file(const std::string& out_file);

  private:
    const std::size_t size;
    static constexpr unsigned nucleobase_count{nucleotide::numberOfBasicSymbols};
    std::vector<std::array<count_type, nucleobase_count * nucleobase_count * nucleobase_count>> data;
};
}

#endif
