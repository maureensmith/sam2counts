#ifndef COUNT_HPP_
#define COUNT_HPP_

#include "reference.hpp"
#include "ref_map.hpp"
#include "utils.hpp"

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
    : data(ref.size())
    {
    }

    void count(const ref::ref_map& read);

    void count(const std::size_t index,
               const char base)
    {
        ++data[index][nucleotid::nucleobase{base}.to_id()];
    }

    void write_to_file(const std::string& out_file);

    static constexpr unsigned nucleobase_count{4u};

  private:
    //containing counts for all 4 nucleotides for each position of the reference (counted from 0, for input/output -1)
    std::vector<std::array<count_type, nucleobase_count>> data;
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

    void write_to_file(const std::string& out_file);

  private:
    const std::size_t size;
    std::vector<std::array<count_type, counter_1::nucleobase_count * counter_1::nucleobase_count>> data;
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

    void write_to_file(const std::string& out_file);

  private:
    const std::size_t size;
    std::vector<std::array<count_type, counter_1::nucleobase_count * counter_1::nucleobase_count * counter_1::nucleobase_count>> data;
};
}

#endif
