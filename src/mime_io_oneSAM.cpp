#include "mime_io.hpp"
#include "io_tools.hpp"
#include "count.hpp"
#include "reference.hpp"
#include "aligner.hpp"


#include <algorithm>
#include <iterator>
#include <utility>
#include <iostream>
#include <chrono>


namespace mime_io
{

namespace
{
void workflow_1(const utils::reference& reference,
                std::ifstream& input_a,
                const std::string& out_file)
{
    std::string line_a = io_tools::get_first_data_line(input_a);
    //the paired reads are consecutively
    std::string line_b = io_tools::get_first_data_line(input_a);

    utils::counter_1 counter{reference};
    utils::ref_map read;

    utils::aligner aligner{reference, read};

    unsigned prep = 0;
    unsigned align_a = 0;
    unsigned align_b = 0;
    unsigned count = 0;

    //unsigned line_no = 0;

    while (input_a.good())
    {
        read.clear();
        //missing check whether the lines belong together
        auto now = std::chrono::high_resolution_clock::now();
        const auto is_prepared = aligner.prepare(line_a, line_b);
        auto diff  = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - now);
        prep += diff.count();
        if (is_prepared)
        {
            now = std::chrono::high_resolution_clock::now();
            aligner.align();
            diff  = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - now);
            align_a += diff.count();

            now = std::chrono::high_resolution_clock::now();
            aligner.align_1(counter);
            //utils::align_to_reference_2(sep_b.read_seq_it, sep_b.pos_in_ref, sep_b.cigar_it, sep_b.cigar_end, read, reference);
            diff  = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - now);
            align_b += diff.count();

            now = std::chrono::high_resolution_clock::now();
            counter.count(read);
            diff  = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - now);
            count += diff.count();
        }

        /*
           ++line_no;
           if (line_no == 10000)
           {
           break;
           }
           */

        std::getline(input_a, line_a);
        std::getline(input_b, line_b);
    }

    counter.write_to_file(out_file);

    std::cout << "Prepare: " << prep << std::endl;
    std::cout << "Align A: " << align_a << std::endl;
    std::cout << "Align B: " << align_b << std::endl;
    std::cout << "  Count: " << count << std::endl;
}

void workflow_2(const utils::reference& reference,
                std::ifstream& input_a,
                std::ifstream& input_b,
                const std::string& out_file)
{
    std::string line_a = io_tools::get_first_data_line(input_a);
    std::string line_b = io_tools::get_first_data_line(input_b);

    utils::counter_2 counter{reference};
    utils::ref_map read;

    utils::aligner aligner{reference, read};

    unsigned align_a = 0;
    unsigned align_b = 0;
    unsigned count = 0;

    while (input_a.good() and input_b.good())
    {
        read.clear();
        //missing check whether the lines belong together
        if (aligner.prepare(line_a, line_b))
        {
            auto now = std::chrono::high_resolution_clock::now();
            aligner.align();
            auto diff  = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - now);
            align_a += diff.count();

            now = std::chrono::high_resolution_clock::now();
            aligner.align_2();
            diff  = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - now);
            align_b += diff.count();

            std::sort(read.begin(), read.end(), [](const auto& lhs, const auto& rhs)
                      {
                      return lhs.first < rhs.first;
                      });

            now = std::chrono::high_resolution_clock::now();
            counter.count(read);
            diff  = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - now);
            count += diff.count();
            }

        std::getline(input_a, line_a);
        std::getline(input_b, line_b);
    }

    counter.write_to_file(out_file);

    std::cout << "Align A: " << align_a << std::endl;
    std::cout << "Align B: " << align_b << std::endl;
    std::cout << "  Count: " << count << std::endl;
}
}

void analyse_positions(const std::string& ref,
                       const std::string& sam_a,
                       const std::string& out_file,
                       const int dimension)
{
    const auto ref_map = io_tools::read_reference(ref);

    std::ifstream input_a(sam_a);

    switch (dimension)
    {
      case 1:
        workflow_1(ref_map, input_a, out_file);
        break;
      case 2:
        workflow_2(ref_map, input_a, out_file);
        break;
      default:
        throw "";
    }
}

}
