#include "mime_io_paired.hpp"
#include "io_tools.hpp"
#include "count.hpp"
#include "reference.hpp"
#include "aligner.hpp"


#include <algorithm>
#include <iostream>
#include <chrono>


namespace mime_io_paired
{

namespace
{
void workflow_1(const ref::reference& reference,
                std::ifstream& input_a,
                std::ifstream& input_b,
                const std::string& out_file,
                const int qualityThreshold)
{
    std::string line_a = io_tools::get_first_data_line(input_a);
    std::string line_b = io_tools::get_first_data_line(input_b);

    count::counter_1 counter{reference};
    ref::ref_map read;

    aligner::aligner aligner{reference, read, qualityThreshold};

    // record to time for each step
    unsigned prep = 0;
    unsigned align_a = 0;
    unsigned align_b = 0;
    unsigned count = 0;

    while (input_a.good() and input_b.good())
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
            //aligner.align_2();
            //utils::align_to_reference_2(sep_b.read_seq_it, sep_b.pos_in_ref, sep_b.cigar_it, sep_b.cigar_end, read, reference);
            diff  = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - now);
            align_b += diff.count();

            now = std::chrono::high_resolution_clock::now();
            counter.count(read);
            diff  = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - now);
            count += diff.count();
        }

        std::getline(input_a, line_a);
        std::getline(input_b, line_b);
    }

    counter.write_to_file(out_file);

    std::cout << "Prepare: " << prep << std::endl;
    std::cout << "Align A: " << align_a << std::endl;
    std::cout << "Align B: " << align_b << std::endl;
    std::cout << "  Count: " << count << std::endl;
}

void workflow_2(const ref::reference& reference,
                std::ifstream& input_a,
                std::ifstream& input_b,
                const std::string& out_file,
                const int qualityThreshold)
{
    std::string line_a = io_tools::get_first_data_line(input_a);
    std::string line_b = io_tools::get_first_data_line(input_b);

    count::counter_2 counter{reference};
    ref::ref_map read;

    aligner::aligner aligner{reference, read, qualityThreshold};

    // save times/duration for the alignments
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

void workflow_3(const ref::reference& reference,
                std::ifstream& input_a,
                std::ifstream& input_b,
                const std::string& out_file,
                const int qualityThreshold)
{
    std::string line_a = io_tools::get_first_data_line(input_a);
    std::string line_b = io_tools::get_first_data_line(input_b);

    count::counter_3 counter{reference};
    ref::ref_map read;

    aligner::aligner aligner{reference, read, qualityThreshold};

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
                       const std::string& sam_b,
                       const std::string& out_file,
                       const int dimension,
                       const int qualityThreshold)
{
    const auto ref_map = io_tools::read_reference(ref);

    std::ifstream input_a(sam_a);
    std::ifstream input_b(sam_b);

    switch (dimension)
    {
      case 1:
        workflow_1(ref_map, input_a, input_b, out_file, qualityThreshold);
        break;
      case 2:
        workflow_2(ref_map, input_a, input_b, out_file, qualityThreshold);
        break;
      case 3:
        workflow_3(ref_map, input_a, input_b, out_file, qualityThreshold);
        break;
      default:
        throw "";
    }
}
}

