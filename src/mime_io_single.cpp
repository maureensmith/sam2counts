//
// Created by Smith, Maureen on 16.04.20.
//

#include <fstream>
#include <chrono>
#include <iostream>

#include "mime_io_single.hpp"
#include "io_tools.hpp"
#include "count.hpp"
#include "aligner.hpp"

namespace mime_io_single {

    namespace {

        //TODO: counter -> interface, call the same routine for different counter
        void workflow_1(const ref::reference &reference,
                        std::ifstream &input,
                        const std::string &out_file,
                        const int qualityThreshold,
                        // bool checkInsertions,
                        bool checkDeletions,
                        const bool ambig) {
            std::string line = io_tools::get_first_data_line(input);

            count::counter_1 counter{reference, checkDeletions, ambig};
            ref::ref_map read;

            aligner::aligner aligner{reference, read, qualityThreshold, 
            // checkInsertions, 
            checkDeletions, ambig};

            // record to time for each step
            unsigned prep = 0;
            unsigned align_a = 0;
            unsigned count = 0;

            while (input.good()) {
                read.clear();
                auto now = std::chrono::high_resolution_clock::now();
                const auto is_prepared = aligner.prepare(line);
                auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::high_resolution_clock::now() - now);
                prep += diff.count();
                if (is_prepared) {
                    now = std::chrono::high_resolution_clock::now();
                    aligner.align();
                    diff = std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::high_resolution_clock::now() - now);
                    align_a += diff.count();

                    now = std::chrono::high_resolution_clock::now();
                    counter.count(read);
                    diff = std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::high_resolution_clock::now() - now);
                    count += diff.count();
                }
                std::getline(input, line);
            }
            counter.write_to_file(out_file);

            std::cout << "Prepare: " << prep << std::endl;
            std::cout << "Align: " << align_a << std::endl;
            std::cout << "Count: " << count << std::endl;
        }


        void workflow_2(const ref::reference &reference,
                        std::ifstream &input,
                        const std::string &out_file,
                        const int qualityThreshold
                        // bool checkInsertions,bool checkDeletions
                        ) {
            std::string line = io_tools::get_first_data_line(input);

            count::counter_2 counter{reference};
            ref::ref_map read;

            aligner::aligner aligner{reference, read, qualityThreshold, 
            // checkInsertions, checkDeletions
            false //<-checkDeletions
            };

            // record to time for each step
            unsigned prep = 0;
            unsigned align_a = 0;
            unsigned count = 0;

            while (input.good()) {
                read.clear();
                auto now = std::chrono::high_resolution_clock::now();
                const auto is_prepared = aligner.prepare(line);
                auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::high_resolution_clock::now() - now);
                prep += diff.count();

                if (is_prepared) {
                    now = std::chrono::high_resolution_clock::now();
                    aligner.align();
                    diff = std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::high_resolution_clock::now() - now);
                    align_a += diff.count();

                    now = std::chrono::high_resolution_clock::now();
                    counter.count(read);
                    diff = std::chrono::duration_cast<std::chrono::nanoseconds>(
                            std::chrono::high_resolution_clock::now() - now);
                    count += diff.count();
                }
                std::getline(input, line);
            }
            counter.write_to_file(out_file);

            std::cout << "Prepare: " << prep << std::endl;
            std::cout << "Align: " << align_a << std::endl;
            std::cout << "Count: " << count << std::endl;
        }
    }

    void workflow_3(const ref::reference &reference,
                    std::ifstream &input,
                    const std::string &out_file,
                    const int qualityThreshold
                    // bool checkInsertions,bool checkDeletions
                    ) {
        std::string line = io_tools::get_first_data_line(input);

        count::counter_3 counter{reference};
        ref::ref_map read;

        aligner::aligner aligner{reference, read, qualityThreshold, 
        // checkInsertions, checkDeletions
        false //<-checkDeletions
        };

        // record to time for each step
        unsigned prep = 0;
        unsigned align_a = 0;
        unsigned count = 0;

        while (input.good()) {
            read.clear();
            auto now = std::chrono::high_resolution_clock::now();
            const auto is_prepared = aligner.prepare(line);
            auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - now);
            prep += diff.count();

            if (is_prepared) {
                now = std::chrono::high_resolution_clock::now();
                aligner.align();
                diff = std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::high_resolution_clock::now() - now);
                align_a += diff.count();

                now = std::chrono::high_resolution_clock::now();
                counter.count(read);
                diff = std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::high_resolution_clock::now() - now);
                count += diff.count();
            }
            std::getline(input, line);
        }
        counter.write_to_file(out_file);

        std::cout << "Prepare: " << prep << std::endl;
        std::cout << "Align: " << align_a << std::endl;
        std::cout << "Count: " << count << std::endl;
    }

    void analyse_positions(const std::string& ref,
                                  const std::string& sam_a,
                                  const std::string& out_file,
                                  const int dimension,
                                  const int qualityThreshold,
                                //   bool checkInsertions,
                                  bool checkDeletions,
                                  bool ambig) {

        std::ifstream input(sam_a);
        //TODO: check ob choose(ref.size(), dim) < max(size_t) und bestimme max dimension die möglich ist
        const auto ref_map = io_tools::read_reference(ref);

        //TODO Umbau: Interface für counter, dann mit count1-3 implementierten
//        ref::ref_map read;
//
//        aligner::aligner aligner{ref_map, read, qualityThreshold};
//
//        // record to time for each step
//        unsigned prep = 0;
//        unsigned align_a = 0;
//        unsigned count = 0;
//
//        std::string line = io_tools::get_first_data_line(input);
//        switch (dimension) {
//            case 1:
//                count::counter_1 counter{ref_map};
//                break;
//            case 2:
//                count::counter_2 counter{ref_map};
//                break;
//            case 3:
//                count::counter_3 counter{ref_map};
//                break;
//            default:
//                throw "";
//        }
//        while (input.good()) {
//            read.clear();
//            auto now = std::chrono::high_resolution_clock::now();
//            const auto is_prepared = aligner.prepare(line);
//            auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(
//                    std::chrono::high_resolution_clock::now() - now);
//            prep += diff.count();
//            if (is_prepared) {
//                now = std::chrono::high_resolution_clock::now();
//                aligner.align();
//                diff = std::chrono::duration_cast<std::chrono::nanoseconds>(
//                        std::chrono::high_resolution_clock::now() - now);
//                align_a += diff.count();
//
//                now = std::chrono::high_resolution_clock::now();
//                counter.count(read);
//                diff = std::chrono::duration_cast<std::chrono::nanoseconds>(
//                        std::chrono::high_resolution_clock::now() - now);
//                count += diff.count();
//            }
//            std::getline(input, line);
//        }
//        counter.write_to_file(out_file);
//
//        std::cout << "Prepare: " << prep << std::endl;
//        std::cout << "Align: " << align_a << std::endl;
//        std::cout << "Count: " << count << std::endl;
        std::cout << "\nStart counting" << std::endl;
        switch (dimension)
        {
            case 1:
                workflow_1(ref_map, input, out_file, qualityThreshold, 
                // checkInsertions, 
                checkDeletions, ambig);
                break;
            case 2:
                workflow_2(ref_map,input, out_file, qualityThreshold
                // checkInsertions, checkDeletions
                );
                break;
            case 3:
                workflow_3(ref_map,input, out_file, qualityThreshold
                // checkInsertions, checkDeletions
                );
                break;
            default:
                throw "";
        }

    }
}
