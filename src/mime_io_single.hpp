//
// Created by Smith, Maureen on 16.04.20.
//

#ifndef SAM2COUNTS_MIME_IO_SINGLE_HPP
#define SAM2COUNTS_MIME_IO_SINGLE_HPP

#include <string>


namespace mime_io_single {
    void analyse_positions(const std::string& ref,
                                  const std::string& sam_a,
                                  const std::string& out_file,
                                  const int dimension,
                                  const int qualityThreshold,
                                //   bool checkInsertions,
                                  bool checkDeletions,
                                  bool ambig);
}


#endif //SAM2COUNTS_MIME_IO_SINGLE_HPP
