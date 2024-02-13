#ifndef MIME_IO_HPP_
#define MIME_IO_HPP_

#include <string>

namespace mime_io_paired
{
    void analyse_positions(const std::string& ref,
                                  const std::string& sam_a,
                                  const std::string& sam_b,
                                  const std::string& out_file,
                                  const int dimension,
                                  const int qualityThreshold,
                                //   bool checkInsertions,
                                  bool checkDeletions);

}

#endif
