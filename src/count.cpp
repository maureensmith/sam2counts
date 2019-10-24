#include "count.hpp"

#include <algorithm>
#include <fstream>
#include <vector>
#include <iterator>

#include <iostream>


namespace count
{

    void counter_1::count(const ref::ref_map& read) {
        counter_1::count(read,1);
    }
    void counter_1::count(const ref::ref_map& read, const unsigned times)
    {
        std::for_each(read.begin(), read.end(), [this, times](const auto& entry)
        {
            //TODO: HIER SCHAUEN WELCHE READS GEZÄHLT WERDEN
           // if(entry.first == 17 && entry.second.get() == 'A')
    //            std::cout << "read " << i << " all from read pair1 " << entry.second.to_id() << std::endl;
            if(entry.second.get() not_eq 'N')
            {
                data[entry.first - 1][entry.second.to_id()] += times;
            }
        });
    }


void counter_1::write_to_file(const std::string& out_file)
{
    std::ofstream outfile(out_file);

    if (outfile.good())
    {
        outfile << "pos1\tA\tC\tG\tT\n";
        for (unsigned i = 0; i < data.size(); ++i)
        {
            outfile << (i + 1);
            std::for_each(data[i].cbegin(), data[i].cend(), [&outfile](const auto& entry)
            {
                outfile << '\t' << entry;
            });
            outfile << '\n';
        }
    }
}

void counter_2::count(const ref::ref_map& read) {
    counter_2::count(read,1);
}

void counter_2::count(const ref::ref_map& read, const unsigned times)
{
    if (read.size() < 2)
    {
        return;
    }

    auto pos1 = read.begin();
    const auto end = std::prev(read.end());
    const auto size_x = size - 1;
    while (pos1 != end)
    {
        const auto pos1_idx = (pos1->first - 1);
        auto nucl1 = pos1->second;
        //pairwise position index start
        const auto i = pos1_idx * size_x - utils::choose(pos1_idx, 2);
        if(pos1->second.get() not_eq 'N')
        {
            auto pos2 = ++pos1;
            while (pos2 != read.end())
            {
                if(pos2->second.get() not_eq 'N')
                {
                    //const auto id = i + (pos2->first - pos1_idx - 2);
                    // pairwise substitution index
                    //TODO prüfen
                   /* if(pos1_idx+1 == 153 && pos2->first == 161)
                    {
                        //std::cout << " count pos1 " << pos1_idx+1 << " " <<  nucl1.get() << " pos2 "<< pos2->first << " " << pos2->second.get() << std::endl;
//                        std::cout << "id " <<  i + (pos2->first - pos1_idx - 2) << std::endl;
                    } */

                    const auto j = counter_1::nucleobase_count * nucl1.to_id();
                    data[i + pos2->first - pos1_idx - 2][j + pos2->second.to_id()] += times;
                }
                ++pos2;
            }
        } else
        {
            ++pos1;
        }
    }
}




void counter_2::write_to_file(const std::string& out_file)
{
    std::ofstream outfile(out_file);

    if (outfile.good())
    {
        outfile << "pos1\tpos2\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\n";

        unsigned line = 0;
        unsigned i = 1;
        unsigned j = 2;
        while (line < data.size())
        {
            outfile << i << '\t' << j;
            std::for_each(data[line].cbegin(), data[line].cend(), [&outfile](const auto& entry)
            {
                outfile << '\t' << entry;
            });
            outfile << '\n';

            ++line;
            if (j == size)
            {
                ++i;
                j = i + 1;
            }
            else
            {
                ++j;
            }
        }
    }
}

//TODO anpassen
void counter_3::count(const ref::ref_map& read) {
    counter_3::count(read,1);
}

void counter_3::count(const ref::ref_map& read, const unsigned times)
{
    if (read.size() < 3)
    {
        return;
    }

    auto pos1 = read.begin();
    const auto end1 = std::prev(read.end(),2);
    //const auto size_x = size - 1;
//    std::cout << "size " << size << std::endl;
//    std::cout << "data size " << data.size() << std::endl;
    while (pos1 != end1)
    {
        //const auto pos1_idx = (pos1->first - 1);

        if(pos1->second.get() not_eq 'N')
        {
            auto p1 = pos1->first;
            auto nucl1 = pos1->second.to_id();
            //index start of the triplet, 1,2,3 is 1, 1,2,4 is 2, ...
            const auto seqpos1 = utils::choose(size,3) - utils::choose(size-p1+1, 3);
//            std::cout << "pos1 " << p1 << " seqpos 1 " << seqpos1 << std::endl;
            const auto end2 = std::prev(read.end());
            auto pos2 = ++pos1;
            while (pos2 != end2)
            {
                if(pos2->second.get() not_eq 'N')
                {
                    auto p2 = pos2->first;
                    auto nucl2 = pos2->second.to_id();
//                    std::cout << "p2 " << p2 << std::endl;
                    // indexstart of position 2
                    auto seqpos2 = (p2-p1-1)*(size-p1) - utils::choose(p2-p1,2);
//                    std::cout << "pos2 " << p2 << " seqpos 2 " << seqpos2 << std::endl;
                    auto pos3 = ++pos2;
                    while(pos3 != read.end())
                    {
                        if(pos3->second.get() not_eq 'N')
                        {
                            // pairwise substitution index
                            const auto mutPos = (counter_1::nucleobase_count^2) * nucl1
                                    + counter_1::nucleobase_count * nucl2
                                    + pos3->second.to_id();
                            //TODO
//                             std::cout << "pos1 " << p1 << " pos2 " << p2 << " pos3 " << pos3->first << std::endl;
//                             std::cout << "3D seqPos " << seqpos1 + seqpos2 + pos3->first - p2 - 1 << std::endl;
                             data[seqpos1 + seqpos2 + pos3->first - p2 - 1][mutPos] += times;
                        }
                        ++pos3;
                    }
                } else {
                    ++pos2;
                }
            }
        } else
        {
            ++pos1;
        }
    }
}

//TODO anpassen
void counter_3::write_to_file(const std::string& out_file)
{

    unsigned line = 0;
    unsigned i = 1;
    unsigned j = 2;
    unsigned k = 3;

    const std::string header = "pos1\tpos2\tpos3\tAAA\tAAC\tAAG\tAAT\tACA\tACC\tACG\tACT\tAGA\tAGC\tAGG\tAGT\tATA\tATC\tATG\tATT"
           "\tCAA\tCAC\tCAG\tCAT\tCCA\tCCC\tCCG\tCCT\tAGA\tCGC\tCGG\tCGT\tCTA\tCTC\tCTG\tCTT"
            "\tGAA\tGAC\tGAG\tGAT\tGCA\tGCC\tGCG\tGCT\tAGA\tGGC\tGGG\tGGT\tGTA\tGTC\tGTG\tGTT"
            "\tTAA\tTAC\tTAG\tTAT\tTCA\tTCC\tTCG\tTCT\tAGA\tTGC\tTGG\tTGT\tTTA\tTTC\tTTG\tTTT\n";

    std::ofstream outfile(out_file+"_"+std::to_string(i) +"_" +std::to_string(j));
    if (outfile.good())
    {
        outfile << header;

        while (line < data.size())
        {
            outfile << i << '\t' << j << '\t' << k;
            std::for_each(data[line].cbegin(), data[line].cend(), [&outfile](const auto& entry)
            {
                outfile << '\t' << entry;
            });
            outfile << '\n';

            ++line;
            if(k == size)
            {
                if(j==size-1){
                    ++i;
                    j = i + 1;
                }
                else {
                    ++j;
                }
                outfile.close();
                outfile.clear();
                outfile.open(out_file+"_"+std::to_string(i) +"_" +std::to_string(j));
                if (outfile.good())
                {
                    outfile << header;
                } else
                {
                    //TODO: fehlerbehandlung?
                    break;
                }

                k = j+1;


            }else
            {
                ++k;
            }
        }
    }
}

}

