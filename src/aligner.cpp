#include "aligner.hpp"
#include "utils.hpp"

#include <algorithm>

#include <iostream>
#include <bitset>

namespace
{
char to_upper(const char in)
{
    static char offsets[] = {32, 0};
    return in - offsets[in <= 'Z'];
}

// convert symbol of fastq quality string to phred+33 quality value
int get_quality(const char symbol)
{
    return ((int)symbol) - 33 ;
}
}


namespace aligner
{

    bool aligner::prepare(std::string& line_a,
                          std::string& line_b)
    {
        aligner::prepare(line_a);

        auto it = utils::find_tab_in_string<1>(line_b.begin());
        bool flag = utils::string_to_unsigned(it);
        std::bitset<3>flag_bin_b(flag);
        //if the flag includes 4 (= third bit set), the read is unmapped
        bool unmapped_b = flag_bin_b.test(2);

        //it = find_tab_in_string<3>(line_b.begin());
        it = utils::find_tab_in_string<2>(it);
        posinref_b = utils::string_to_unsigned(it);
        if (posinref_b == 0)
        {
            return false;
        }

        // getcigar containing the softclipped regions, indels and aligned regions (6th entry in SAM row)
        cigar_it_b = utils::find_tab_in_string<2>(it);
        if (unmapped_b or *cigar_it_b == '*')
        {
            return false;
        }

        cigar_end_b = utils::find_tab_in_string<1>(cigar_it_b) - 1;
        read_seq_b = utils::find_tab_in_string<4>(cigar_end_b);
        quality_seq_b = utils::find_tab_in_string<1>(read_seq_b);

        return true;
    }

    bool aligner::prepare(std::string& line)
    {
        //SAM doc: "Bit 0x4 (of the FLAG) is the only reliable place to tell whether the read is unmapped.
        //If 0x4 is set, no assumptions can be made about the rest.."
        auto it = utils::find_tab_in_string<1>(line.begin());
        unsigned flag = utils::string_to_unsigned(it);
        std::bitset<3>flag_bin(flag);
        //if the flag includes 4 (= third bit set), the read is unmapped
        bool unmapped = flag_bin.test(2);

        // get left most position in reference where read maps (4th entry in SAM row)
        //auto it = find_tab_in_string<3>(line_a.begin());
        it = utils::find_tab_in_string<2>(it);
        posinref_a = utils::string_to_unsigned(it);
        if (posinref_a == 0)
        {
            return false;
        }

        // getcigar containg the sofclipped regions, indels and aligned regions (6th entry in SAM row)
        cigar_it_a = utils::find_tab_in_string<2>(it);
        if (unmapped or *cigar_it_a == '*')
        {
            return false;
        }

        cigar_end_a = utils::find_tab_in_string<1>(cigar_it_a) - 1;
        read_seq_a = utils::find_tab_in_string<4>(cigar_end_a);
        quality_seq_a = utils::find_tab_in_string<1>(read_seq_a);

        return true;
    }


void aligner::align()
{
    aligning_started = false;

    while (cigar_it_a not_eq cigar_end_a)
    {
        const auto num = utils::string_to_unsigned(cigar_it_a);
        const auto c = to_upper(*cigar_it_a);
        ++cigar_it_a;

        if (c == 'M')
        {
            read.reserve(read.size() + num);
            aligning_started = true;
            for (unsigned i = 0; i < num; ++i)
            {
                const int quali = get_quality(*quality_seq_a);
                const char base = (*quality_seq_a != '*' &&  quali < quality_threshold) ? 'X' : to_upper(*read_seq_a);

                //TODO: erstmal rausnehmen, da wir für den Vergleich mit Red auch die rausschmeißen wo bei Überlappung utnerschiede sind
                // und wenn N und irgendeine Mutation vorhanden ist, soll sie ja sowie rausgeschmissen werden
//                if (base_id not_eq 'N')
//                {
                    read.add({posinref_a, nucleotide::nucleobase{base}});
//                }
                ++posinref_a;
                ++read_seq_a;
                //if only "*" is given, no quality sequence is available
                if(*quality_seq_a != '*')
                    ++quality_seq_a;
            }

        }
        else if (c == 'D')
        {
            posinref_a += num;
        } 
        else if (c == 'I' )
        {
            read_seq_a += num;
            //if only "*" is given, no quality sequence is available
            if(*quality_seq_a != '*')
                quality_seq_a += num;
        } 
        else if (c == 'S') 
        {
            read_seq_a += num;
            //if only "*" is given, no quality sequence is available
            if(*quality_seq_a != '*')
                quality_seq_a += num;
            if (aligning_started)
            {
                posinref_a += num;
            }
        }
    }
}

void aligner::align_1(count::counter_1& count_obj)
{
    aligning_started = false;

    auto pos = read.begin();

    while (cigar_it_b not_eq cigar_end_b)
    {
        const auto num = utils::string_to_unsigned(cigar_it_b);
        const char c = to_upper(*cigar_it_b);
        ++cigar_it_b;

        if (c == 'M')
        {
            aligning_started = true;
            for (unsigned i = 0; i < num; ++i)
            {
                const int quali = get_quality(*quality_seq_b);
                const char base = (*quality_seq_b != '*' && quali < quality_threshold) ? 'X' : to_upper(*read_seq_b);


               //if (base_id not_eq 'N')
               //{
                    //maybe two from last_pos to end, and then from begin to las_pos

                    // check if the paired reads overlap...
                    // ...for the the case the second read is behind the first one
                    pos = std::find_if(pos, read.end(), [this](const auto& val)
                                       {
                                       return val.first == this->posinref_b;
                                       });
                    // ...for the the case the second read is before the first one
                    if (pos == read.end())
                    {
                        pos = std::find_if(read.begin(), pos, [this](const auto& val)
                                           {
                                           return val.first == this->posinref_b;
                                           });
                    }



                    if (pos == read.end())
                    {
                        //no overlapping: count all nucleotides for the second read
                        //(is not stored in the read for later counting)
                        //if(base not_eq 'N')
                        if(nucleotide::isValidNucl(base, ambig))
                        {
                            count_obj.count(posinref_b - 1,
                                            base);
                        }
                    }
                    else
                    {
                        //overlapping: if the bases are equal-> ignore (will be counted later,
                        //as already stored in the read)
                        //if (pos->second.get() not_eq base or base=='N' or pos->second.get()=='N')
                        if(pos->second.get() not_eq base or !nucleotide::isValidNucl(base, ambig)
                        or !nucleotide::isValidNucl(pos->second.get(), ambig))
                        {
                            //if not equal-> check if one of the bases is equal to ref and count ref
                            // if different mutations: ignore
                            // if both have N, remove
                            const auto ref_base = ref.get(posinref_b - 1).get();
                            if (ref_base == base or ref_base == pos->second.get())
                            {
                                count_obj.count(posinref_b - 1, ref_base);
                            }
                            // remove to not count twice (as it's stored in the read)
                            read.remove(pos);
                        }
                    }
                //}
                ++posinref_b;
                ++read_seq_b;
                //if only "*" is given, no quality sequence is available
                if(*quality_seq_b != '*')
                    ++quality_seq_b;
            }

        }
        else if (c == 'D')
        {
            posinref_b += num;
        } 
        else if (c == 'I' )
        {
            read_seq_b += num;
            //if only "*" is given, no quality sequence is available
            if(*quality_seq_b != '*')
                quality_seq_b += num;
        } 
        else if (c == 'S') 
        {
            read_seq_b += num;
            //if only "*" is given, no quality sequence is available
            if(*quality_seq_b != '*')
                quality_seq_b += num;
            if (aligning_started)
            {
                posinref_b += num;
            }
        }
    }
}

void aligner::align_2()
{
    aligning_started = false;

    while (cigar_it_b not_eq cigar_end_b)
    {
        const auto num = utils::string_to_unsigned(cigar_it_b);
        const char c = to_upper(*cigar_it_b);
        ++cigar_it_b;

        if (c == 'M')
        {
            read.reserve(read.size() + num);

            auto pos = read.begin();

            aligning_started = true;
            for (unsigned i = 0; i < num; ++i)
            {
                const char base = to_upper(*read_seq_b);
//                if (base_id not_eq 'N')
//                {
                    //maybe two from last_pos to end, and then from begin to las_pos

                    // check if the paired reads overlap...
                    // ...for the case the second read is behind the first one
                    pos = std::find_if(pos, read.end(), [this](const auto& val)
                                       {
                                       return val.first == this->posinref_b;
                                       });

                    //...for the the case the second read is before the first one
                    if (pos == read.end())
                    {
                        pos = std::find_if(read.begin(), pos, [this](const auto& val)
                                           {
                                           return val.first == this->posinref_b;
                                           });
                    }

                    if (pos == read.end())
                    {
                        //no overlapping: count all nucleotides for the second read
                        //(is not stored in the read for later counting)
                        //if(base not_eq 'N')
                        if(nucleotide::isValidNucl(base, ambig))
                        {
                            read.add({posinref_b, nucleotide::nucleobase{base}});
                        }
                    }
                    else
                    {
                        //overlapping: if the bases are equal-> ignore (will be counted later,
                        //as already stored in the read)
                        //if (pos->second.get() not_eq base or base=='N' or pos->second.get()=='N')
                        if(pos->second.get() not_eq base or !nucleotide::isValidNucl(base, ambig)
                               or !nucleotide::isValidNucl(pos->second.get(), ambig))
                        {
                            //if not equal-> check if one of the bases is equal to ref and count ref
                            // if different mutations: ignore
                            const auto ref_base = ref.get(posinref_b - 1).get();
                            if (ref_base == base or ref_base == pos->second.get())
                            {
                                pos->second = ref.get(posinref_b - 1);
                            }else
                            {
                                // if different mutations or mutation and N: remove from read
                                read.remove(pos);
                                pos = read.begin();
                            }
                        }
                    }
//                }
                ++posinref_b;
                ++read_seq_b;
            }

        }
        else if (c == 'D')
        {
            posinref_b += num;
        } 
        else if (c == 'I' )
        {
            read_seq_b += num;
        } 
        else if (c == 'S') 
        {
            read_seq_b += num;
            if (aligning_started)
            {
                posinref_b += num;
            }
        }
    }
}

}

