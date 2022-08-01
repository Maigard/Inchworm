#ifndef __FASTA_READER__
#define __FASTA_READER__

#include <iostream>
#include <fstream>
#include <string>
#include "Fasta_entry.hpp"
#include <map>

using namespace std;

class Fasta_reader {
    
public:
    Fasta_reader(string filename); // constructor
    Fasta_reader(string filename, long start_reading, long end_reading); // constructor
    
    bool hasNext();  // prime for next 
    Fasta_entry getNext(); // retrieves next fasta entry, or NULL if EOF.

    bool hasNext_mt();  // hasNext without locking
    Fasta_entry getNext_mt(); // getNext without locking

    unsigned long countSequences();
    
    map<string,string> retrieve_all_seqs_hash();

    long getFilelength();
    
    
    
private:
    ifstream _filereader;
    //bool _hasNext;
    string _lastline;
    long start_reading; // optional file position to stop reading.
    long end_reading; // optional file position to stop reading.
    long _file_length;

    long file_byte_pos;


    void _init_reader();

    ofstream debugfile;
};


#endif
