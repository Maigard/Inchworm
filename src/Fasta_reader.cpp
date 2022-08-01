#include "Fasta_reader.hpp"
#include "sequenceUtil.hpp"
#include "stacktrace.hpp"
#include <algorithm>
#include <iostream>
#include <string>
#include <omp.h>

//constructor
Fasta_reader::Fasta_reader (string filename) {
  
    //this->_hasNext = false;
    
    this->end_reading = -1; // turn off
    this->file_byte_pos = 0; // init
    this->_file_length = -1;

    if (filename == "-") {
        filename = "/dev/fd/0"; // read from stdin
    }
    this->_filereader.open(filename.c_str());
    if (! _filereader.is_open()) {
        throw(stacktrace() + "\n\nError, cannot open file " + filename );
    }
    if (filename != "-") {
        this->_filereader.seekg(0, ios_base::end);
        this->_file_length = this->_filereader.tellg();
        this->_filereader.seekg(this->file_byte_pos, _filereader.beg);
    }
        
    this->_init_reader();
   
}

Fasta_reader::Fasta_reader(string filename, long start_reading, long end_reading) {

    this->start_reading = this->file_byte_pos = start_reading;
    this->end_reading = end_reading;
    
    this->_filereader.open(filename.c_str());
    if (! _filereader.is_open()) {
        throw(stacktrace() + "\n\nError, cannot open file " + filename );
    }

    this->_filereader.seekg(0, ios_base::end);
    this->_file_length = this->_filereader.tellg();
    this->_filereader.seekg(start_reading, _filereader.beg);
    this->file_byte_pos = start_reading;
    
    this->_init_reader();
    
    // this->debugfile.open("/tmp/inchworm." + to_string(omp_get_thread_num()));
    // int count = this->countSequences();
    // this->debugfile << count << endl;
}

void Fasta_reader::_init_reader() {

    // primer reader to first fasta header
    getline(this->_filereader, this->_lastline);
    this->file_byte_pos += this->_lastline.length() + 1;
    this->debugfile << this->_lastline << endl;

    while ((! this->_filereader.eof()) && this->_lastline[0] != '>') {
        getline(this->_filereader, this->_lastline);
        this->file_byte_pos += this->_lastline.length() + 1;
        this->debugfile << this->_lastline << endl;
    }

}


bool Fasta_reader::hasNext() {
    bool ret;
    #pragma omp critical (FileReader)
    ret = this->hasNext_mt();
    return ret;
}

bool Fasta_reader::hasNext_mt() {
    bool ret;
    
    ret = !(this->_filereader.eof());
    
    // see if we're reading only part of the file.
    if (ret && this->end_reading > 0) {
        if (this->file_byte_pos >= end_reading) { // bad:  this->_filereader.tellg() >= end_reading) {
            // force it to go to the end of the file
            this->_filereader.seekg(0, this->_filereader.end);
            ret = false;
        }
    }
    return ret;
}

Fasta_entry Fasta_reader::getNext() {
    Fasta_entry fe;
    #pragma omp critical (FileReader)
    fe = this->getNext_mt();
    return fe;
}

// 1: this->_lastline = ">3"
// $6 = {_filereader = <incomplete type>, _lastline = ">3", end_reading = 5641, _file_length = 14524, eof = false, file_byte_pos = 5614}
// 1: this->_lastline = "CATCTACTTCCTCGAGCAGACAAAG"
// 1: this->_lastline = ">1"
// $4 = {_filereader = <incomplete type>, _lastline = ">1", end_reading = 5641, _file_length = 14524, eof = false, file_byte_pos = 5643}
Fasta_entry Fasta_reader::getNext_mt() {
    
    string sequence;
    string header;
    bool ret;

    header = this->_lastline;
    
    // ret == true
    ret = !(this->_filereader.eof());
    if (ret == true)
    {
        this->_lastline = "";
        while ((! this->_filereader.eof()) && this->_lastline[0] != '>') {
            getline(this->_filereader, this->_lastline);
            this->debugfile << this->_lastline << endl;
            
            // check if only reading section of a file
            if (this->end_reading > 0) {
                if (this->file_byte_pos >= end_reading && !this->_filereader.eof()) { // bad: this->_filereader.tellg() >= end_reading) {
                    // force it to go to the end of the file
                    this->_filereader.seekg(0, this->_filereader.end);
                }
            }
            this->file_byte_pos += this->_lastline.length() + 1;

            // cerr << "Comparing mine: " << this->file_byte_pos << " to " << this->_filereader.tellg() << endl;
            
            if (this->_lastline[0] != '>') {
                sequence += this->_lastline;
            }
        }
    }

    if (ret == true)
    {
        sequence = remove_whitespace(sequence);
        transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
        Fasta_entry fe(header, sequence);
        return(fe);
    } else {
        Fasta_entry fe("", "");
        return(fe);
    }
}

unsigned long Fasta_reader::countSequences() {
    int myTid = omp_get_thread_num();
    streampos savepos = _filereader.tellg();
    streampos startpos = this->start_reading;
    _filereader.seekg(startpos);
    streampos curpos = startpos;
    // ofstream outfile("/tmp/trinity" + std::to_string(omp_get_thread_num()) + ".txt"); 
    char buffer[1024*1024];
    long unsigned int to_read;
    long unsigned int sum = 0;
    int i=0;
    while(curpos < end_reading) {
        to_read = end_reading - curpos;
        if (sizeof(buffer) < to_read) {
            to_read = sizeof(buffer);
        }
        _filereader.read(buffer, to_read);
        // outfile << buffer;
        sum += count(buffer, buffer + to_read, '>');
        curpos += to_read;
    }
    _filereader.seekg(savepos);
    // sum += count(this->_lastline.begin(), this->_lastline.end(), '>');
    return sum;
}

map<string,string> Fasta_reader::retrieve_all_seqs_hash() {
    
    map<string,string> all_seqs_hash;
    
    while (this->hasNext()) {
        Fasta_entry f = this->getNext();
        string acc = f.get_accession();
        string seq = f.get_sequence();
        
        all_seqs_hash[acc] = seq;
    }
    
    return(all_seqs_hash);
}

long Fasta_reader::getFilelength() {
    return _file_length;
}
