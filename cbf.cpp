/*
    Counting Bloom Filter C++ Library (libCBF-cpp)

    Copyright (C) 2020 Lorenzo Pellegrini
    University of Bologna

    Based on Spatial Bloom Filter C++ Library (https://github.com/spatialbloomfilter/libSBF-cpp)
    Copyright (C) 2017  Luca Calderoni, Dario Maio,
    University of Bologna
    Copyright (C) 2017  Paolo Palmieri,
    Cranfield University

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define CBF_DLL

#include "cbf.h"

#include <iostream>
#include <stdexcept>
#include<bits/stdc++.h>

#include <openssl/md4.h>
#include <openssl/md5.h>
#include <openssl/rand.h>
#include <openssl/sha.h>


namespace cbf {

/* **************************** PRIVATE METHODS **************************** */


    // Sets the hash digest length depending on the selected hash function
    void CBF::SetHashDigestLength() {
        switch (this->HASH_family) {
            case 1:
                this->HASH_digest_length = SHA_DIGEST_LENGTH;
                break;
            case 4:
                this->HASH_digest_length = MD4_DIGEST_LENGTH;
                break;
            case 5:
                this->HASH_digest_length = MD5_DIGEST_LENGTH;
                break;
            default:
                this->HASH_digest_length = MD4_DIGEST_LENGTH;
                break;
        }
    }


    // Computes the hash digest, calling the selected hash function
    // char *d            is the input of the hash value
    // size_t n           is the input length
    // unsigned char *md  is where the output should be written
    void CBF::Hash(char *d, size_t n, unsigned char *md) const {
        switch (this->HASH_family) {
            case 1:
                SHA1((unsigned char *) d, n, (unsigned char *) md);
                break;
            case 4:
                MD4((unsigned char *) d, n, (unsigned char *) md);
                break;
            case 5:
                MD5((unsigned char *) d, n, (unsigned char *) md);
                break;
            default:
                MD4((unsigned char *) d, n, (unsigned char *) md);
                break;
        }
    }


    // Stores a hash salt byte array for each hash (the number of hashes is
    // HASH_number). Each input element will be combined with the salt via XOR, by
    // the Insert and Check methods. The length of salts is MAX_INPUT_SIZE bytes.
    // Hashes are stored encoded in base64.
    void CBF::CreateHashSalt(const std::string &path) {
        BYTE buffer[CBF::MAX_INPUT_SIZE];
        int rc;
        std::ofstream myfile;

        myfile.open(path.c_str());

        for (int i = 0; i < this->HASH_number; i++) {
            rc = RAND_bytes(buffer, sizeof(buffer));
            if (rc != 1) {
                throw std::runtime_error("Failed to generate hash salt");
            }

            // Fills hash salt matrix
            memcpy(this->HASH_salt[i], buffer, CBF::MAX_INPUT_SIZE);
            // Writes hash salt to disk to the path given in input
            std::string encoded = cbf::base64_encode(reinterpret_cast<const unsigned char *>(this->HASH_salt[i]),
                                                     CBF::MAX_INPUT_SIZE);
            myfile << encoded << std::endl;
        }

        myfile.close();
    }


    // Loads from the path in input a hash salt byte array, one line per hash.
    // Hashes are stored encoded in base64, and need to be decoded.
    void CBF::LoadHashSalt(const std::string &path) {
        std::ifstream myfile;
        std::string line;

        myfile.open(path.c_str());

        for (int i = 0; i < this->HASH_number; i++) {
            // Reads one base64 hash salt from file (one per line)
            getline(myfile, line);

            //decode and fill hash salt matrix
            memcpy(this->HASH_salt[i], cbf::base64_decode(line).c_str(), CBF::MAX_INPUT_SIZE);
        }

        myfile.close();
    }


    // Sets the cell by incrementing the cell counter. This method is called
    // by Insert with the cell index and the multiplicity. It manages the two
    // different possible cell sizes (one or two bytes) automatically set during
    // filter construction.
    void CBF::SetCell(unsigned int index, int multiplicity) {
        int max_multiplicity = this->cell_size == 1 ? 255 : 65535;
        int cell_value = this->GetCell(index);
        int new_cell_value = cell_value + multiplicity;
        int n_overflows = new_cell_value - max_multiplicity;

        if ((multiplicity > max_multiplicity) || (multiplicity < 0)) {
            std::string error_message = "Multiplicity must be in [1, ";
            error_message += std::to_string(max_multiplicity);
            error_message += "]\n";
            throw std::invalid_argument(error_message);
        }

        if (n_overflows > 0) {
            overflows[index] += n_overflows;
        }

        new_cell_value = std::min(max_multiplicity, new_cell_value);

        switch (this->cell_size) {
            // 1-byte cell size
            case 1:
                this->filter[index] = new_cell_value;
                break;
                // 2-bytes cell size. Writing values over the two bytes is managed
                // manually, by copying byte per byte
            case 2:

                this->filter[2 * index] = (BYTE) (new_cell_value >> 8);
                this->filter[(2 * index) + 1] = (BYTE) new_cell_value;
                break;
            default:
                break;
        }
    }


    // Returns the counter stored at the specified index
    int CBF::GetCell(unsigned int index) const {
        int counter;
        switch (this->cell_size) {
            // 1-byte cell size
            case 1:
                return (int) this->filter[index];
                break;
                // 2-bytes cell size:
                // here we need to copy values byte per byte
            case 2:
                // Copies counter (one byte at a time) in two adjacent bytes
                counter = (int) ((this->filter[2 * index] << 8) | this->filter[(2 * index) + 1]);
                return counter;
                break;
            default:
                // This condition should never be reached
                return -1;
                break;
        }
    }


/* ***************************** PUBLIC METHODS ***************************** */


    // Prints the filter and related statistics to the standard output
    // mode: 0    prints CBF stats only
    // mode: 1    prints CBF information and the full CBF content
    void CBF::PrintFilter(const int mode) const {
        int potential_elements;

        printf("Counting Bloom Filter details:\n\n");

        printf("HASH details:\n");
        printf("Hash family: %d\n", this->HASH_family);
        printf("Number of hash runs: %d\n\n", this->HASH_number);

        printf("Filter details:\n");
        printf("Number of cells: %d\n", this->cells);
        printf("Size in Bytes: %d\n", this->size);
        printf("Filter sparsity: %.5f\n", this->GetFilterSparsity());
        printf("Filter a-priori fpp: %.5f\n", this->GetFilterAPrioriFpp());
        printf("Filter fpp: %.5f\n", this->GetFilterFpp());
        printf("Number of mapped elements: %d\n", this->members);
        printf("Number of unique elements: %d\n", this->unique_members);
        printf("Cell a-priori overflow probability: %Le\n", this->GetCellAPrioriOverflow());
        printf("Number of overflows: %d\n", this->GetOverallOverflows());
        printf("Number of overflown cells: %d\n", this->GetOverflownCells());

        if (mode == 1) {
            printf("\nFilter cells content:");
            for (int i = 0; i < this->size; i += this->cell_size) {
                // For readability purposes, we print a line break after 32 cells
                if (i % (32 * this->cell_size) == 0)printf("\n");
                switch (this->cell_size) {
                    case 1:
                        std::cout << (unsigned int) (filter[i]) << "|";
                        break;
                    case 2:
                        unsigned long value;
                        value = (filter[i] << 8) | (filter[i + 1]);
                        std::cout << (unsigned int) value << "|";
                        break;
                    default:
                        break;
                }
            }
            printf("\n\n");
        } else printf("\n");

        /*std::cout << "Overflows:" << std::endl;
        for (int i = 0; i < this->cells; i++) {
            std::cout << "Cell " << i << ": " << this->overflows[i] << std::endl;
        }*/

        printf("\n");
    }


    // Prints the filter and related statistics onto a CSV file (path)
    // mode: 1    writes CBF metadata (CSV: key;value)
    // mode: 0    writes CBF cells (CSV: value)
    void CBF::SaveToDisk(const std::string &path, int mode) {
        std::ofstream myfile;

        myfile.open(path.c_str());

        myfile.setf(std::ios_base::fixed, std::ios_base::floatfield);
        myfile.precision(5);

        if (mode) {

            myfile << "hash_family" << ";" << this->HASH_family << std::endl;
            myfile << "hash_number" << ";" << this->HASH_number << std::endl;
            myfile << "max_multiplicity" << ";" << this->MULTIPLICITY_max << std::endl;
            myfile << "bit_mapping" << ";" << this->bit_mapping << std::endl;
            myfile << "cells_number" << ";" << this->cells << std::endl;
            myfile << "cell_size" << ";" << this->cell_size << std::endl;
            myfile << "byte_size" << ";" << this->size << std::endl;
            myfile << "members" << ";" << this->members << std::endl;
            myfile << "unique_members" << ";" << this->unique_members << std::endl;
            myfile << "overflows" << ";" << this->GetOverallOverflows() << std::endl;
            myfile << "overflown_cells" << ";" << this->GetOverflownCells() << std::endl;
            myfile << "sparsity" << ";" << this->GetFilterSparsity() << std::endl;
            myfile << "a-priori fpp" << ";" << this->GetFilterAPrioriFpp() << std::endl;
            myfile << "fpp" << ";" << this->GetFilterFpp() << std::endl;
            myfile << "a-priori overflow" << ";" << this->GetCellAPrioriOverflow() << std::endl;

            for (int i = 0; i < this->cells; i++) {
                myfile << "overflows_" << i << ";" << this->overflows[i] << std::endl;
            }
        } else {
            for (int i = 0; i < this->size; i += this->cell_size) {
                switch (this->cell_size) {
                    case 1:
                        myfile << (unsigned int) (filter[i]) << std::endl;
                        break;
                    case 2:
                        unsigned long value;
                        value = (filter[i] << 8) | (filter[i + 1]);
                        myfile << (unsigned int) value << std::endl;
                        break;
                    default:
                        break;
                }
            }
        }

        myfile.close();
    }


    // Maps a single element (passed as a char array) to the CBF. For each hash
    // function, internal method SetCell is called, passing elements coupled with
    // its multiplicity.
    // char *string     element to be mapped
    // int size         length of the element
    // int multiplicity the element multiplicity
    void CBF::Insert(const char *string, const int size, const int multiplicity) {
        char *buffer = new char[size];

        // We allow a maximum CBF mapping of 32 bit (resulting in 2^32 cells).
        // Thus, the hash digest is limited to the first four bytes.
        unsigned char digest32[CBF::MAX_BYTE_MAPPING];

        unsigned char *digest = new unsigned char[this->HASH_digest_length];

        // Computes the hash digest of the input 'HASH_number' times; each
        // iteration combines the input char array with a different hash salt
        for (int k = 0; k < this->HASH_number; k++) {

            for (int j = 0; j < size; j++) {
                buffer[j] = (char) (string[j] ^ this->HASH_salt[k][j]);
            }

            this->Hash(buffer, size, (unsigned char *) digest);

            // Truncates the digest after the first 32 bits (see above)
            for (int i = 0; i < CBF::MAX_BYTE_MAPPING; i++) {
                digest32[i] = digest[i];
            }

            // Copies the truncated digest (one byte at a time) in an integer
            // variable (endian independent)
            unsigned int digest_index;
            if (this->BIG_end) {
                digest_index = (digest32[0] << 24) | (digest32[1] << 16) | (digest32[2] << 8) | digest32[3];
            } else {
                digest_index = (digest32[3] << 24) | (digest32[2] << 16) | (digest32[1] << 8) | digest32[0];
            }


            // Shifts bits in order to preserve only the first 'bit_mapping'
            // least significant bits
            digest_index >>= (CBF::MAX_BIT_MAPPING - this->bit_mapping);

            this->SetCell(digest_index, multiplicity);
        }

        this->unique_members++;
        this->members += multiplicity;

        delete[] buffer;
        delete[] digest;
    }

    // Verifies weather the input element belongs to the set.
    // Returns the counter (i.e. the minimum cell number) if the element belongs to a set, 0 otherwise.
    // char *string     the element to be verified
    // int size         length of the element
    int CBF::Check(const char *string, const int size) const {
        std::vector<char> buffer(size);
        int counter = INT_MAX;
        int current_counter = 0;

        // We allow a maximum CBF mapping of 32 bit (resulting in 2^32 cells).
        // Thus, the hash digest is limited to the first four bytes.
        unsigned char digest32[CBF::MAX_BYTE_MAPPING];

        std::vector<unsigned char> digest(this->HASH_digest_length);

        // Computes the hash digest of the input 'HASH_number' times; each
        // iteration combines the input char array with a different hash salt
        for (int k = 0; k < this->HASH_number; k++) {

            for (int j = 0; j < size; j++) {
                buffer[j] = (char) (string[j] ^ this->HASH_salt[k][j]);
            }

            this->Hash(buffer.data(), size, digest.data());

            // Truncates the digest to the first 32 bits
            for (int i = 0; i < CBF::MAX_BYTE_MAPPING; i++) {
                digest32[i] = digest[i];
            }

            // Copies the truncated digest (one byte at a time) in an integer
            // variable (endian independent)
            unsigned int digest_index;
            if (this->BIG_end) {
                digest_index = (digest32[0] << 24) | (digest32[1] << 16) | (digest32[2] << 8) | digest32[3];
            } else {
                digest_index = (digest32[3] << 24) | (digest32[2] << 16) | (digest32[1] << 8) | digest32[0];
            }

            // Shifts bits in order to preserve only the first 'bit_mapping' least
            // significant bits
            digest_index >>= (CBF::MAX_BIT_MAPPING - this->bit_mapping);

            current_counter = this->GetCell(digest_index);

            counter = std::min(counter, current_counter);
            // If one hash points to an empty cell, the element does not belong
            // to any set.
            if (counter == 0) break;
        }

        return counter;
    }

    // Returns the sparsity of the entire CBF
    float CBF::GetFilterSparsity() const {
        float ret;
        int e = 0;
        for (int i = 1; i < this->cells; i++) {
            if (this->GetCell(i) != 0) {
                e++;
            }
        }
        ret = ((float) e / (float) this->cells);

        return ret;
    }

    //https://www.geeksforgeeks.org/binomial-coefficient-dp-9/
    long binomialCoeff(int n, int k)
    {
        long C[k+1];
        memset(C, 0, sizeof(C));

        C[0] = 1;  // nC0 is 1

        for (int i = 1; i <= n; i++)
        {
            // Compute next row of pascal triangle using
            // the previous row
            for (int j = std::min(i, k); j > 0; j--)
                C[j] = C[j] + C[j-1];
        }
        return C[k];
    }

    // Returns the a-priori overflow probability of a cell
    long double CBF::GetCellAPrioriOverflow() const {
        int j = 255;
        if (this->cell_size == 2) {
            j = 65535;
        }

        int m = this->cells;
        int k = this->HASH_number;
        int n = this->members;

        /* Only for testing purpose
         * See: Ficara et al. "Multilayer Compressed Counting Bloom Filters"
            int k = 10;
            int n = 1000;
            int m = 14427;
            j = 8;
        */

        int kn = k * n;

        long double p = exp(1.0);
        p *= (long double) kn;
        p /= (long double) m * j;
        p = std::pow(p, j);
        return p;
    }


    // Returns the a-priori false positive probability over the entire filter
    float CBF::GetFilterAPrioriFpp() const {
        double p;

        p = (double) (1 - 1 / (double) this->cells);
        p = (double) (1 - (double) pow(p, this->HASH_number * this->unique_members));
        p = (double) pow(p, this->HASH_number);

        return (float) p;
    }


    // Returns the a-posteriori false positive probability over the entire filter
    float CBF::GetFilterFpp() const {
        double p;
        int c = 0;
        // Counts non-zero cells
        for (int i = 1; i < this->cells; i++) {
            if (this->GetCell(i) > 0) {
                c++;
            }
        }
        p = (double) c / (double) this->cells;

        p = (double) (pow(p, this->HASH_number));

        return (float) p;
    }

    // Returns the overall number of overflows
    int CBF::GetOverallOverflows() const {
        int total = 0;
        for (auto const &value: this->overflows) {
            total += value;
        }

        return total;
    }

    // Returns the number of overflown cells
    int CBF::GetOverflownCells() const {
        int total = 0;
        for (auto const &value: this->overflows) {
            if (value != 0) {
                total++;
            }
        }

        return total;
    }

} //namespace cbf