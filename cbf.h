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

#pragma once

#ifndef CBF_H
#define CBF_H

// OS specific headers
#if defined(__MINGW32__) || defined(__MINGW64__)
#include <windef.h>
#include "win/libexport.h"
#elif defined(_MSC_VER)
#include <windows.h>
#include "win/libexport.h"
#define WIN32_LEAN_AND_MEAN
#elif __GNUC__
#include "linux/lindef.h"
#include "linux/libexport.h"
#endif

#include "end.h"

#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <vector>

#include "base64.h"




namespace cbf {
    long binomialCoeff(int n, int k);

	// The CBF class implementing the Spatial Bloom FIlters
	class DLL_PUBLIC CBF
	{

	private:
		BYTE *filter;
		BYTE ** HASH_salt;
		int bit_mapping;
		int cells;
		int cell_size;
		int size;
		int HASH_family;
		int HASH_number;
		int HASH_digest_length;
		int members;
        int unique_members;
		int MULTIPLICITY_max;
		std::vector<int> overflows;
		int BIG_end;

		// Private methods (commented in the cbf.cpp)
		void SetCell(unsigned int index, int area);
		int GetCell(unsigned int index) const;
		void CreateHashSalt(const std::string& path);
		void LoadHashSalt(const std::string& path);
		void SetHashDigestLength();
		void Hash(char *d, size_t n, unsigned char *md) const;


	public:
		// The maximum string (as a char array) length in bytes of each element
		// given as input to be mapped in the CBF
		const static int MAX_INPUT_SIZE = 128;
		// This value defines the maximum size (as in number of cells) of the CBF:
		// MAX_BIT_MAPPING = 32 states that the CBF will be composed at most by
		// 2^32 cells. The value is the number of bits used for CBF indexing.
		const static int MAX_BIT_MAPPING = 32;
		// Utility byte value of the above MAX_BIT_MAPPING
		const static int MAX_BYTE_MAPPING = MAX_BIT_MAPPING / 8;
		// The maximum value of the counters. This way, we limit the memory size
		// (which is the memory size of each cell) to 2 bytes
		const static int MAX_MULTIPLICITY = 65535;
		// The maximum number of allowed digests
		const static int MAX_HASH_NUMBER = 1024;

		// CBF class constructor
		// Arguments:
		// bit_mapping      actual size of the filter (as in number of cells): for
		//                  instance, bit_mapping = 10 states that the CBF will be
		//                  composed by 2^10 cells. As such, the size can only be
		//                  defined as a power of 2. This value is bounded by the
		//                  MAX_BIT_MAPPING constant.
		// HASH_family      specifies the hash function to be used. Currently available:
		//                  1: SHA1
		//                  4: MD4
		//                  5: MD5
		// HASH_number      number of digests to be produced (by running the hash function
		//                  specified above using different salts) in the insertion and
		//                  check phases.
		// MULTIPLICITY_max maximum multiplicity found in the construction dataset.
		// salt_path        path to the file where to read from/write to the hash salts
		//                  (which are used to produce different digests).
		//                  If the file exists, reads one salt per line.
		//                  If the file doesn't exist, the salts are randomly generated
		//                  during the filter creation phase
		CBF(int bit_mapping, int HASH_family, int HASH_number, int MULTIPLICITY_max,
		        const std::string& salt_path, int forced_cell_size=0)
		{

			// Argument validation
			if (bit_mapping <= 0 || bit_mapping > MAX_BIT_MAPPING) throw std::invalid_argument("Invalid bit mapping.");
			if (MULTIPLICITY_max <= 0 || MULTIPLICITY_max > MAX_MULTIPLICITY) throw std::invalid_argument("Invalid multipliciy value.");
			if (HASH_number <= 0 || HASH_number > MAX_HASH_NUMBER) throw std::invalid_argument("Invalid number of hash runs.");
			if (salt_path.length() == 0) throw std::invalid_argument("Invalid hash salt path.");

			// Checks whether the execution is being performed on a big endian or little endian machine
			this->BIG_end = cbf::is_big_endian();

			// Defines the number of bytes required for each cell depending on MULTIPLICITY_max
			// In order to reduce the memory footprint of the filter, we use 1 byte 
			// for a maximum multiplicity <= 255, 2 bytes for a up to MAX_MULTIPLICITY
			if(forced_cell_size > 0) {
			    if(forced_cell_size > 2) {
			        throw std::invalid_argument("Forced cell size must be 1 or 2");
                }

                this->cell_size = forced_cell_size;
            } else {
                if (MULTIPLICITY_max <= 255) this->cell_size = 1;
                else if (MULTIPLICITY_max > 255) this->cell_size = 2;
			}


			// Sets the type of hash function to be used
			this->HASH_family = HASH_family;
			this->SetHashDigestLength();
			// Sets the number of digests
			this->HASH_number = HASH_number;

			// Initializes the HASH_salt matrix
			this->HASH_salt = new BYTE*[HASH_number];
			for (int j = 0; j<HASH_number; j++) {
				this->HASH_salt[j] = new BYTE[CBF::MAX_INPUT_SIZE];
			}

			// Creates the hash salts or loads them from the specified file
			std::ifstream my_file(salt_path.c_str());
			if (my_file.good()) this->LoadHashSalt(salt_path);
			else this->CreateHashSalt(salt_path);

			// Defines the number of cells in the filter
			this->cells = (int)pow(2, bit_mapping);
			this->bit_mapping = bit_mapping;

			// Defines the total size in bytes of the filter
			this->size = this->cell_size*this->cells;

			// Memory allocation for the CBF array
			this->filter = new BYTE[this->size];

			// Initializes the cells to 0
			for (int i = 0; i < this->size; i++) {
				this->filter[i] = 0;
			}

			// Initializes the overflow counter
            overflows = std::vector<int>(this->cells, 0);

            // Initializes the members counters
            this->members = 0;
            this->unique_members = 0;

			// Sets the maximum multiplicity found in the construction dataset
			this->MULTIPLICITY_max = MULTIPLICITY_max;
		}

		// CBF class destructor
		~CBF()
		{
			// Frees the allocated memory
			delete[] filter;
			for (int j = 0; j<this->HASH_number; j++) {
				delete[] HASH_salt[j];
			}
			delete[] HASH_salt;
		}


		// Public methods (commented in the cbf.cpp)
		void PrintFilter(int mode) const;
		void SaveToDisk(const std::string& path, int mode);
		void Insert(const char *string, int size, int area);
		int Check(const char *string, int size) const;
		float GetFilterSparsity() const;
		float GetFilterFpp() const;
		float GetFilterAPrioriFpp() const;
        long double GetCellAPrioriOverflow() const;
		int GetOverallOverflows() const;
        int GetOverflownCells() const;
	};

} //namespace cbf

#endif /* CBF_H */