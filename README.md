## libCBF-cpp ##
libCBF - implementation of the Counting Bloom Filters in C++.

The Couting Bloom Filters (CBF) extend the original Bloom filter concept by allowing multiplicity queries.

## Implementation ##
This library is based on the [Spatial Bloom Filters implementation](https://github.com/spatialbloomfilter/libSBF-cpp) by Luca Calderoni, Dario Maio and Paolo Palmieri.

Please refer to the linked repository for more details.

## License ##
The source code is released under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. Please refer to the included files [COPYING](COPYING) and [COPYING.LESSER](COPYING.LESSER).

## Acknowledgements ##
This project uses, where possible, existing libraries. In particular:
- We use the OpenSSL implementation of hash functions, copyright of The OpenSSL Project. Please refer to https://www.openssl.org/source/license.txt for OpenSSL project licence details.
- The functions implementing base64 encoding and decoding (provided in base64.cpp and base64.h) was written by [René Nyffenegger](mailto:rene.nyffenegger@adp-gmbh.ch). A full copyright statement is provided within each file.
