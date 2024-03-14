/*
* This file is part of paraprobe-toolbox.
*
* paraprobe-toolbox is free software: you can redistribute it and/or modify it
* under the terms of the GNU General Public License as published by the
* Free Software Foundation, either version 3 of the License,
*  or (at your option) any later version.
*
* paraprobe-toolbox is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
* See the GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with paraprobe-toolbox. If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef __HDFFIVECORE_H__
#define __HDFFIVECORE_H__


//C++ STL
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <stdexcept>
#include <string>
#include <iomanip>
#include <limits>
#include <vector>
#include <algorithm>
#include <cmath>
#include <list>

#include <cassert>
#include <map>
#include <iterator>
#include <utility>
#include <set>
#include <unordered_set>
#include <type_traits>

#include <random>
#include <complex>

//low level stuff for interacting with the operating system
//#define NDEBUG
#include <assert.h>
#include <stdint.h>
#include <ios>

//C header stuff to pull low level system and process status pieces of information
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>

//##MK::for querying system maximum physical memory see
//https://stackoverflow.com/questions/2513505/how-to-get-available-memory-c-g
//https://stackoverflow.com/questions/850774/how-to-determine-the-hardware-cpu-and-ram-on-a-machine
//https://stackoverflow.com/questions/349889/how-do-you-determine-the-amount-of-linux-system-ram-in-c
//https://stackoverflow.com/questions/5553665/get-ram-system-size


//forward declaration for global scope
using namespace std;

//use functions from the HDF5 library
#include "hdf5.h"


#define MYHDF5_COMPRESSION_NONE		(0x00)
#define MYHDF5_COMPRESSION_DEFAULT	(0x01) //0x01 fastest but often still some head room, 0x09 slowest but higher
#define MYHDF5_COMPRESSION_GZIP		(0x01)


enum MYHDF5_ACCESSTYPES {
	READONLY,
	WRITE_NEW,
	WRITE_APPEND
};


enum MYHDF5_DATATYPES {
	H5_NONE,
	H5_U8,
	H5_I8,
	H5_U16,
	H5_I16,
	H5_U32,
	H5_I32,
	H5_U64,
	H5_I64,
	H5_F32,
	H5_F64,
	H5_STR
};

//comment in to obtain detailed reportings of the HDF5 library calls to hunt down bugs in the HDF5 wrapper
//#define MYHDF5_VERBOSE_LEVEL_TWO
//#define MYHDF5_VERBOSE_LEVEL_ONE

//return codes for wrapped H5 functions
#define MYHDF5_SUCCESS								+1 //MK::following the HDF5 convention that error values are positiv in case of success or negative else
#define MYHDF5_ALLOCERR								-1
#define MYHDF5_OUTOFBOUNDS							-2
#define MYHDF5_ARGINCONSISTENT						-3
#define MYHDF5_EXECUTIONORDERISSUE					-4
#define MYHDF5_INCORRECTLOGIC						-5
#define MYHDF5_INCORRECT_DIMS						-6
#define MYHDF5_INCORRECTOFFSETS						-7
#define MYHDF5_FAILED								-8
#define MYHDF5_RBUFALLOC_FAILED						-9
#define MYHDF5_READ_FAILED							-10
#define MYHDF5_DOESNOTEXIST							-11
#define MYHDF5_FILEACCESS_FAILED					-12
//#define MYHDF5_DSETACCESS_FAILED					-13
#define MYHDF5_MSPACE_FAILED						-14
#define MYHDF5_DSPACE_FAILED						-15
#define MYHDF5_HYPERSLAB_FAILED						-16
#define MYHDF5_GRPACCESS_FAILED						-17

#define MYHDF5_FCREATE_FAILED						-20
#define MYHDF5_FOPEN_FAILED							-21
#define MYHDF5_FTYPCLOSE_FAILED						-28
#define MYHDF5_FCLOSE_FAILED						-29

#define MYHDF5_GCREATE_FAILED						-30
#define MYHDF5_GOPEN_FAILED							-31
#define MYHDF5_GCLOSE_FAILED						-39

#define MYHDF5_DSPCCREATE_FAILED					-40
#define MYHDF5_DSPCOPEN_FAILED						-41
#define MYHDF5_DSPCACCESS_FAILED					-42
#define MYHDF5_DSPCREAD_FAILED						-43
#define MYHDF5_DSPCWRITE_FAILED						-48
#define MYHDF5_DSPCCLOSE_FAILED						-49

#define MYHDF5_DSETCREATE_FAILED					-50
#define MYHDF5_DSETOPEN_FAILED						-51
#define MYHDF5_DSETACCESS_FAILED					-52
#define MYHDF5_DSETREAD_FAILED						-53
#define MYHDF5_DSETWRITE_FAILED						-58
#define MYHDF5_DSETCLOSE_FAILED						-59

#define MYHDF5_MSPCCREATE_FAILED					-60
#define MYHDF5_MSPCOPEN_FAILED						-61
#define MYHDF5_MSPCCLOSE_FAILED						-69

#define MYHDF5_HSLBCREATE_FAILED					-70
#define MYHDF5_HSLBOPEN_FAILED						-71
#define MYHDF5_HSLBWRITE_FAILED						-75
#define MYHDF5_HSLBREAD_FAILED						-76
#define MYHDF5_HSLBCLOSE_FAILED						-79

#define MYHDF5_ATTRCREATE_FAILED					-70
#define MYHDF5_ATTROPEN_FAILED						-71
#define MYHDF5_ATTRACCESS_FAILED					-72
#define MYHDF5_ATTRREAD_FAILED						-73
#define MYHDF5_ATTRWCOLL_FAILED						-74
#define MYHDF5_ATTRWRITE_FAILED						-78
#define MYHDF5_ATTRCLOSE_FAILED						-79

#define MYHDF5_PLSTCREATE_FAILED					-80
#define MYHDF5_PLSTOPEN_FAILED						-81
#define MYHDF5_PLSTACCESS_FAILED					-82
#define MYHDF5_PLSTFILTER_FAILED					-83
#define MYHDF5_PLSTCLOSE_FAILED						-89

#define MYHDF5_ERROR_OBJ_EXISTS						-90
#define MYHDF5_ERROR_INPUT_SIZE						-91

#define MYHDF5_CMPGZIP_FAILED						-92
#define MYHDF5_CHKSET_FAILED						-93
#define MYHDF5_NOT_IMPLEMENTED						-94
#define MYHDF5_NOT_EXPECTED							-666


struct io_info
{
	//details
	vector<size_t> shape;				//d-dimensional shape, shape.size() == 0 for scalars, shape.size() == d for d-dimensional vector, matrix, tensor
	vector<size_t> chunk;				//if chunk.size() == 0 means dataset will use contiguous layout no chunking, if chunk.size() > 0 chunk.size() has to be == shape.size()
	unsigned char compression;
	unsigned char compression_opts;
	bool is_valid;
	io_info() : shape(vector<size_t>()), chunk(vector<size_t>()),
			compression(MYHDF5_COMPRESSION_NONE), compression_opts(0x00),
			is_valid(true) {}
	io_info( vector<size_t> const & _shp );
	io_info( vector<size_t> const & _shp, vector<size_t> const & _chk );
	io_info( vector<size_t> const & _shp, vector<size_t> const & _chk,
			const unsigned char _cmp, const unsigned int _cmp_opts );
};

ostream& operator << (ostream& in, io_info const & val);


class ioAttributes
{
public:
	//individual "dictionaries"
	//scalar
	map<string, unsigned char> u8;
	map<string, char> i8;
	map<string, unsigned short> u16;
	map<string, short> i16;
	map<string, unsigned int> u32;
	map<string, int> i32;
	map<string, unsigned long> u64;
	map<string, long> i64;
	map<string, float> f32;
	map<string, double> f64;
	map<string, string> str;
	//list_single, or list_long
	map<string, vector<unsigned char>> u8_arr;
	map<string, vector<char>> i8_arr;
	map<string, vector<unsigned short>> u16_arr;
	map<string, vector<short>> i16_arr;
	map<string, vector<unsigned int>> u32_arr;
	map<string, vector<int>> i32_arr;
	map<string, vector<unsigned long>> u64_arr;
	map<string, vector<long>> i64_arr;
	map<string, vector<float>> f32_arr;
	map<string, vector<double>> f64_arr;
	map<string, vector<string>> str_arr;

	map<string, bool> exists;

	ioAttributes();
	~ioAttributes();

	//handle scalar
	template<typename T> void add( const string keyword, const T value );
	void add( const string keyword, const string value );
	/*
	template<typename T> T query( const string keyword );
	string query( const string keyword );
	*/

	//handle list
	template<typename T> void add( const string keyword, vector<T> values );
	void add( const string keyword, vector<string> values );
	/*
	template<typename T> vector<T> query( const string keyword );
	vector<string> query( const string keyword );
	*/

	void report();
};


#define MYHDF5_SCALAR		0x00
#define MYHDF5_LIST_SINGLE	0x01
#define MYHDF5_LIST_LONG	0x02

#define MYHDF5_IS_FIXED		0x00
#define MYHDF5_IS_VARIABLE	0x01

#define MYHDF5_NULLTERM		0x00
#define MYHDF5_NULLPAD		0x01
#define MYHDF5_SPACEPAD		0x02

#define MYHDF5_ASCII		0x00
#define MYHDF5_UTF8			0x01


class HdfFiveSeqHdl
{
public:
	HdfFiveSeqHdl();
	HdfFiveSeqHdl( const string h5fn );
	~HdfFiveSeqHdl();

	//low-level, should not be called in tool executable code directly use high-level functions instead !
	vector<string> split_absolute_path( const string h5_absolute_path );
	int nexus_create();
	bool nexus_path_exists( const string h5_absolute_path );
	int nexus_open( const unsigned flags );
	int nexus_attribute_close( const int previous_status );
	int nexus_close( const int previous_status );

	//attributes
	//high-level
	int nexus_write_attributes( const string loc_name, ioAttributes const & attrs );
	//low-level
	int nexus_write_basic_scalar_attributes( const string loc_name, ioAttributes const & attrs );
	int nexus_write_string_scalar_attributes( const string loc_name, ioAttributes const & attrs );
	int nexus_write_basic_array_attributes( const string loc_name, ioAttributes const & attrs );
	int nexus_write_string_array_attributes( const string loc_name, ioAttributes const & attrs );

	/* int nexus_read_attributes( const string dsnm, ioAttributes & attrs ); */
	//low-level
	template<typename T> int nexus_write_attribute_value( const string keyword, const T val );
	int nexus_write_attribute_value( const string keyword, const string val );
	template<typename T> int nexus_write_attribute_value( const string keyword, vector<T> & val );
	int nexus_write_attribute_value( const string keyword, vector<string> & val );

	//datasets
	//high-level
	template<typename T> int nexus_read( const string dsnm, T & retval );
	int nexus_read( const string dsnm, string & retval );
	template<typename T> int nexus_read( const string dsnm, vector<T> & retval );
	int nexus_read( const string dsnm, vector<string> & retval );

	int nexus_write_group( const string grpnm, ioAttributes const & attrs );
	template<typename T> int nexus_write( const string dsnm, const T val, ioAttributes const & attrs );
	int nexus_write( const string dsnm, const string val, ioAttributes const & attrs );

	template<typename T> int nexus_write( const string dsnm, io_info const & ifo, vector<T> & val, ioAttributes const & attrs );
	int nexus_write( const string dsnm, io_info const & ifo, vector<string> & val, ioAttributes const & attrs );
	int nexus_write( const string dsnm, vector<string> & val, const unsigned char dimensionality,
			const unsigned char memory, const unsigned char terminator, const unsigned char encoding, ioAttributes const & attrs );
	int nexus_write_string_dataset( const string dsnm, vector<string> & val, const size_t max_characters,
			const unsigned char dimensionality, const unsigned char memory,
			const unsigned char terminator, const unsigned char encoding, ioAttributes const & attrs );

	hid_t adtypid;
	hid_t aspcid;
	hid_t attrid;

	hid_t mtypid;
	hid_t dtypid;
	hid_t dspcid;
	hid_t plistid;
	hid_t dsetid;
	hid_t fileid;
	string h5resultsfn;
	//ioAttributes attr_buf;
};

#endif
