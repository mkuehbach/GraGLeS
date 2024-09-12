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


#include "PARAPROBE_HDF5Core.h"

io_info::io_info( vector<size_t> const & _shp )
{
	this->shape = vector<size_t>();
	this->chunk = vector<size_t>();
	this->compression = MYHDF5_COMPRESSION_NONE;
	this->compression_opts = 0x00;
	this->is_valid = true;
	if ( _shp.size() == 0 ) { //scalar when _shp.size() == 0
		this->shape = vector<size_t>();
	}
	else if ( 1 <= _shp.size() && _shp.size() <= 2 ) { //list
		this->shape = _shp;
	}
	else {
		this->is_valid = false;
	}
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
	cout << this << "\n";
#endif
}


io_info::io_info( vector<size_t> const & _shp, vector<size_t> const & _chk )
{
	this->shape = vector<size_t>();
	this->chunk = vector<size_t>();
	this->compression = MYHDF5_COMPRESSION_NONE;
	this->compression_opts = 0x00;
	this->is_valid = true;
	if ( _shp.size() == 0 ) { //scalar
		this->shape = vector<size_t>();
	}
	else if ( 1 <= _shp.size() && _shp.size() <= 2 ) { //list
		this->shape = _shp;
	}
	else {
		this->is_valid = false;
		return;
	}
	//compression only available with chunked data
	//otherwise uncompressed contiguous design as the fallback
	if ( _shp.size() > 0 && _chk.size() == _shp.size() ) {
		//no chunking for scalars and small vectors
		for( size_t i = 0; i < _shp.size(); i++ ) {
			if ( _shp[i] % _chk[i] != 0 ) {
				this->chunk.push_back( _shp[i] );
			}
			else {
				this->chunk.push_back( _chk[i] );
			}
		}
		if ( this->chunk.size() > 0 ) {
			this->compression = MYHDF5_COMPRESSION_GZIP;
			this->compression_opts = MYHDF5_COMPRESSION_DEFAULT;
		}
	}
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
	cout << this << "\n";
#endif
}


io_info::io_info( vector<size_t> const & _shp, vector<size_t> const & _chk,
			const unsigned char _cmp, const unsigned int _cmp_opts )
{
	this->shape = vector<size_t>();
	this->chunk = vector<size_t>();
	this->compression = MYHDF5_COMPRESSION_NONE;
	this->compression_opts = 0x00;
	this->is_valid = true;
	if ( _shp.size() == 0 ) { //scalar
		this->shape = vector<size_t>();
	}
	else if ( 1 <= _shp.size() && _shp.size() <= 2 ) { //list
		this->shape = _shp;
	}
	else {
		this->is_valid = false;
		return;
	}
	//compression only available with chunked data, otherwise uncompressed contiguous design as the fallback
	if ( _shp.size() > 0 && _chk.size() == _shp.size() ) {
		//no chunking for scalars and small vectors
		for( size_t i = 0; i < _shp.size(); i++ ) {
			if ( _shp[i] % _chk[i] != 0 ) {
				this->chunk.push_back( _shp[i] );
			}
			else {
				this->chunk.push_back( _chk[i] );
			}
		}
		if ( _cmp == MYHDF5_COMPRESSION_GZIP ) {
			this->compression = MYHDF5_COMPRESSION_GZIP;
		}
		if ( _cmp_opts <= 9 ) {
			vector<unsigned char> available_strength = {
					0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09 };
			unsigned char desired_strength = (unsigned char) _cmp_opts;
			for( auto it = available_strength.begin(); it != available_strength.end(); it++ ) {
				if ( desired_strength == *it ) {
					this->compression_opts = *it;
					break;
				}
			}
		}
		else {
			this->compression_opts = 0x00;
		}
	}
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
	cout << this << "\n";
#endif
}


ostream& operator << (ostream& in, io_info const & val)
{
	in << "io_info" << "\n";
	in << "shape.size() " << val.shape.size() << "\n";
	in << "chunk.size() " << val.chunk.size() << "\n";
	if ( val.compression == MYHDF5_COMPRESSION_GZIP )
		in << "compression " << "deflate" << "\n";
	else
		in << "compression " << "none" << "\n";
	in << "compression_opts " << (int) val.compression_opts << "\n";
	in << "is_valid " << (int) val.is_valid << "\n";
	return in;
}


ioAttributes::ioAttributes()
{
	u8 = map<string, unsigned char>();
	i8 = map<string, char>();
	u16 = map<string, unsigned short>();
	i16 = map<string, short>();
	u32 = map<string, unsigned int>();
	i32 = map<string, int>();
	u64 = map<string, unsigned long>();
	i64 = map<string, long>();
	f32 = map<string, float>();
	f64 = map<string, double>();
	str = map<string, string>();

	u8_arr = map<string, vector<unsigned char>>();
	i8_arr = map<string, vector<char>>();
	u16_arr = map<string, vector<unsigned short>>();
	i16_arr = map<string, vector<short>>();
	u32_arr = map<string, vector<unsigned int>>();
	i32_arr = map<string, vector<int>>();
	u64_arr = map<string, vector<unsigned long>>();
	i64_arr = map<string, vector<long>>();
	f32_arr = map<string, vector<float>>();
	f64_arr = map<string, vector<double>>();
	str_arr = map<string, vector<string>>();

	exists = map<string, bool>();
}


ioAttributes::~ioAttributes()
{
}


//scalar attributes
template<typename T> void ioAttributes::add( const string keyword, const T value )
{
	map<string, bool>::iterator thisone = exists.find( keyword );
	if ( thisone == exists.end() ) { //in HDF5, attribute keyword names have to be unique for an object
		if constexpr(std::is_same<T, unsigned char>::value) {
			u8[keyword] = value; exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
			cout << "added attribute u08 " << keyword << " " << (int) value << "\n";
#endif
		}
		else if constexpr(std::is_same<T, char>::value) {
			i8[keyword] = value; exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
			cout << "added attribute i08 " << keyword << " " << (int) value << "\n";
#endif
		}
		else if constexpr(std::is_same<T, unsigned short>::value) {
			u16[keyword] = value; exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
			cout << "added attribute u16 " << keyword << " " << value << "\n";
#endif
		}
		else if constexpr(std::is_same<T, short>::value) {
			i16[keyword] = value; exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
			cout << "added attribute i16 " << keyword << " " << value << "\n";
#endif
		}
		else if constexpr(std::is_same<T, unsigned int>::value) {
			u32[keyword] = value; exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
			cout << "added attribute u32 " << keyword << " " << value << "\n";
#endif
		}
		else if constexpr(std::is_same<T, int>::value) {
			i32[keyword] = value; exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
			cout << "added attribute i32 " << keyword << " " << value << "\n";
#endif
		}
		else if constexpr(std::is_same<T, unsigned long>::value) {
			u64[keyword] = value; exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
			cout << "added attribute u64 " << keyword << " " << value << "\n";
#endif
		}
		else if constexpr(std::is_same<T, long>::value) {
			i64[keyword] = value; exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
			cout << "added attribute i64 " << keyword << " " << value << "\n";
#endif
		}
		else if constexpr(std::is_same<T, float>::value) {
			f32[keyword] = value; exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
			cout << "added attribute f32 " << keyword << " " << value << "\n";
#endif
		}
		else if constexpr(std::is_same<T, double>::value) {
			f64[keyword] = value; exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
			cout << "added attribute f64 " << keyword << " " << value << "\n";
#endif
		}
		else {
			return;
		}
		return;
	}
	cerr << "Attribute named " << keyword << " exists!" << "\n";
	return;
}

template void ioAttributes::add( const string keyword, const unsigned char value );
template void ioAttributes::add( const string keyword, const char value );
template void ioAttributes::add( const string keyword, const unsigned short value );
template void ioAttributes::add( const string keyword, const short value );
template void ioAttributes::add( const string keyword, const unsigned int value );
template void ioAttributes::add( const string keyword, const int value );
template void ioAttributes::add( const string keyword, const unsigned long value );
template void ioAttributes::add( const string keyword, const long value );
template void ioAttributes::add( const string keyword, const float value );
template void ioAttributes::add( const string keyword, const double value );

void ioAttributes::add( const string keyword, const string value )
{
	map<string, bool>::iterator thisone = exists.find( keyword );
	if ( thisone == exists.end() ) { //in HDF5, attribute keyword names have to be unique for an object
		str[keyword] = value; exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "added attribute str " << keyword << " " << value << "\n";
#endif
		return;
	}
	cerr << "Attribute named " << keyword << " exists!" << "\n";
	return;
}


//list attributes
template<typename T> void ioAttributes::add( const string keyword, vector<T> values )
{
	map<string, bool>::iterator thisone = exists.find( keyword );
	if ( thisone == exists.end() ) { //in HDF5, attribute keyword names have to be unique for an object
		if constexpr(std::is_same<T, unsigned char>::value) {
			u8_arr[keyword] = vector<unsigned char>();
			for( auto kt = values.begin(); kt != values.end(); kt++ ) { u8_arr[keyword].push_back( *kt ); }
			exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "added attribute u8_arr " << keyword << " " << u8_arr[keyword].size() << "\n";
#endif
		}
		else if constexpr(std::is_same<T, char>::value) {
			i8_arr[keyword] = vector<char>();
			for( auto kt = values.begin(); kt != values.end(); kt++ ) { i8_arr[keyword].push_back( *kt ); }
			exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "added attribute i8_arr " << keyword << " " << i8_arr[keyword].size() << "\n";
#endif
		}
		else if constexpr(std::is_same<T, unsigned short>::value) {
			u16_arr[keyword] = vector<unsigned short>();
			for( auto kt = values.begin(); kt != values.end(); kt++ ) { u16_arr[keyword].push_back( *kt ); }
			exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "added attribute u16_arr " << keyword << " " << u16_arr[keyword].size() << "\n";
#endif
		}
		else if constexpr(std::is_same<T, short>::value) {
			i16_arr[keyword] = vector<short>();
			for( auto kt = values.begin(); kt != values.end(); kt++ ) { i16_arr[keyword].push_back( *kt ); }
			exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "added attribute i16_arr " << keyword << " " << i16_arr[keyword].size() << "\n";
#endif
		}
		else if constexpr(std::is_same<T, unsigned int>::value) {
			u32_arr[keyword] = vector<unsigned int>();
			for( auto kt = values.begin(); kt != values.end(); kt++ ) { u32_arr[keyword].push_back( *kt ); }
			exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "added attribute u32_arr " << keyword << " " << u32_arr[keyword].size() << "\n";
#endif
		}
		else if constexpr(std::is_same<T, int>::value) {
			i32_arr[keyword] = vector<int>();
			for( auto kt = values.begin(); kt != values.end(); kt++ ) { i32_arr[keyword].push_back( *kt ); }
			exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "added attribute i32_arr " << keyword << " " << i32_arr[keyword].size() << "\n";
#endif
		}
		else if constexpr(std::is_same<T, unsigned long>::value) {
			u64_arr[keyword] = vector<unsigned long>();
			for( auto kt = values.begin(); kt != values.end(); kt++ ) { u64_arr[keyword].push_back( *kt ); }
			exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "added attribute u64_arr " << keyword << " " << u64_arr[keyword].size() << "\n";
#endif
		}
		else if constexpr(std::is_same<T, long>::value) {
			i64_arr[keyword] = vector<long>();
			for( auto kt = values.begin(); kt != values.end(); kt++ ) { i64_arr[keyword].push_back( *kt ); }
			exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "added attribute i64_arr " << keyword << " " << i64_arr[keyword].size() << "\n";
#endif
		}
		else if constexpr(std::is_same<T, float>::value) {
			f32_arr[keyword] = vector<float>();
			for( auto kt = values.begin(); kt != values.end(); kt++ ) { f32_arr[keyword].push_back( *kt ); }
			exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "added attribute f32_arr " << keyword << " " << f32_arr[keyword].size() << "\n";
#endif
		}
		else if constexpr(std::is_same<T, double>::value) {
			f64_arr[keyword] = vector<double>();
			for( auto kt = values.begin(); kt != values.end(); kt++ ) { f64_arr[keyword].push_back( *kt ); }
			exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "added attribute f64_arr " << keyword << " " << f64_arr[keyword].size() << "\n";
#endif
		}
		else {
			return;
		}
		return;
	}
	cerr << "Attribute named " << keyword << " exists!" << "\n";
	return;
}

template void ioAttributes::add( const string keyword, vector<unsigned char> values );
template void ioAttributes::add( const string keyword, vector<char> values );
template void ioAttributes::add( const string keyword, vector<unsigned short> values );
template void ioAttributes::add( const string keyword, vector<short> values );
template void ioAttributes::add( const string keyword, vector<unsigned int> values );
template void ioAttributes::add( const string keyword, vector<int> values );
template void ioAttributes::add( const string keyword, vector<unsigned long> values );
template void ioAttributes::add( const string keyword, vector<long> values );
template void ioAttributes::add( const string keyword, vector<float> values );
template void ioAttributes::add( const string keyword, vector<double> values );

void ioAttributes::add( const string keyword, vector<string> values )
{
	map<string, bool>::iterator thisone = exists.find( keyword );
	if ( thisone == exists.end() ) { //in HDF5, attribute keyword names have to be unique for an object
		str_arr[keyword] = vector<string>();
		for( auto kt = values.begin(); kt != values.end(); kt++ ) { str_arr[keyword].push_back( *kt ); }
		exists[keyword] = true;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "added attribute str_arr " << keyword << " " << str_arr[keyword].size() << "\n";
#endif
		return;
	}
	cerr << "Attribute named " << keyword << " exists!" << "\n";
	return;
}

//##MK::this will not work as intended due to rvalue issues and compiler cannot deduce the rvalue type

/*
template<typename T> T ioAttributes::query( const string keyword )
{
	map<string, bool>::iterator thisone = exists.find( keyword );
	if ( thisone != exists.end() ) {
		if constexpr(std::is_same<T, unsigned char>::value) {
			map<string, unsigned char>::iterator jt = u8.find( keyword );
			if ( jt != u8.end() )
				return (unsigned char) jt->second;
		}
		else if constexpr(std::is_same<T, char>::value) {
			map<string, char>::iterator jt = i8.find( keyword );
			if ( jt != i8.end() )
				return (char) jt->second;
		}
		else if constexpr(std::is_same<T, unsigned short>::value) {
			map<string, unsigned short>::iterator jt = u16.find( keyword );
			if ( jt != u16.end() )
				return (unsigned short) jt->second;
		}
		else if constexpr(std::is_same<T, short>::value) {
			map<string, short>::iterator jt = i16.find( keyword );
			if ( jt != i16.end() )
				return (short) jt->second;
		}
		else if constexpr(std::is_same<T, unsigned int>::value) {
			map<string, unsigned int>::iterator jt = u32.find( keyword );
			if ( jt != u32.end() )
				return (unsigned int) jt->second;
		}
		else if constexpr(std::is_same<T, int>::value) {
			map<string, int>::iterator jt = i32.find( keyword );
			if ( jt != i32.end() )
				return (int) jt->second;
		}
		else if constexpr(std::is_same<T, unsigned long>::value) {
			map<string, unsigned long>::iterator jt = u64.find( keyword );
			if ( jt != u64.end() )
				return (unsigned long) jt->second;
		}
		else if constexpr(std::is_same<T, long>::value) {
			map<string, long>::iterator jt = i64.find( keyword );
			if ( jt != i64.end() )
				return (long) jt->second;
		}
		else if constexpr(std::is_same<T, float>::value) {
			map<string, float>::iterator jt = f32.find( keyword );
			if ( jt != f32.end() )
				return (float) jt->second;
		}
		else if constexpr(std::is_same<T, double>::value) {
			map<string, double>::iterator jt = f64.find( keyword );
			if ( jt != f64.end() )
				return (double) jt->second;
		}
		else {
			return (T) 0;
		}
		return (T) 0;
	}
	return (T) 0;
}

template unsigned char ioAttributes::query( const string keyword );
template char ioAttributes::query( const string keyword );
template unsigned short ioAttributes::query( const string keyword );
template short ioAttributes::query( const string keyword );
template unsigned int ioAttributes::query( const string keyword );
template int ioAttributes::query( const string keyword );
template unsigned long ioAttributes::query( const string keyword );
template long ioAttributes::query( const string keyword );
template float ioAttributes::query( const string keyword );
template double ioAttributes::query( const string keyword );

string ioAttributes::query( const string keyword )
{
	map<string, string>::iterator jt = str.find( keyword );
	if ( jt != str.end() ) {
		return jt->second;
	}
	return "";
}


template<typename T> vector<T> ioAttributes::query( const string keyword )
{
	map<string, bool>::iterator thisone = exists.find( keyword );
	if ( thisone != exists.end() ) {
		if constexpr(std::is_same<T, unsigned char>::value) {
			map<string, vector<unsigned char>>::iterator jt = u8_arr.find( keyword );
			if ( jt != u8_arr.end() )
				return jt->second;
		}
		else if constexpr(std::is_same<T, char>::value) {
			map<string, vector<char>>::iterator jt = i8_arr.find( keyword );
			if ( jt != i8_arr.end() )
				return jt->second;
		}
		else if constexpr(std::is_same<T, unsigned short>::value) {
			map<string, vector<unsigned short>>::iterator jt = u16_arr.find( keyword );
			if ( jt != u16_arr.end() )
				return jt->second;
		}
		else if constexpr(std::is_same<T, short>::value) {
			map<string, vector<short>>::iterator jt = i16_arr.find( keyword );
			if ( jt != i16_arr.end() )
				return jt->second;
		}
		else if constexpr(std::is_same<T, unsigned int>::value) {
			map<string, vector<unsigned int>>::iterator jt = u32_arr.find( keyword );
			if ( jt != u32_arr.end() )
				return jt->second;
		}
		else if constexpr(std::is_same<T, int>::value) {
			map<string, vector<int>>::iterator jt = i32_arr.find( keyword );
			if ( jt != i32_arr.end() )
				return jt->second;
		}
		else if constexpr(std::is_same<T, unsigned long>::value) {
			map<string, vector<unsigned long>>::iterator jt = u64_arr.find( keyword );
			if ( jt != u64_arr.end() )
				return jt->second;
		}
		else if constexpr(std::is_same<T, long>::value) {
			map<string, vector<long>>::iterator jt = i64_arr.find( keyword );
			if ( jt != i64_arr.end() )
				return jt->second;
		}
		else if constexpr(std::is_same<T, float>::value) {
			map<string, vector<float>>::iterator jt = f32_arr.find( keyword );
			if ( jt != f32_arr.end() )
				return jt->second;
		}
		else if constexpr(std::is_same<T, double>::value) {
			map<string, vector<double>>::iterator jt = f64_arr.find( keyword );
			if ( jt != f64_arr.end() )
				return jt->second;
		}
		else {
			return vector<T>();
		}
		return vector<T>();
	}
	return vector<T>();
}

template vector<unsigned char> ioAttributes::query( const string keyword );
template vector<char> ioAttributes::query( const string keyword );
template vector<unsigned short> ioAttributes::query( const string keyword );
template vector<short> ioAttributes::query( const string keyword );
template vector<unsigned int> ioAttributes::query( const string keyword );
template vector<int> ioAttributes::query( const string keyword );
template vector<unsigned long> ioAttributes::query( const string keyword );
template vector<long> ioAttributes::query( const string keyword );
template vector<float> ioAttributes::query( const string keyword );
template vector<double> ioAttributes::query( const string keyword );

vector<string> ioAttributes::query( const string keyword )
{
	map<string, vector<string>>::iterator jt = str_arr.find( keyword );
	if ( jt != str_arr.end() ) {
		return jt->second;
	}
	return vector<string>();
}
*/


void ioAttributes::report()
{
	cout << "Named attributes in the dictionary" << "\n";
	for( auto it = exists.begin(); it != exists.end(); it++ ) {
		if ( it->second == true ) {
			cout << it->first << "\n";
		}
	}

	cout << "Named scalar attributes in the dictionary" << "\n";
	cout << "Named u08 scalar attribute in the dictionary" << "\n";
	for( auto it = u8.begin(); it != u8.end(); it++ ) {
		cout << it->first << ", " << (int) it->second << "\n";
	}
	cout << "Named i08 scalar attribute in the dictionary" << "\n";
	for( auto it = i8.begin(); it != i8.end(); it++ ) {
		cout << it->first << ", " << (int) it->second << "\n";
	}
	cout << "Named u16 scalar attribute in the dictionary" << "\n";
	for( auto it = u16.begin(); it != u16.end(); it++ ) {
		cout << it->first << ", " << it->second << "\n";
	}
	cout << "Named i16 scalar attribute in the dictionary" << "\n";
	for( auto it = i16.begin(); it != i16.end(); it++ ) {
		cout << it->first << ", " << it->second << "\n";
	}
	cout << "Named u32 scalar attribute in the dictionary" << "\n";
	for( auto it = u32.begin(); it != u32.end(); it++ ) {
		cout << it->first << ", " << it->second << "\n";
	}
	cout << "Named i32 scalar attribute in the dictionary" << "\n";
	for( auto it = i32.begin(); it != i32.end(); it++ ) {
		cout << it->first << ", " << it->second << "\n";
	}
	cout << "Named u64 scalar attribute in the dictionary" << "\n";
	for( auto it = u64.begin(); it != u64.end(); it++ ) {
		cout << it->first << ", " << it->second << "\n";
	}
	cout << "Named i64 scalar attribute in the dictionary" << "\n";
	for( auto it = i64.begin(); it != i64.end(); it++ ) {
		cout << it->first << ", " << it->second << "\n";
	}
	cout << "Named f32 scalar attribute in the dictionary" << "\n";
	for( auto it = f32.begin(); it != f32.end(); it++ ) {
		cout << it->first << ", " << it->second << "\n";
	}
	cout << "Named f64 scalar attribute in the dictionary" << "\n";
	for( auto it = f64.begin(); it != f64.end(); it++ ) {
		cout << it->first << ", " << it->second << "\n";
	}
	cout << "Named string scalar attribute in the dictionary" << "\n";
	for( auto it = str.begin(); it != str.end(); it++ ) {
		cout << it->first << ", " << it->second << "\n";
	}

	cout << "Named list attributes in the dictionary" << "\n";
	cout << "Named u08 list attribute in the dictionary" << "\n";
	for( auto it = u8_arr.begin(); it != u8_arr.end(); it++ ) {
		cout << it->first << ", " << it->second.size() << "\n";
	}
	cout << "Named i08 list attribute in the dictionary" << "\n";
	for( auto it = i8_arr.begin(); it != i8_arr.end(); it++ ) {
		cout << it->first << ", " << it->second.size() << "\n";
	}
	cout << "Named u16 list attribute in the dictionary" << "\n";
	for( auto it = u16_arr.begin(); it != u16_arr.end(); it++ ) {
		cout << it->first << ", " << it->second.size() << "\n";
	}
	cout << "Named i16 list attribute in the dictionary" << "\n";
	for( auto it = i16_arr.begin(); it != i16_arr.end(); it++ ) {
		cout << it->first << ", " << it->second.size() << "\n";
	}
	cout << "Named u32 list attribute in the dictionary" << "\n";
	for( auto it = u32_arr.begin(); it != u32_arr.end(); it++ ) {
		cout << it->first << ", " << it->second.size() << "\n";
	}
	cout << "Named i32 list attribute in the dictionary" << "\n";
	for( auto it = i32_arr.begin(); it != i32_arr.end(); it++ ) {
		cout << it->first << ", " << it->second.size() << "\n";
	}
	cout << "Named u64 list attribute in the dictionary" << "\n";
	for( auto it = u64_arr.begin(); it != u64_arr.end(); it++ ) {
		cout << it->first << ", " << it->second.size() << "\n";
	}
	cout << "Named i64 list attribute in the dictionary" << "\n";
	for( auto it = i64_arr.begin(); it != i64_arr.end(); it++ ) {
		cout << it->first << ", " << it->second.size() << "\n";
	}
	cout << "Named f32 list attribute in the dictionary" << "\n";
	for( auto it = f32_arr.begin(); it != f32_arr.end(); it++ ) {
		cout << it->first << ", " << it->second.size() << "\n";
	}
	cout << "Named f64 list attribute in the dictionary" << "\n";
	for( auto it = f64_arr.begin(); it != f64_arr.end(); it++ ) {
		cout << it->first << ", " << it->second.size() << "\n";
	}
	cout << "Named string list attribute in the dictionary" << "\n";
	for( auto it = str_arr.begin(); it != str_arr.end(); it++ ) {
		cout << it->first << ", " << it->second.size() << "\n";
	}
}


HdfFiveSeqHdl::HdfFiveSeqHdl()
{
	h5resultsfn = "PARAPROBE.Test.nxs";
	fileid = 0; //H5I_INVALID_HID;
	mtypid = 0; //H5I_INVALID_HID;
	dspcid = 0; //H5I_INVALID_HID;
	plistid = 0; //H5I_INVALID_HID;
	dsetid = 0; //H5I_INVALID_HID;
	dtypid = 0; //H5I_INVALID_HID;
	//attr_buf = ioAttributes();
}


HdfFiveSeqHdl::HdfFiveSeqHdl( const string h5fn )
{
	h5resultsfn = h5fn;
	fileid = 0; //H5I_INVALID_HID;
	mtypid = 0; //H5I_INVALID_HID;
	dspcid = 0; //H5I_INVALID_HID;
	plistid = 0; //H5I_INVALID_HID;
	dsetid = 0; //H5I_INVALID_HID;
	dtypid = 0; //H5I_INVALID_HID;
	//attr_buf = ioAttributes();
}


HdfFiveSeqHdl::~HdfFiveSeqHdl()
{
}


//low-level functions which must not be called alone but only from within a function
int HdfFiveSeqHdl::nexus_open( const unsigned flags )
{
	//open an existent file in a mode given by flags
	fileid = H5Fopen( h5resultsfn.c_str(), flags, H5P_DEFAULT);
	if ( fileid != H5I_INVALID_HID ) {
		//nop
	}
	else {
		cerr << "Opening " << h5resultsfn << " failed !" << "\n";
	}
	return fileid;
}


int HdfFiveSeqHdl::nexus_attribute_close( const int previous_status )
{
	//low-level function to release dangling object handles correctly to not leave the handle manager polluted after an access to an object
	int current_status = previous_status;
	htri_t state = 0;

	state = H5Iis_valid(adtypid);
	if ( state > 0 ) {
		if ( H5Tclose(adtypid) >= 0 ) {
		}
		else {
			cerr << "H5Tclose error on " << h5resultsfn << "\n";
			current_status = MYHDF5_FAILED;
		}
	}
	else if ( state == 0 ) {
	}
	else {
		cerr << "H5Iis_valid error on " << h5resultsfn << "\n";
		current_status = MYHDF5_FAILED;
	}

	state = H5Iis_valid(aspcid);
	if ( state > 0 ) { //identifier is still valid and thus should be closed
		if ( H5Sclose(aspcid) >= 0 ) {
		}
		else {
			cerr << "H5Sclose error on " << h5resultsfn << "\n";
			current_status = MYHDF5_FAILED;
		}
	}
	else if ( state == 0 ) { //identifer is no longer valid or call was errorneous, identifier must not be closed again
	}
	else {
		cerr << "H5Iis_valid error on " << h5resultsfn << "\n";
		current_status = MYHDF5_FAILED;
	}

	state = H5Iis_valid(attrid);
	if ( state > 0 ) {
		if ( H5Aclose(attrid) >= 0 ) {
		}
		else {
			cerr << "H5Aclose error on " << h5resultsfn << "\n";
			current_status = MYHDF5_FAILED;
		}
	}
	else if ( state == 0 ) {
	}
	else {
		cerr << "H5Iis_valid error on " << h5resultsfn << "\n";
		current_status = MYHDF5_FAILED;
	}

#ifdef MYHDF5_VERBOSE_LEVEL_TWO
	cout << "nexus_attribute_close " << current_status << "\n";
#endif
	return current_status;
}


int HdfFiveSeqHdl::nexus_close( const int previous_status )
{
	//low-level function to release dangling object handles correctly to not leave the handle manager polluted after an access to an object
	int current_status = previous_status;
	htri_t state = 0;

	state = H5Iis_valid(dsetid);
	if ( state > 0 ) {
		if ( H5Dclose(dsetid) >= 0 ) {
		}
		else {
			cerr << "H5Dclose error on " << h5resultsfn << "\n";
			current_status = MYHDF5_DSETCLOSE_FAILED;
		}
	}
	else if ( state == 0 ) {
	}
	else {
		cerr << "H5Iis_valid error on " << h5resultsfn << "\n";
		current_status = MYHDF5_FAILED;
	}

	state = H5Iis_valid(dtypid);
	if ( state > 0 ) { //identifier is still valid and thus should be closed
		if ( H5Tclose(dtypid) >= 0 ) {
		}
		else {
			cerr << "H5Tclose error on " << h5resultsfn << "\n";
			current_status = MYHDF5_MSPCCLOSE_FAILED;
		}
	}
	else if ( state == 0 ) { //identifer is no longer valid or call was errorneous, identifier must not be closed again
	}
	else {
		cerr << "H5Iis_valid error on " << h5resultsfn << "\n";
		current_status = MYHDF5_FAILED;
	}

	state = H5Iis_valid(dspcid);
	if ( state > 0 ) {
		if ( H5Sclose(dspcid) >= 0 ) {
		}
		else {
			cerr << "H5Sclose error on " << h5resultsfn << "\n";
			current_status = MYHDF5_DSPCCLOSE_FAILED;
		}
	}
	else if ( state == 0 ) {
	}
	else {
		cerr << "H5Iis_valid error on " << h5resultsfn << "\n";
		current_status = MYHDF5_FAILED;
	}

	state = H5Iis_valid(mtypid);
	if ( state > 0 ) {
		if ( H5Tclose(mtypid) >= 0 ) {
		}
		else {
			cerr << "H5Tclose error on " << h5resultsfn << "\n";
			current_status = MYHDF5_FAILED;
		}
	}
	else if ( state == 0 ) {
	}
	else {
		cerr << "H5Iis_valid error on " << h5resultsfn << "\n";
		current_status = MYHDF5_FAILED;
	}

	state = H5Iis_valid(plistid);
	if ( state > 0 ) {
		if ( H5Pclose(plistid) >= 0 ) {
		}
		else {
			cerr << "H5Pclose error on " << h5resultsfn << "\n";
			current_status = MYHDF5_PLSTCLOSE_FAILED;
		}
	}
	else if ( state == 0 ) {
	}
	else {
		cerr << "H5Iis_valid error on " << h5resultsfn << "\n";
		current_status = MYHDF5_FAILED;
	}

	state = H5Iis_valid(fileid);
	if ( state > 0 ) {
		if ( H5Fclose( fileid ) >= 0 ) {
		}
		else {
			cerr << "H5Fclose error on " << h5resultsfn << "\n";
			current_status = MYHDF5_FCLOSE_FAILED;
		}
	}
	else if ( state == 0 ) {
	}
	else {
		cerr << "H5Iis_valid error on " << h5resultsfn << "\n";
		current_status = MYHDF5_FAILED;
	}
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
	cout << "nexus_close " << current_status << "\n";
#endif
	return current_status;
}


//high-level function
int HdfFiveSeqHdl::nexus_create()
{
	 //create an HDF5 file truncating ALL eventually existent content when a file with the same name already exists
	fileid = H5Fcreate( h5resultsfn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
	if ( fileid != H5I_INVALID_HID ) {
		if ( H5Fclose(fileid) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
			cout << h5resultsfn << " created successfully" << "\n";
#endif
			return MYHDF5_SUCCESS;
		}
		else {
			cerr << "H5Fclose failed on " << h5resultsfn << "\n";
			return MYHDF5_FCLOSE_FAILED;
		}
	}
	//dont use an invalid fileid to attempt closing the file
	return MYHDF5_FILEACCESS_FAILED;
}


bool HdfFiveSeqHdl::nexus_path_exists( const string h5_absolute_path )
{
	//check if a given path in an HDF5 file exists
	//mainly relevant to recursively instantiate a chain of nested groups to create a dataset deep down in the hierarchy
	if ( nexus_open( H5F_ACC_RDONLY ) == H5I_INVALID_HID ) {
		cerr << "Opening file during link exists returned invalid fileid!" << "\n";
		return false;
	}
	//now the file is open otherwise we wouldnt have gotten here

	//check if the path in the open file pointed to by fileid exists
	//cout << "absolute_path __" << h5_absolute_path << "__" << "\n";
	vector<string> level_by_level = split_absolute_path( h5_absolute_path );

	//start at the file, root
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
	cout << "Checking if each group exists level-by-level..." << "\n";
#endif
	hid_t curr_loc_id = fileid;
	string curr_group = "";
	string current_path = "";
	string absolute_path = "/" + h5_absolute_path;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
	cout << "absolute_path __" << absolute_path << "__" << "\n";
#endif

	//##MK::https://gist.github.com/jzrake/3025642
	for( size_t i = 0; i < level_by_level.size(); i++ ) {
		curr_group = level_by_level[i];
		current_path += "/" + level_by_level[i];
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "Currently visiting " << curr_group << "\n";
		cout << "current_path __" << current_path << "__" << "\n";
#endif
		htri_t status = H5Lexists( curr_loc_id, curr_group.c_str(), H5P_DEFAULT );
		//htri_t three-valued outcome, exists (>0), does not exist (==0), error (<0)
		if( status > 0 ) { //exists
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
			cout << curr_group << " --> exists" << "\n";
#endif
			//only really open the group if we have not yet arrived at location named h5_absolute_path
			//take two examples
			///entry/group, this case will work
			///entry/dset, this case will throw because yes /entry/dset may exists as reported by H5Lexists but H5Gopen cannot be applied to it dset because it is
			if ( current_path.compare( absolute_path ) != 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
				cout << "Not yet arrived at absolute path" << "\n";
#endif
				hid_t next_loc_id = H5Gopen( curr_loc_id, curr_group.c_str(), H5P_DEFAULT );
				if ( next_loc_id == H5I_INVALID_HID ) { //cant traverse further, need to close file...
					if ( H5Fclose( fileid ) >= 0 ) {
						cerr << "Successful H5Fclose after failed H5Gopen on " << curr_group << "\n";
						return false;
					}
					cerr << "Failed H5Fclose after failed H5Gopen on " << curr_group << "\n";
					return false;
				}
				else { //traverse further
					if ( curr_loc_id != fileid ) {
						if ( H5Gclose( curr_loc_id ) >= 0 ) {
							//nop
						}
						else {
							if ( H5Fclose( fileid ) >= 0 ) {
								cerr << "Successful H5Fclose after a failed H5Gclose while preparing walk after " << curr_group << "\n";
								return  false;
							}
							cerr << "Failed H5Fclose after a failed H5Gclose while preparing walk after " << curr_group << "\n";
							return false;
						}
					}
					curr_loc_id = next_loc_id;
					continue;
				}
			} //not yet arrived at the target absolute path
			else {
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
				cout << "Arrived at absolute path, stopping traversal" << "\n";
#endif
				if ( H5Fclose( fileid ) >= 0 ) {
					return true;
				}
				cerr << "Failed H5Fclose after a successful H5Gclose while stopping traversal " << curr_group << "\n";
				return false;
			}
		}
		else { //status <= 0
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
			cout << curr_group << " --> does not exist" << "\n";
#endif
			if ( curr_loc_id != fileid ) {
				if ( H5Gclose( curr_loc_id ) >= 0 ) {
					if ( H5Fclose( fileid ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
						cout << "Successful H5Fclose after a successful H5Gclose when link did not exist " << curr_group << "\n";
#endif
						return false;
					}
					cerr << "Failed H5Fclose after a successful H5Gclose when link did not exist " << curr_group << "\n";
					return false;
				}
				else {
					if ( H5Fclose( fileid ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
						cout << "Successful H5Fclose after a failed H5Gclose when link did not exist " << curr_group << "\n";
#endif
						return false;
					}
					cerr << "Failed H5Fclose after a failed H5Gclose when link did not exist " << curr_group << "\n";
					return false;
				}
			}
			else { //no need to close a group because we are still at fileid, i.e. root
				if ( H5Fclose( fileid ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
					cout << "Successful H5Fclose when link did not exist " << curr_group << "\n";
#endif
					return false;
				}
				cerr << "Failed H5Fclose when link did not exist " << curr_group << "\n";
				return false;
			}
		}
	}

	//success is achieved only when you have not previous run into a failure
	//such as having to stop the walk because of a missing path
	if ( H5Gclose( curr_loc_id ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
		cout << "Successful H5Gclose during cleanup" << "\n";
#endif
		if ( H5Fclose( fileid) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
			cout << "Successful H5Fclose during cleanup" << "\n";
#endif
			return true;
		}
		else {
			cerr << "Failed H5Fclose during cleanup" << "\n";
			return false;
		}
	}
	else {
		cerr << "Failed H5Gclose during cleanup" << "\n";
		if ( H5Fclose( fileid ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
			cout << "Successful H5Fclose after failed H5Gclose during cleanup" << "\n";
#endif
			return false;
		}
		else {
			cerr << "Failed H5Fclose after failed H5Gclose during cleanup" << "\n";
			return false;
		}
	}
}


/*
int HdfFiveSeqHdl::nexus_read_attributes( const string dsnm, ioAttributes & attrs )
{
	//MK::this function must only be called from within a nexus_read template instance as it
	//assumes/reuses open(ed) resources on an existent file and dataset to reduce costs
	//for checks and unnecessary interaction with the HDF5 library !

	//assume fileid, dset are set to valid resources on an existent dataset entry dsnm in an open file
	//this function extracts all scalar attributes for elementary types
	//https://fossies.org/linux/hdf5/examples/h5_attribute.c

	//clear the class-internal attribute buffer
	//attr_buf = ioAttributes();
	attrs = ioAttributes();

	if ( nexus_open( H5F_ACC_RDONLY ) != H5I_INVALID_HID ) {
		//if ( nexus_path_exists( dsnm ) == true ) { //H5Gget_objinfo( fileid, dsnm.c_str(), 0, NULL) ) >= 0 ) {
			dsetid = H5Dopen2( fileid, dsnm.c_str(), H5P_DEFAULT );
			if ( dsetid != H5I_INVALID_HID ) {
				H5O_info2_t oinfo;
				if ( H5Oget_info3(dsetid, &oinfo, H5O_INFO_NUM_ATTRS) >= 0 ) {

					for ( unsigned i = 0; i < (unsigned) oinfo.num_attrs; i++ ) {
						hid_t attrid = H5Aopen_by_idx( dsetid, ".", H5_INDEX_CRT_ORDER, H5_ITER_INC, (hsize_t) i, H5P_DEFAULT, H5P_DEFAULT );
						if ( attrid >= 0 ) {
							hid_t dtype = H5Aget_type(attrid);
							if ( dtype != H5I_INVALID_HID ) {
								//parse details of the attribute to identify which datatype was used
								char* rbuf = new char[U16MX];
								ssize_t nbytes = H5Aget_name(attrid, (size_t) U16MX, rbuf );
								H5T_class_t dtype_class_id = H5Tget_class(dtype);
								size_t dtype_size = H5Tget_size(dtype);
								H5T_sign_t dtype_sign = H5Tget_sign(dtype);

								stringstream linestream;
								for( ssize_t j = 0; j < nbytes; j++  ) {
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
									cout << "-->__" << rbuf[j] << "__" << "\n";
#endif
									linestream << rbuf[j];
								}
								delete [] rbuf; rbuf = NULL;
								string keyword = linestream.str();

#ifdef MYHDF5_VERBOSE_LEVEL_TWO
								cout << "attribute name " << keyword << "\n";
								cout << "attribute dtype class id " << dtype_class_id << "\n";
								cout << "attribute dtype class size bytes " << dtype_size << "\n";
								cout << "attribute dtype sgn " << dtype_sign << "\n";
#endif

								//extract interpretable scalar or string attribute values of specific type
								//assuming that paraprobe uses the H5T specific little endian variants
								switch(dtype_class_id)
								{
									case H5T_INTEGER:
									{
										switch(dtype_size)
										{
											case 1:
											{
												if ( dtype_sign == H5T_SGN_NONE ) { //H5T_STD_U8LE
													unsigned char value = 0x00;
													if ( H5Aread(attrid, dtype, &value ) >= 0 ) {
														attrs.add( keyword, value );
													}
												}
												else if ( dtype_sign == H5T_SGN_2 ) { //H5T_STD_I8LE
													char value = 0x00;
													if ( H5Aread(attrid, dtype, &value ) >= 0 ) {
														attrs.add( keyword, value );
													}
												}
												else {
													//H5T_SGN_ERROR
												}
												break;
											}
											case 2:
											{
												if ( dtype_sign == H5T_SGN_NONE ) { //H5T_STD_U16LE
													unsigned short value = 0;
													if ( H5Aread(attrid, dtype, &value ) >= 0 ) {
														attrs.add( keyword, value );
													}
												}
												else if ( dtype_sign == H5T_SGN_2 ) { //H5T_STD_I16LE
													short value = 0;
													if ( H5Aread(attrid, dtype, &value ) >= 0 ) {
														attrs.add( keyword, value );
													}
												}
												else {
													//H5T_SGN_ERROR
												}
												break;
											}
											case 4:
											{
												if ( dtype_sign == H5T_SGN_NONE ) { //H5T_STD_U32LE
													unsigned int value = 0;
													if ( H5Aread(attrid, dtype, &value ) >= 0 ) {
														attrs.add( keyword, value );
													}
												}
												else if ( dtype_sign == H5T_SGN_2 ) { //H5T_STD_I32LE
													int value = 0;
													if ( H5Aread(attrid, dtype, &value ) >= 0 ) {
														attrs.add( keyword, value );
													}
												}
												else {
													//H5T_SGN_ERROR
												}
												break;
											}
											case 8:
											{
												if ( dtype_sign == H5T_SGN_NONE ) { //H5T_STD_U64LE
													unsigned long value = 0;
													if ( H5Aread(attrid, dtype, &value ) >= 0 ) {
														attrs.add( keyword, value );
													}
												}
												else if ( dtype_sign == H5T_SGN_2 ) { //H5T_STD_I64LE
													long value = 0;
													if ( H5Aread(attrid, dtype, &value ) >= 0 ) {
														attrs.add( keyword, value );
													}
												}
												else {
													//H5T_SGN_ERROR
												}
												break;
											}
											default:
											{
												break;
											}
										}
										break;
									}
									case H5T_FLOAT:
									{
										if ( dtype_size == 4 ) {
											float value = MYZERO;
											if ( H5Aread(attrid, dtype, &value ) >= 0 ) {
												attrs.add( keyword, value );
											}
										}
										else if ( dtype_size == 8 ) {
											double value = MYZERO;
											if ( H5Aread(attrid, dtype, &value ) >= 0 ) {
												attrs.add( keyword, value );
											}
										}
										else {
											//nope
										}
										break;
									}
									case H5T_STRING:
									{
										hid_t atype_mem = H5Tget_native_type(dtype, H5T_DIR_ASCEND);
										if ( atype_mem != H5I_INVALID_HID ) {
											char* rbuf = new char[U16MX];
											if ( H5Aread(attrid, atype_mem, rbuf ) >= 0 ) {
												stringstream tmpstream;
												for( size_t j = 0; j < dtype_size; j++  ) {
													if ( rbuf[j] != 0x00 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
														cout << "------>__" << rbuf[j] << "__" << "\n";
#endif
														tmpstream << rbuf[j];
													}
													else {
														break;
													}
												}
												string value = tmpstream.str();
												delete [] rbuf; rbuf = NULL;
												if ( value.length() > 0 ) {
													attrs.add( keyword, value );
												}
											}
											delete [] rbuf; rbuf = NULL;
										}
										break;
									}
									default:
									{
										break;
									}
								}
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
								cout << "Reading " << dsnm << " attribute keyword " << keyword << " success" << "\n";
#endif
							}
							else {
								if ( H5Aclose(attrid) >= 0 ) {
									cerr << "Reading " << dsnm << " attribute failed, dtype invalid will continue with next attribute !" << "\n";
									continue;
								}
								else {
									cerr << "Reading " << dsnm << " attribute failed, unable to close opened attribute, stopping reading further attributes !" << "\n";
									break;
								}
							}
						}
						else {
							cerr << "Reading " << dsnm << " attribute failed, attribute not openable, stopping reading further attributes !" << "\n";
							break;
						}
					} //inspect the next attribute
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
					cout << "Reading " << dsnm << " all attributes visited" << "\n";
#endif
				}
				if ( H5Dclose( dsetid ) >= 0 ) {
					if ( H5Fclose( fileid ) >= 0 ) {
						return MYHDF5_SUCCESS;
					}
					cerr << "H5Fclose error on " << h5resultsfn << "\n";
					return MYHDF5_FCLOSE_FAILED;
				}
				cerr << "H5Dclose error on " << h5resultsfn << "\n";
				if ( H5Fclose( fileid ) >= 0 ) {
					return MYHDF5_DSETCLOSE_FAILED;
				}
				cerr << "H5Fclose error on " << h5resultsfn << "\n";
				return MYHDF5_FCLOSE_FAILED;
			}
			cerr << "H5Dopen error on " << h5resultsfn << "\n";
			if ( H5Fclose( fileid ) >= 0 ) {
				return MYHDF5_DSETOPEN_FAILED;
			}
			cerr << "H5Fclose error on " << h5resultsfn << "\n";
			return MYHDF5_FCLOSE_FAILED;
		//}
	}
	return MYHDF5_NOT_EXPECTED;

	//// get attribute info using iteration function
	//if ( H5Aiterate2(dsetid, H5_INDEX_NAME, H5_ITER_INC, NULL, attr_info, NULL) >= 0 ) {
	//	//copy over the attributes from the buffer and clear that buffer
	//	attrs = attr_buf;
	//	attr_buf = ioAttributes();
	//	return MYHDF5_SUCCESS;
	//}
	//return MYHDF5_ATTRREAD_FAILED;
}
*/

/*
//deprecated:
template<typename T> int HdfFiveSeqHdl::nexus_read( const string dsnm, T & retval )
{
	vector<T> tmp = vector<T>();
	if ( nexus_read( dsnm, tmp ) == MYHDF5_SUCCESS ) {
		if ( tmp.size() >= 1 ) {
			retval = tmp[0];
			return MYHDF5_SUCCESS;
		}
		cerr << "Reading " << dsnm << " scalar, failed, tmp.size() == 0 !" << "\n";
		return MYHDF5_NOT_EXPECTED;
	}
	cerr << "Reading " << dsnm << " scalar, failed, output left unchanged !" << "\n";
	return MYHDF5_DSETREAD_FAILED;
}
*/

template <typename T> int HdfFiveSeqHdl::nexus_read( const string dsnm, T & retval )
{
	//open file in read-only mode to not modified the file last access time stamp
	//reads list_single and list_long of base C/C++ data types, for vector<string> there is an own function
	dtypid = H5I_INVALID_HID;
	string msg = "";
	if constexpr(std::is_same<T, unsigned char>::value) { dtypid = H5T_STD_U8LE; msg = "u8"; }
	else if constexpr(std::is_same<T, char>::value) { dtypid = H5T_STD_I8LE; msg = "i8"; }
	else if constexpr(std::is_same<T, unsigned short>::value) { dtypid = H5T_STD_U16LE; msg = "u16"; }
	else if constexpr(std::is_same<T, short>::value) { dtypid = H5T_STD_I16LE; msg = "i16"; }
	else if constexpr(std::is_same<T, unsigned int>::value) { dtypid = H5T_STD_U32LE; msg = "u32"; }
	else if constexpr(std::is_same<T, int>::value) { dtypid = H5T_STD_I32LE; msg = "i32"; }
	else if constexpr(std::is_same<T, unsigned long>::value) { dtypid = H5T_STD_U64LE; msg = "u64"; }
	else if constexpr(std::is_same<T, long>::value) { dtypid = H5T_STD_I64LE; msg = "i64"; }
	else if constexpr(std::is_same<T, float>::value) { dtypid = H5T_IEEE_F32LE; msg = "f32"; }
	else if constexpr(std::is_same<T, double>::value) { dtypid = H5T_IEEE_F64LE; msg = "f64"; }
	else { /*nop*/ }
	#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "dtypid " << msg << " " << dtypid << "\n";
	#endif

	retval = (T) 0;
	if ( nexus_open( H5F_ACC_RDONLY ) != H5I_INVALID_HID ) {
		//if ( nexus_path_exists( dsnm ) == true ) { //H5Gget_objinfo( fileid, dsnm.c_str(), 0, NULL) ) >= 0 ) {
			dsetid = H5Dopen2( fileid, dsnm.c_str(), H5P_DEFAULT );
			if ( dsetid != H5I_INVALID_HID ) {
				//no plist expected for scalars
				dspcid = H5Dget_space( dsetid );
				if ( dspcid != H5I_INVALID_HID ) {
					const int ndims = H5Sget_simple_extent_ndims( dspcid );
					if ( ndims == 0 ) {
						if ( H5Dread( dsetid, dtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &retval ) >= 0 ) {
	#ifdef MYHDF5_VERBOSE_LEVEL_ONE
							cout << "Reading " << dsnm << " success" << "\n";
	#endif
							return nexus_close( MYHDF5_SUCCESS );
						}
						cerr << "Reading " << dsnm << " failed, read !" << "\n";
						return nexus_close( MYHDF5_DSETREAD_FAILED );
					}
					cerr << "Reading " << dsnm << " failed, ndims != 0 !" << "\n";
					return nexus_close( MYHDF5_DSPACE_FAILED );
				}
				cerr << "Reading " << dsnm << " failed, dspc !" << "\n";
				return nexus_close( MYHDF5_DSPCACCESS_FAILED );
			}
			cerr << "Reading " << dsnm << " failed, dset !" << "\n";
			return nexus_close( MYHDF5_DSETACCESS_FAILED );
		//}
		//cerr << "Reading " << dsnm << " failed, does not exists !" << "\n";
		//return nexus_close( MYHDF5_DSETACCESS_FAILED );
	}
	cerr << "Opening file " << h5resultsfn << " failed, output cleared !" << "\n";
	return nexus_close( MYHDF5_FOPEN_FAILED );
}

template int HdfFiveSeqHdl::nexus_read( const string dsnm, unsigned char & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, char & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, unsigned short & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, short & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, unsigned int & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, int & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, unsigned long & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, long & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, float & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, double & retval );

int HdfFiveSeqHdl::nexus_read( const string dsnm, string & retval )
{
	//open file in read-only mode to not modify the file access timestamps
	retval = "";
	if ( nexus_open( H5F_ACC_RDONLY ) != H5I_INVALID_HID ) {
		//if ( nexus_path_exists( dsnm ) == true ) { //H5Gget_objinfo( fileid, dsnm.c_str(), 0, NULL) ) >= 0 ) {
		dsetid = H5Dopen2( fileid, dsnm.c_str(), H5P_DEFAULT );
		if ( dsetid != H5I_INVALID_HID ) {
			dtypid = H5Dget_type(dsetid);
			if ( dtypid != H5I_INVALID_HID) {
				//https://docs.hdfgroup.org/hdf5/develop/group___h5_t.html#ga364545c053f925fec65880b235e37898
				H5T_class_t which_type = H5Tget_class(dtypid);
				if ( which_type == H5T_STRING ) {
					dspcid = H5Dget_space(dsetid);
					if ( dspcid != H5I_INVALID_HID ) {
						//analyze memory layout, encoding, and formatting of the string
						int status = 0;
						//H5S_class_t extent = H5Sget_simple_extent_type(dspcid);
						hsize_t dims[1] = {0};
						hsize_t maxdims[1] = {0};
						int ndims = H5Sget_simple_extent_dims(dspcid, dims, maxdims);
						htri_t tribool = H5Tis_variable_str(dtypid);
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
						H5T_str_t strpad = H5Tget_strpad(dtypid);
#endif
						H5T_cset_t cset = H5Tget_cset(dtypid);
						size_t max_size = H5Tget_size(dtypid);
						//##MK::improve fault tolerance

#ifdef MYHDF5_VERBOSE_LEVEL_TWO
						cout << "ndims " << ndims << "\n";
						for ( int dim = 0; dim < ndims; dim++ ) {
							cout << "dim[" << dim << "] " << dims[dim] << "\n";
						}
						//maxdims
						if ( tribool > 0 ) { cout << "variable" << "\n"; }
						else if ( tribool == 0 ) { cout << "fixed" << "\n"; }
						else { cerr << "string memory error!" << "\n"; }
						switch(strpad)
						{
							case 0: { cout << "strpad nullterm" << "\n"; break; }
							case 1: { cout << "strpad nullpad" << "\n"; break; }
							case 2: { cout << "strpad spacepad" << "\n"; break; }
							default: { cerr << "strpad error!" << "\n"; break; }
						}
						switch(cset)
						{
							case 0: { cout << "cset ascii" << "\n"; break; }
							case 1: { cout << "cset_utf8" << "\n"; break; }
							default: { cerr << "cset_error" << "\n"; break; }
						}
						cout << "max_size " << max_size << "\n";
#endif
						//read
						if ( tribool > 0 ) { //MYHDF5_IS_VARIABLE
							//https://git.rwth-aachen.de/lukas.weber2/hdf5/-/blob/0d1b43c79c479b1f88eff1cc39baaf037efccfc8/test/tvlstr.c
							if ( ndims == 0 ) { //SCALAR
								mtypid = H5Tcopy(H5T_C_S1);
								status = H5Tset_size(mtypid, H5T_VARIABLE);
								if ( status >= 0 ) {
									if ( cset == 1 ) {
										status = H5Tset_cset(mtypid, H5T_CSET_UTF8);
									}
									if ( status >= 0 ) {
										char** rbuf = (char **) malloc(1 * sizeof(char *));
										if ( rbuf != NULL ) {
											if ( H5Dread(dsetid, mtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, rbuf) >= 0 ) {
												string value = string(rbuf[0]);
												retval = value;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
												cout << "__" << value << "__" << "\n";
#endif
												status = H5Treclaim(mtypid, dspcid, H5P_DEFAULT, rbuf);
												if ( status >= 0 ) {
													free(rbuf); rbuf = NULL;
													return nexus_close( MYHDF5_SUCCESS );
												}
												free(rbuf); rbuf = NULL;
												return nexus_close( MYHDF5_DSETREAD_FAILED );
											}
											cerr << "Reading " << dsnm << " buffer reclaim failed!" << "\n";
											//status = H5Treclaim(mtypid, dspcid, H5P_DEFAULT, rbuf);
											//if ( status >= 0 ) {
											//	free(rbuf); rbuf = NULL;
											//	return nexus_close( MYHDF5_DSETREAD_FAILED );
											//}
											free(rbuf); rbuf = NULL;
											return nexus_close( MYHDF5_DSETREAD_FAILED );
										}
										cerr << "Reading " << dsnm << " buffer allocation failed!" << "\n";
										return nexus_close( MYHDF5_DSETREAD_FAILED );
									}
									cerr << "Reading " << dsnm << " failed, set_cset !" << "\n";
									return nexus_close( MYHDF5_DSETREAD_FAILED );
								}
								cerr << "Reading " << dsnm << " failed, set_size !" << "\n";
								return nexus_close( MYHDF5_DSETREAD_FAILED );
							} //scalar case handled
							else { //list_single, list_long
								cerr << "Reading " << dsnm << " dataset is not a scalar use list array read function!" << "\n";
								return nexus_close( MYHDF5_FAILED );
							}
						} //variable case handled
						else if ( tribool == 0 ) { //MYHDF5_IS_FIXED
							if ( ndims == 0 ) { //SCALAR
								mtypid = H5Tcopy(H5T_C_S1);
								max_size++;
								status = H5Tset_size(mtypid, max_size);
								if ( status >= 0 ) {
									if ( cset == 1 ) {
										status = H5Tset_cset(mtypid, H5T_CSET_UTF8);
									}
									if ( status >= 0 ) {
										char* rbuf = (char *) malloc(max_size * sizeof(char));
										if (rbuf != NULL) {
											if ( H5Dread(dsetid, mtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, rbuf) >= 0 ) {
												string value = string(rbuf);
												retval = value;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
												cout << "__" << value << "__" << "\n";
#endif
												free(rbuf); rbuf = NULL;
												return nexus_close( MYHDF5_SUCCESS );
											}
											cerr << "Reading " << dsnm << " failed!" << "\n";
											free(rbuf); rbuf = NULL;
											return nexus_close( MYHDF5_DSETREAD_FAILED );
										}
										cerr << "Reading " << dsnm << " failed, buffer allocation!" << "\n";
										return nexus_close( MYHDF5_DSETREAD_FAILED );
									}
									cerr << "Reading " << dsnm << " failed, set_cset !" << "\n";
									return nexus_close( MYHDF5_DSETREAD_FAILED );
								}
								cerr << "Reading " << dsnm << " failed, set_size !" << "\n";
								return nexus_close( MYHDF5_DSETREAD_FAILED );
							} //scalar case handled
							else { //list_single, list_long
								cerr << "Reading " << dsnm << " dataset is not a scalar use list array read function!" << "\n";
								return nexus_close( MYHDF5_FAILED );
							}
						}
						else {
							cerr << "Reading " << dsnm << " unexpected case!" << "\n";
							return nexus_close( MYHDF5_FAILED );
						}
					}
					cerr << "Reading " << dsnm << " failed, dspc !" << "\n";
					return nexus_close( MYHDF5_DSPCOPEN_FAILED );
				}
				cerr << "Reading " << dsnm << " is not a string!" << "\n";
				return nexus_close( MYHDF5_FAILED );
			}
			cerr << "Reading " << dsnm << " failed, dtypid !" << "\n";
			return nexus_close( MYHDF5_DSPCREAD_FAILED );
		}
		cerr << "Reading " << dsnm << " failed, dsetid !" << "\n";
		return nexus_close( MYHDF5_DSETREAD_FAILED );
	}
	cerr << "Opening file " << h5resultsfn << " failed!" << "\n";
	return nexus_close( MYHDF5_FOPEN_FAILED );
}

/*
//deprecated, was refactored!
template<typename T> int HdfFiveSeqHdl::nexus_read( const string dsnm, vector<T> & retval )
{
	//open file in read-only mode to not modified the file access timestamps
	hid_t dtype = 0; //H5T_STD_U32LE;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
	cout << "dtype " << dtype << "\n";
#endif
	size_t max_characters = 0;
	if constexpr(std::is_same<T, unsigned char>::value) {
		dtype = H5T_STD_U8LE;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "dtype u8 " << dtype << "\n";
#endif
	}
	else if constexpr(std::is_same<T, char>::value) {
		dtype = H5T_STD_I8LE;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "dtype i8 " << dtype << "\n";
#endif
	}
	else if constexpr(std::is_same<T, unsigned short>::value) {
		dtype = H5T_STD_U16LE;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "dtype u16 " << dtype << "\n";
#endif
	}
	else if constexpr(std::is_same<T, short>::value) {
		dtype = H5T_STD_I16LE;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "dtype i16 " << dtype << "\n";
#endif
	}
	else if constexpr(std::is_same<T, unsigned int>::value) {
		dtype = H5T_STD_U32LE;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "dtype u32 " << dtype << "\n";
#endif
	}
	else if constexpr(std::is_same<T, int>::value) {
		dtype = H5T_STD_I32LE;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "dtype i32 " << dtype << "\n";
#endif
	}
	else if constexpr(std::is_same<T, unsigned long>::value) {
		dtype = H5T_STD_U64LE;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "dtype u64 " << dtype << "\n";
#endif
	}
	else if constexpr(std::is_same<T, long>::value) {
		dtype = H5T_STD_I64LE;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "dtype i64 " << dtype << "\n";
#endif
	}
	else if constexpr(std::is_same<T, float>::value) {
		dtype = H5T_IEEE_F32LE;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "dtype f32 " << dtype << "\n";
#endif
	}
	else if constexpr(std::is_same<T, double>::value) {
		dtype = H5T_IEEE_F64LE;
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "dtype f64 " << dtype << "\n";
#endif
	}
	else {
		//nop
	}

	retval = vector<T>();
	if ( nexus_open( H5F_ACC_RDONLY ) != H5I_INVALID_HID ) {
		//if ( nexus_path_exists( dsnm ) == true ) { //H5Gget_objinfo( fileid, dsnm.c_str(), 0, NULL) ) >= 0 ) {
			dsetid = H5Dopen2( fileid, dsnm.c_str(), H5P_DEFAULT );
			if ( dsetid != H5I_INVALID_HID ) {

				if constexpr(std::is_same<T, string>::value) {
					max_characters = H5Dget_storage_size( dsetid );
				}

				plistid = H5Dget_create_plist( dsetid );
				if ( plistid != H5I_INVALID_HID ) {
					unsigned int flags;
					unsigned filter_info;
					H5Z_filter_t filter_type;
					int n_filters = H5Pget_nfilters(plistid);
					for ( int i = 0; i < n_filters; i++) {
						size_t n_elements = 0;
						filter_type = H5Pget_filter2( plistid, (unsigned) i, &flags, &n_elements, NULL, 0, NULL, &filter_info);
						if ( filter_type >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
							cout << "Detected the following filter: ";
							switch(filter_type)
							{
								case H5Z_FILTER_DEFLATE: { cout << "H5Z_FILTER_DEFLATE" << "\n"; break; }
								case H5Z_FILTER_SZIP: { cout << "H5Z_FILTER_SZIP" << "\n"; break; }
								default: { cout << "OTHER" << "\n"; break; }
							}
#endif
						}
						else {
							cerr << "Reading " << dsnm << " failed, reading filter !" << "\n";
							return nexus_close( MYHDF5_PLSTFILTER_FAILED ); //, true, true, true, true );
						}
					} //all filters inspected

					dspcid = H5Dget_space( dsetid );
					if ( dspcid != H5I_INVALID_HID ) {
						size_t n_values = 1;
						const int ndims = H5Sget_simple_extent_ndims( dspcid );
						if ( ndims == 0 ) {
							//nop, scalar
						}
						else if ( ndims == 1 ) {
							hsize_t dims[1] = { 0 };
							if ( H5Sget_simple_extent_dims(dspcid, dims, NULL) == ndims ) {
								n_values = dims[0];
							}
							else {
								return nexus_close( MYHDF5_DSETACCESS_FAILED );
							}
						}
						else if ( ndims == 2 ) {
							hsize_t dims[2] = { 0, 0 };
							if ( H5Sget_simple_extent_dims(dspcid, dims, NULL) == ndims ) {
								n_values = (size_t) dims[0] * (size_t) dims[1];
							}
							else {
								return nexus_close( MYHDF5_DSETACCESS_FAILED );
							}
						}
						else {
							cerr << "Reading " << dsnm << " failed, ndims is not 0, 1, or 2 but ndims is " << ndims << "\n";
							return nexus_close( MYHDF5_NOT_EXPECTED);
						}

						if constexpr(std::is_same<T, string>::value) {
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
							cout << "ndims " << ndims << " n_values " << n_values << "\n";
#endif
							vector<char> tmp = vector<char>( max_characters, 0x00 );
							//https://github.com/HDFGroup/hdf5/issues/544
							hid_t memtype = H5Tcopy(H5T_C_S1);
							if ( memtype >= 0 ) {
								if ( H5Tset_size( memtype, H5T_VARIABLE ) >= 0 ) {
									//if ( H5Tset_strpad(memtype, H5T_STR_NULLTERM) >= 0 ) {
									//	if ( H5Tset_cset(memtype, H5T_CSET_UTF8) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
									cout << "memtype str " << memtype << "\n";
#endif
									if ( H5Dread( dsetid, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp.data() ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
										cout << "Reading " << dsnm << " success" << "\n";
#endif
										retval.push_back( string( tmp.begin(), tmp.end() ) );
										if ( H5Tclose( memtype ) >= 0 ) {
											return nexus_close( MYHDF5_SUCCESS );
										}
										cerr << "Reading " << dsnm << " failed, memtype close !" << "\n";
										return MYHDF5_NOT_EXPECTED;
									}
									cerr << "Reading " << dsnm << " failed !" << "\n";
									if ( H5Tclose( memtype ) >= 0 )
										return nexus_close( MYHDF5_NOT_EXPECTED );
									return nexus_close( MYHDF5_DSETREAD_FAILED );
								}
								cerr << "memtype str set size " << memtype << "\n";
								if ( H5Tclose( memtype ) >= 0 )
									return nexus_close( MYHDF5_NOT_EXPECTED );
								return nexus_close( MYHDF5_NOT_EXPECTED );

										//cerr << "Reading " << dsnm << " failed !" << "\n";
										//if ( H5Tclose( memtype ) >= 0 )
										//	return nexus_close( MYHDF5_NOT_EXPECTED );
										//return nexus_close( MYHDF5_NOT_EXPECTED );
									//}
									//cerr << "memtype str strp set " << memtype << "\n";
									//if ( H5Tclose( memtype ) >= 0 )
									//	return nexus_close( MYHDF5_NOT_EXPECTED );
									//return nexus_close( MYHDF5_NOT_EXPECTED );
							}
							cerr << "Reading " << dsnm << " failed, set size !" << "\n";
							if ( H5Tclose( memtype ) >= 0 )
								return nexus_close( MYHDF5_NOT_EXPECTED );
							return nexus_close( MYHDF5_NOT_EXPECTED );
						}
						else {
							retval = vector<T>( n_values, (T) 0 );
							if ( H5Dread( dsetid, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, retval.data() ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
								cout << "Reading " << dsnm << " success" << "\n";
#endif
								return nexus_close( MYHDF5_SUCCESS );
							}
						}

						cerr << "Reading " << dsnm << " failed, read !" << "\n";
						return nexus_close( MYHDF5_DSETREAD_FAILED );
					}
					cerr << "Reading " << dsnm << " failed, dspc !" << "\n";
					return nexus_close( MYHDF5_DSPCACCESS_FAILED );
				}
				cerr << "Reading " << dsnm << " failed, read plist !" << "\n";
				return nexus_close( MYHDF5_PLSTACCESS_FAILED );
			}
			cerr << "Reading " << dsnm << " failed, dset !" << "\n";
			return nexus_close( MYHDF5_DSETACCESS_FAILED );
		//}
		//cerr << "Reading " << dsnm << " failed, does not exists !" << "\n";
		//return nexus_close( MYHDF5_DSETACCESS_FAILED );
	}
	cerr << "Opening file " << h5resultsfn << " failed, output cleared !" << "\n";
	return MYHDF5_FOPEN_FAILED;
}
*/

template<typename T> int HdfFiveSeqHdl::nexus_read( const string dsnm, vector<T> & retval )
{
	//open file in read-only mode to not modified the file last access time stamp
	//reads list_single and list_long of base C/C++ data types, for vector<string> there is an own function
	dtypid = H5I_INVALID_HID;
	string msg = "";
	if constexpr(std::is_same<T, unsigned char>::value) { dtypid = H5T_STD_U8LE; msg = "u8"; }
	else if constexpr(std::is_same<T, char>::value) { dtypid = H5T_STD_I8LE; msg = "i8"; }
	else if constexpr(std::is_same<T, unsigned short>::value) { dtypid = H5T_STD_U16LE; msg = "u16"; }
	else if constexpr(std::is_same<T, short>::value) { dtypid = H5T_STD_I16LE; msg = "i16"; }
	else if constexpr(std::is_same<T, unsigned int>::value) { dtypid = H5T_STD_U32LE; msg = "u32"; }
	else if constexpr(std::is_same<T, int>::value) { dtypid = H5T_STD_I32LE; msg = "i32"; }
	else if constexpr(std::is_same<T, unsigned long>::value) { dtypid = H5T_STD_U64LE; msg = "u64"; }
	else if constexpr(std::is_same<T, long>::value) { dtypid = H5T_STD_I64LE; msg = "i64"; }
	else if constexpr(std::is_same<T, float>::value) { dtypid = H5T_IEEE_F32LE; msg = "f32"; }
	else if constexpr(std::is_same<T, double>::value) { dtypid = H5T_IEEE_F64LE; msg = "f64"; }
	else { /*nop*/ }
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "dtypid " << msg << " " << dtypid << "\n";
#endif

	retval = vector<T>();
	if ( nexus_open( H5F_ACC_RDONLY ) != H5I_INVALID_HID ) {
		//if ( nexus_path_exists( dsnm ) == true ) { //H5Gget_objinfo( fileid, dsnm.c_str(), 0, NULL) ) >= 0 ) {
			dsetid = H5Dopen2( fileid, dsnm.c_str(), H5P_DEFAULT );
			if ( dsetid != H5I_INVALID_HID ) {
				plistid = H5Dget_create_plist( dsetid );
				if ( plistid != H5I_INVALID_HID ) {
					unsigned int flags;
					unsigned filter_info;
					H5Z_filter_t filter_type;
					int n_filters = H5Pget_nfilters(plistid);
					for ( int i = 0; i < n_filters; i++) {
						size_t n_elements = 0;
						filter_type = H5Pget_filter2( plistid, (unsigned) i, &flags, &n_elements, NULL, 0, NULL, &filter_info);
						if ( filter_type >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
							cout << "Detected the following filter: ";
							switch(filter_type)
							{
								case H5Z_FILTER_DEFLATE: { cout << "H5Z_FILTER_DEFLATE" << "\n"; break; }
								case H5Z_FILTER_SZIP: { cout << "H5Z_FILTER_SZIP" << "\n"; break; }
								default: { cout << "OTHER" << "\n"; break; }
							}
#endif
						}
						else {
							cerr << "Reading " << dsnm << " failed, reading filter !" << "\n";
							return nexus_close( MYHDF5_PLSTFILTER_FAILED );
						}
					} //all filters inspected

					dspcid = H5Dget_space( dsetid );
					if ( dspcid != H5I_INVALID_HID ) {
						size_t n_values = 1;
						const int ndims = H5Sget_simple_extent_ndims( dspcid );
						if ( ndims == 1 ) {
							hsize_t dims[1] = { 0 };
							if ( H5Sget_simple_extent_dims(dspcid, dims, NULL) == ndims ) {
								n_values = dims[0];
							}
							else {
								return nexus_close( MYHDF5_DSETACCESS_FAILED );
							}
						}
						else if ( ndims == 2 ) {
							hsize_t dims[2] = { 0, 0 };
							if ( H5Sget_simple_extent_dims(dspcid, dims, NULL) == ndims ) {
								n_values = (size_t) dims[0] * (size_t) dims[1];
							}
							else {
								return nexus_close( MYHDF5_DSETACCESS_FAILED );
							}
						}
						/*
						else if ( ndims == 3 ) {
							hsize_t dims[3] = { 0, 0, 0 };
							if ( H5Sget_simple_extent_dims(dspcid, dims, NULL) == ndims ) {
								n_values = (size_t) dims[0] * (size_t) dims[1] * (size_t) dims[2];
							}
							else {
								return nexus_close( MYHDF5_DSETACCESS_FAILED );
							}
						}
						*/
						else {
							cerr << "Reading " << dsnm << " failed, ndims is " << ndims << " but only 1, and 2 is supported!" << "\n";
							return nexus_close( MYHDF5_NOT_EXPECTED);
						}

						retval = vector<T>( n_values, (T) 0 );
						if ( H5Dread( dsetid, dtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, retval.data() ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
							cout << "Reading " << dsnm << " success" << "\n";
#endif
							return nexus_close( MYHDF5_SUCCESS );
						}
						cerr << "Reading " << dsnm << " failed, read !" << "\n";
						return nexus_close( MYHDF5_DSETREAD_FAILED );
					}
					cerr << "Reading " << dsnm << " failed, dspc !" << "\n";
					return nexus_close( MYHDF5_DSPCACCESS_FAILED );
				}
				cerr << "Reading " << dsnm << " failed, read plist !" << "\n";
				return nexus_close( MYHDF5_PLSTACCESS_FAILED );
			}
			cerr << "Reading " << dsnm << " failed, dset !" << "\n";
			return nexus_close( MYHDF5_DSETACCESS_FAILED );
		//}
		//cerr << "Reading " << dsnm << " failed, does not exists !" << "\n";
		//return nexus_close( MYHDF5_DSETACCESS_FAILED );
	}
	cerr << "Opening file " << h5resultsfn << " failed, output cleared !" << "\n";
	return nexus_close( MYHDF5_FOPEN_FAILED );
}

template int HdfFiveSeqHdl::nexus_read( const string dsnm, vector<unsigned char> & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, vector<char> & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, vector<unsigned short> & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, vector<short> & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, vector<unsigned int> & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, vector<int> & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, vector<unsigned long> & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, vector<long> & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, vector<float> & retval );
template int HdfFiveSeqHdl::nexus_read( const string dsnm, vector<double> & retval );

int HdfFiveSeqHdl::nexus_read( const string dsnm, vector<string> & retval )
{
	//open file in read-only mode to not modify the file access timestamps
	retval = vector<string>();
	if ( nexus_open( H5F_ACC_RDONLY ) != H5I_INVALID_HID ) {
		//if ( nexus_path_exists( dsnm ) == true ) { //H5Gget_objinfo( fileid, dsnm.c_str(), 0, NULL) ) >= 0 ) {
		dsetid = H5Dopen2( fileid, dsnm.c_str(), H5P_DEFAULT );
		if ( dsetid != H5I_INVALID_HID ) {
			dtypid = H5Dget_type(dsetid);
			if ( dtypid != H5I_INVALID_HID) {
				//https://docs.hdfgroup.org/hdf5/develop/group___h5_t.html#ga364545c053f925fec65880b235e37898
				H5T_class_t which_type = H5Tget_class(dtypid);
				if ( which_type == H5T_STRING ) {
					dspcid = H5Dget_space(dsetid);
					if ( dspcid != H5I_INVALID_HID ) {
						//analyze memory layout, encoding, and formatting of the string
						int status = 0;
						//H5S_class_t extent = H5Sget_simple_extent_type(dspcid);
						hsize_t dims[1] = {0};
						hsize_t maxdims[1] = {0};
						int ndims = H5Sget_simple_extent_dims(dspcid, dims, maxdims);
						htri_t tribool = H5Tis_variable_str(dtypid);
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
						H5T_str_t strpad = H5Tget_strpad(dtypid);
#endif
						H5T_cset_t cset = H5Tget_cset(dtypid);
						size_t max_size = H5Tget_size(dtypid);
						//##MK::improve fault tolerance

#ifdef MYHDF5_VERBOSE_LEVEL_TWO
						cout << "ndims " << ndims << "\n";
						for ( int dim = 0; dim < ndims; dim++ ) {
							cout << "dim[" << dim << "] " << dims[dim] << "\n";
						}
						//maxdims
						if ( tribool > 0 ) { cout << "variable" << "\n"; }
						else if ( tribool == 0 ) { cout << "fixed" << "\n"; }
						else { cerr << "string memory error!" << "\n"; }
						switch(strpad)
						{
							case 0: { cout << "strpad nullterm" << "\n"; break; }
							case 1: { cout << "strpad nullpad" << "\n"; break; }
							case 2: { cout << "strpad spacepad" << "\n"; break; }
							default: { cerr << "strpad error!" << "\n"; break; }
						}
						switch(cset)
						{
							case 0: { cout << "cset ascii" << "\n"; break; }
							case 1: { cout << "cset_utf8" << "\n"; break; }
							default: { cerr << "cset_error" << "\n"; break; }
						}
						cout << "max_size " << max_size << "\n";
#endif
						//read
						if ( tribool > 0 ) { //MYHDF5_IS_VARIABLE
							//https://git.rwth-aachen.de/lukas.weber2/hdf5/-/blob/0d1b43c79c479b1f88eff1cc39baaf037efccfc8/test/tvlstr.c
							if ( ndims == 0 ) { //SCALAR
								mtypid = H5Tcopy(H5T_C_S1);
								status = H5Tset_size(mtypid, H5T_VARIABLE);
								if ( status >= 0 ) {
									if ( cset == 1 ) {
										status = H5Tset_cset(mtypid, H5T_CSET_UTF8);
									}
									if ( status >= 0 ) {
										char** rbuf = (char **) malloc(1 * sizeof(char *));
										if ( rbuf != NULL ) {
											if ( H5Dread(dsetid, mtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, rbuf) >= 0 ) {
												string value = string(rbuf[0]);
												retval.push_back(value);
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
												cout << "__" << value << "__" << "\n";
#endif
												status = H5Treclaim(mtypid, dspcid, H5P_DEFAULT, rbuf);
												if ( status >= 0 ) {
													free(rbuf); rbuf = NULL;
													return nexus_close( MYHDF5_SUCCESS );
												}
												free(rbuf); rbuf = NULL;
												return nexus_close( MYHDF5_DSETREAD_FAILED );
											}
											cerr << "Reading " << dsnm << " buffer reclaim failed!" << "\n";
											//status = H5Treclaim(mtypid, dspcid, H5P_DEFAULT, rbuf);
											//if ( status >= 0 ) {
											//	free(rbuf); rbuf = NULL;
											//	return nexus_close( MYHDF5_DSETREAD_FAILED );
											//}
											free(rbuf); rbuf = NULL;
											return nexus_close( MYHDF5_DSETREAD_FAILED );
										}
										cerr << "Reading " << dsnm << " buffer allocation failed!" << "\n";
										return nexus_close( MYHDF5_DSETREAD_FAILED );
									}
									cerr << "Reading " << dsnm << " failed, set_cset !" << "\n";
									return nexus_close( MYHDF5_DSETREAD_FAILED );
								}
								cerr << "Reading " << dsnm << " failed, set_size !" << "\n";
								return nexus_close( MYHDF5_DSETREAD_FAILED );
							} //scalar case handled
							else { //LIST_SINGLE or LIST_LONG
								mtypid = H5Tcopy(H5T_C_S1);
								status = H5Tset_size(mtypid, H5T_VARIABLE);
								if ( status >= 0 ) {
									if ( cset == 1 ) {
										status = H5Tset_cset(mtypid, H5T_CSET_UTF8);
									}
									if ( status >= 0 ) {
										char** rbuf = (char **) malloc(dims[0] * sizeof(char *));
										if (rbuf != NULL) {
											if ( H5Dread(dsetid, mtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, rbuf) >= 0 ) {
												for ( size_t i = 0; i < dims[0]; i++ ) {
													string value = string(rbuf[i]);
													retval.push_back(value);
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
													cout << "__" << value << "__" << "\n";
#endif
												}
												status = H5Treclaim(mtypid, dspcid, H5P_DEFAULT, rbuf);
												if ( status >= 0 ) {
													free(rbuf); rbuf = NULL;
													return nexus_close( MYHDF5_SUCCESS );
												}
												cerr << "Reading " << dsnm << " buffer reclaim failed!" << "\n";
												free(rbuf); rbuf = NULL;
												return nexus_close( MYHDF5_DSETREAD_FAILED );
												//although retval might have the right values
												//we need to be picky to avoid a resource leakage
											}
											cerr << "Reading " << dsnm << " buffer reclaim failed!" << "\n";
											//status = H5Treclaim(mtypid, dspcid, H5P_DEFAULT, rbuf);
											//if ( status >= 0 ) {
											//	free(rbuf); rbuf = NULL;
											//	return nexus_close( MYHDF5_DSETREAD_FAILED );
											//}
											free(rbuf); rbuf = NULL;
											return nexus_close( MYHDF5_DSETREAD_FAILED );
										}
										cerr << "Reading " << dsnm << " failed, buffer allocation !" << "\n";
										return nexus_close( MYHDF5_DSETREAD_FAILED );
									}
									cerr << "Reading " << dsnm << " failed, set_cset !" << "\n";
									return nexus_close( MYHDF5_DSETREAD_FAILED );
								}
								cerr << "Reading " << dsnm << " failed, set_size !" << "\n";
								return nexus_close( MYHDF5_DSETREAD_FAILED );
							} //list case handled
						} //variable case handled
						else if ( tribool == 0 ) { //MYHDF5_IS_FIXED
							if ( ndims == 0 ) { //SCALAR
								mtypid = H5Tcopy(H5T_C_S1);
								max_size++;
								status = H5Tset_size(mtypid, max_size);
								if ( status >= 0 ) {
									if ( cset == 1 ) {
										status = H5Tset_cset(mtypid, H5T_CSET_UTF8);
									}
									if ( status >= 0 ) {
										char* rbuf = (char *) malloc(max_size * sizeof(char));
										if (rbuf != NULL) {
											if ( H5Dread(dsetid, mtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, rbuf) >= 0 ) {
												string value = string(rbuf);
												retval.push_back(value);
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
												cout << "__" << value << "__" << "\n";
#endif
												free(rbuf); rbuf = NULL;
												return nexus_close( MYHDF5_SUCCESS );
											}
											cerr << "Reading " << dsnm << " failed!" << "\n";
											free(rbuf); rbuf = NULL;
											return nexus_close( MYHDF5_DSETREAD_FAILED );
										}
										cerr << "Reading " << dsnm << " failed, buffer allocation!" << "\n";
										return nexus_close( MYHDF5_DSETREAD_FAILED );
									}
									cerr << "Reading " << dsnm << " failed, set_cset !" << "\n";
									return nexus_close( MYHDF5_DSETREAD_FAILED );
								}
								cerr << "Reading " << dsnm << " failed, set_size !" << "\n";
								return nexus_close( MYHDF5_DSETREAD_FAILED );
							} //scalar case handled
							else { //LIST_SINGLE or LIST_LONG
								mtypid = H5Tcopy(H5T_C_S1);
								max_size++;
								status = H5Tset_size(mtypid, max_size);
								if ( status >= 0 ) {
									if ( cset == 1 ) {
										status = H5Tset_cset(mtypid, H5T_CSET_UTF8);
									}
									if ( status >= 0 ) {
										char** rbuf = (char **) malloc(dims[0] * sizeof(char *));
										if (rbuf != NULL) {
											rbuf[0] = (char *) malloc(dims[0] * max_size * sizeof(char));
											for( size_t i = 1; i < dims[0]; i++ ) {
												rbuf[i] = rbuf[0] + i * max_size;
											}
											if ( rbuf[0] != NULL ) {
												if ( H5Dread(dsetid, mtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, rbuf[0]) >= 0 ) {
													for ( size_t i = 0; i < dims[0]; i++ ) {
														string value = string(rbuf[i]);
														retval.push_back(value);
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
														cout << "__" << value << "__" << "\n";
#endif
													}
													free(rbuf[0]); rbuf[0] = NULL;
													free(rbuf); rbuf = NULL;
													return nexus_close( MYHDF5_SUCCESS );
												}
												cerr << "Reading " << dsnm << " failed!" << "\n";
												free(rbuf[0]); rbuf[0] = NULL;
												free(rbuf); rbuf = NULL;
												return nexus_close( MYHDF5_DSETREAD_FAILED );
											}
											cerr << "Reading " << dsnm << " failed!" << "\n";
											free(rbuf); rbuf = NULL;
											return nexus_close( MYHDF5_DSETREAD_FAILED );
										}
										cerr << "Reading " << dsnm << " failed, buffer allocation !" << "\n";
										return nexus_close( MYHDF5_DSETREAD_FAILED );
									}
									cerr << "Reading " << dsnm << " failed, set_cset !" << "\n";
									return nexus_close( MYHDF5_DSETREAD_FAILED );
								}
								cerr << "Reading " << dsnm << " failed, set_size !" << "\n";
								return nexus_close( MYHDF5_DSETREAD_FAILED );
							} //list case handled
						}
						else {
							cerr << "Reading " << dsnm << " unexpected case!" << "\n";
							return nexus_close( MYHDF5_FAILED );
						}
					}
					cerr << "Reading " << dsnm << " failed, dspc !" << "\n";
					return nexus_close( MYHDF5_DSPCOPEN_FAILED );
				}
				cerr << "Reading " << dsnm << " is not a string!" << "\n";
				return nexus_close( MYHDF5_FAILED );
			}
			cerr << "Reading " << dsnm << " failed, dtypid !" << "\n";
			return nexus_close( MYHDF5_DSPCREAD_FAILED );
		}
		cerr << "Reading " << dsnm << " failed, dsetid !" << "\n";
		return nexus_close( MYHDF5_DSETREAD_FAILED );
	}
	cerr << "Opening file " << h5resultsfn << " failed!" << "\n";
	return nexus_close( MYHDF5_FOPEN_FAILED );
}


//low-level, scalar
template<typename T> int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, const T val )
{
	//MK::this function may only be called from nexus_write_attribute as it adds a typed and named attribute
	//to an (assumed existent) dataset referred to by dsetid inside an RDWR-opened HDF5 file referred to by fileid
	//for a specific type with the name in the keyword and the value in the val variable
	//this function will release attribute related resources
	//scalar basic type attribute
	adtypid = H5I_INVALID_HID;
	string msg = "";
	if constexpr(std::is_same<T, unsigned char>::value) { adtypid = H5T_STD_U8LE; msg = "u8"; }
	else if constexpr(std::is_same<T, char>::value) { adtypid = H5T_STD_I8LE; msg = "i8"; }
	else if constexpr(std::is_same<T, unsigned short>::value) { adtypid = H5T_STD_U16LE; msg = "u16"; }
	else if constexpr(std::is_same<T, short>::value) { adtypid = H5T_STD_I16LE; msg = "i16"; }
	else if constexpr(std::is_same<T, unsigned int>::value) { adtypid = H5T_STD_U32LE; msg = "u32"; }
	else if constexpr(std::is_same<T, int>::value) { adtypid = H5T_STD_I32LE; msg = "i32"; }
	else if constexpr(std::is_same<T, unsigned long>::value) { adtypid = H5T_STD_U64LE; msg = "u64"; }
	else if constexpr(std::is_same<T, long>::value) { adtypid = H5T_STD_I64LE; msg = "i64"; }
	else if constexpr(std::is_same<T, float>::value) { adtypid = H5T_IEEE_F32LE; msg = "f32"; }
	else if constexpr(std::is_same<T, double>::value) { adtypid = H5T_IEEE_F64LE; msg = "f64"; }
	else { /*nop*/ }
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "adtypid " << msg << " " << adtypid << "\n";
#endif

	if ( adtypid != H5I_INVALID_HID ) {
		aspcid = H5Screate(H5S_SCALAR);
		if ( aspcid != H5I_INVALID_HID ) {
			attrid = H5Acreate(dsetid, keyword.c_str(), adtypid, aspcid, H5P_DEFAULT, H5P_DEFAULT);
			if ( attrid >= 0 ) {
				vector<T> adata;
				adata.push_back(val);
				if ( H5Awrite(attrid, adtypid, adata.data() ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
					cout << "Writing attribute " << keyword << " success" << "\n";
#endif
					return nexus_attribute_close( MYHDF5_SUCCESS );
				}
				cerr << "Writing attribute " << keyword << " failed!" << "\n";
				return nexus_attribute_close( MYHDF5_ATTRWRITE_FAILED );
			}
			cerr << "Creating attribute " << keyword << " failed, attrid !" << "\n";
			return nexus_attribute_close( MYHDF5_ATTRCREATE_FAILED );
		}
		cerr << "Creating attribute " << keyword << " failed, aspc !" << "\n";
		return nexus_attribute_close( MYHDF5_ATTRCREATE_FAILED );
	}
	cerr << "Formatting cset attribute " << keyword << " failed, adtypid !" << "\n";
	return nexus_attribute_close( MYHDF5_ATTRCREATE_FAILED );
}

template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, const unsigned char val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, const char val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, const unsigned short val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, const short val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, const unsigned int val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, const int val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, const unsigned long val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, const long val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, const float val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, const double val );

int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, const string val )
{
	//this function will release attribute related resources
	//scalar string attribute
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
	unsigned char dimensionality = MYHDF5_SCALAR;
#endif
	//using the h5py3.9 defaults
	unsigned char memory = MYHDF5_IS_VARIABLE;
	unsigned char terminator = MYHDF5_NULLTERM;
	unsigned char encoding = MYHDF5_UTF8;
	size_t max_characters = val.length();
	max_characters++;
	//we use Cstyle arrays which is why we need space for the terminator

#ifdef MYHDF5_VERBOSE_LEVEL_TWO
	cout << "val.length() " << val.length() << "\n";
	cout << "val.size() " << val.size() << "\n";
	cout << "val.c_str() __" << val.c_str() << "__" << "\n";
	cout << "dimensionality " << (int) dimensionality << "\n";
	cout << "memory " << (int) memory << "\n";
	cout << "terminator " << (int) terminator << "\n";
	cout << "encoding " << (int) encoding << "\n";
	cout << "ndims " << 0 << "\n";
	cout << "sdims " << max_characters << "\n";
#endif
	//MK::this function may only be called from nexus_write_attribute as it adds a typed and named attribute
	//to an (assumed existent) dataset referred to by dsetid inside an RDWR-opened HDF5 file referred to by fileid
	//for a specific type with the name in the keyword and the value in the val variable
	//see also an example for string arrays here https://svn.ssec.wisc.edu/repos/geoffc/C/HDF5/Examples_by_API/h5ex_t_stringatt.c

	int status = 0;
	adtypid = H5Tcopy(H5T_C_S1);
	if ( adtypid != H5I_INVALID_HID ) {
		if ( memory == MYHDF5_IS_VARIABLE ) {
			status = H5Tset_size(adtypid, H5T_VARIABLE);
		}
		else { //MYHDF5_IS_FIXED
			status = H5Tset_size(adtypid, max_characters);
		}
		if ( status >= 0 ) {
			switch(terminator)
			{
				case MYHDF5_NULLTERM: { status = H5Tset_strpad(adtypid, H5T_STR_NULLTERM); break; }
				case MYHDF5_NULLPAD: { status = H5Tset_strpad(adtypid, H5T_STR_NULLPAD); break; }
				case MYHDF5_SPACEPAD: { status = H5Tset_strpad(adtypid, H5T_STR_SPACEPAD); break; }
				default: { status = -1; break; }
			}
			if ( status >= 0 ) {
				switch(encoding)
				{
					case MYHDF5_ASCII: { status = H5Tset_cset(adtypid, H5T_CSET_ASCII); break; }
					case MYHDF5_UTF8: { status = H5Tset_cset(adtypid, H5T_CSET_UTF8); break; }
					default: { status = -1; break; }
				}
				if ( status >= 0 ) {
					aspcid = H5Screate(H5S_SCALAR);
					if ( aspcid != H5I_INVALID_HID ) {
						attrid = H5Acreate2(dsetid, keyword.c_str(), adtypid, aspcid, H5P_DEFAULT, H5P_DEFAULT);
						if ( attrid >= 0 ) {
							if ( memory == MYHDF5_IS_VARIABLE ) {
								vector<string> adata;
								adata.push_back(val);
								//char* adata = new char[max_characters];
								//if ( adata != NULL ) {
								//	strcpy(adata, val.c_str());
								if ( H5Awrite(attrid, adtypid, adata.data() ) >= 0 ) { //adata ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
									cout << "Writing attribute " << keyword << " success" << "\n";
#endif
									//delete[] adata; adata = NULL;
									return nexus_attribute_close( MYHDF5_SUCCESS );
								}
								//cerr << "Writing attribute " << keyword << " failed!" << "\n";
								//delete[] adata; adata = NULL;
								//return nexus_attribute_close( MYHDF5_ATTRWRITE_FAILED );
								//}
								cerr << "Writing attribute " << keyword << " failed!" << "\n";
							    //delete[] adata; adata = NULL;
								return nexus_attribute_close( MYHDF5_ATTRWRITE_FAILED );
							}
							else { //MYHDF5_IS_FIXED
								char* adata = new char[max_characters];
								if ( adata != NULL ) {
									strcpy(adata, val.c_str());
									if ( H5Awrite(attrid, adtypid, adata ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
										cout << "Writing attribute " << keyword << " success" << "\n";
#endif
										delete[] adata; adata = NULL;
										return nexus_attribute_close( MYHDF5_SUCCESS );
									}
									cerr << "Writing attribute " << keyword << " failed!" << "\n";
									delete[] adata; adata = NULL;
									return nexus_attribute_close( MYHDF5_ATTRWRITE_FAILED );
								}
								cerr << "Writing attribute " << keyword << " failed!" << "\n";
								delete[] adata; adata = NULL;
								return nexus_attribute_close( MYHDF5_ATTRWRITE_FAILED );
							}
						}
						cerr << "Creating attribute " << keyword << " failed!" << "\n";
						return nexus_attribute_close( MYHDF5_ATTRCREATE_FAILED );
					}
					cerr << "Dspacing attribute " << keyword << " failed!" << "\n";
					return nexus_attribute_close( MYHDF5_ATTRCREATE_FAILED );
				}
				cerr << "Formatting cset attribute " << keyword << " failed!" << "\n";
				return nexus_attribute_close( MYHDF5_FAILED );
			}
			cerr << "Formatting str_strpad attribute " << keyword << " failed!" << "\n";
			return nexus_attribute_close( MYHDF5_FAILED );
		}
		cerr << "Formatting size attribute " << keyword << " failed!" << "\n";
		return nexus_attribute_close( MYHDF5_FAILED );
	}
	cerr << "Creating datatype attribute " << keyword << " failed!" << "\n";
	return nexus_attribute_close( MYHDF5_FAILED );
}

//low-level, list
template<typename T> int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, vector<T> & val )
{
	//MK::this function may only be called from nexus_write_attribute as it adds a typed and named attribute
	//to an (assumed existent) dataset referred to by dsetid inside an RDWR-opened HDF5 file referred to by fileid
	//for a specific type with the name in the keyword and the value in the val variable
	if ( val.size() >= 1 ) {
		adtypid = H5I_INVALID_HID;
		string msg = "";
		if constexpr(std::is_same<T, unsigned char>::value) { adtypid = H5T_STD_U8LE; msg = "u8"; }
		else if constexpr(std::is_same<T, char>::value) { adtypid = H5T_STD_I8LE; msg = "i8"; }
		else if constexpr(std::is_same<T, unsigned short>::value) { adtypid = H5T_STD_U16LE; msg = "u16"; }
		else if constexpr(std::is_same<T, short>::value) { adtypid = H5T_STD_I16LE; msg = "i16"; }
		else if constexpr(std::is_same<T, unsigned int>::value) { adtypid = H5T_STD_U32LE; msg = "u32"; }
		else if constexpr(std::is_same<T, int>::value) { adtypid = H5T_STD_I32LE; msg = "i32"; }
		else if constexpr(std::is_same<T, unsigned long>::value) { adtypid = H5T_STD_U64LE; msg = "u64"; }
		else if constexpr(std::is_same<T, long>::value) { adtypid = H5T_STD_I64LE; msg = "i64"; }
		else if constexpr(std::is_same<T, float>::value) { adtypid = H5T_IEEE_F32LE; msg = "f32"; }
		else if constexpr(std::is_same<T, double>::value) { adtypid = H5T_IEEE_F64LE; msg = "f64"; }
		else { /*nop*/ }
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
			cout << "adtypid " << msg << " " << adtypid << "\n";
#endif

		hsize_t dims[1] = { val.size() };
		aspcid = H5Screate_simple( 1, dims, NULL );
		if ( aspcid != H5I_INVALID_HID ) {
			attrid = H5Acreate2(dsetid, keyword.c_str(), adtypid, aspcid, H5P_DEFAULT, H5P_DEFAULT);
			if ( attrid >= 0 ) {
				if ( H5Awrite(attrid, adtypid, val.data() ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
					cout << "Writing attribute " << keyword << " success" << "\n";
#endif
					return nexus_attribute_close( MYHDF5_SUCCESS );
				}
				cerr << "Writing attribute " << keyword << " failed, write !" << "\n";
				return nexus_attribute_close( MYHDF5_ATTRWCOLL_FAILED );
			}
			cerr << "Writing attribute " << keyword << " failed, attrid !" << "\n";
			return nexus_attribute_close( MYHDF5_ATTRWCOLL_FAILED );
		}
		cerr << "Writing attribute " << keyword << " failed, aspc !" << "\n";
		return nexus_attribute_close( MYHDF5_ATTRCREATE_FAILED );
	}
	cerr << "Writing attribute " << keyword << " failed, attribute val.size() == 0 !" << "\n";
	return nexus_attribute_close( MYHDF5_FAILED );
}

template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, vector<unsigned char> & val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, vector<char> & val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, vector<unsigned short> & val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, vector<short> & val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, vector<unsigned int> & val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, vector<int> & val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, vector<unsigned long> & val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, vector<long> & val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, vector<float> & val );
template int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, vector<double> & val );

int HdfFiveSeqHdl::nexus_write_attribute_value( const string keyword, vector<string> & val )
{
	if ( val.size() > 0 ) {
		//this function will release attribute related resources

		//always a one-d array with at least one string
		unsigned char dimensionality = MYHDF5_LIST_SINGLE;
		if ( val.size() > 1 ) {
			dimensionality = MYHDF5_LIST_LONG;
		}
		//using the h5py3.9 defaults
		unsigned char memory = MYHDF5_IS_VARIABLE;
		unsigned char terminator = MYHDF5_NULLTERM;
		unsigned char encoding = MYHDF5_UTF8;
		size_t max_characters = 0;
		for( auto it = val.begin(); it != val.end(); it++ ) {
			if ( it->length() < max_characters )
				continue;
			else
				max_characters = it->length();
		}
		max_characters++;
		//we use Cstyle arrays which is why we need space for the terminator

#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "dimensionality " << (int) dimensionality << "\n";
		cout << "memory " << (int) memory << "\n";
		cout << "terminator " << (int) terminator << "\n";
		cout << "encoding " << (int) encoding << "\n";
		cout << "ndims " << val.size() << "\n";
		cout << "sdims " << max_characters << "\n";
#endif
		//MK::this function may only be called from nexus_write_attribute as it adds a typed and named attribute
		//to an (assumed existent) dataset referred to by dsetid inside an RDWR-opened HDF5 file referred to by fileid
		//for a specific type with the name in the keyword and the value in the val variable
		//see also an example for string arrays here https://svn.ssec.wisc.edu/repos/geoffc/C/HDF5/Examples_by_API/h5ex_t_stringatt.c

		int status = 0;
		adtypid = H5Tcopy(H5T_C_S1);
		if ( adtypid != H5I_INVALID_HID ) {
			if ( memory == MYHDF5_IS_VARIABLE ) {
				status = H5Tset_size(adtypid, H5T_VARIABLE);
			}
			else { //MYHDF5_IS_FIXED
				status = H5Tset_size(adtypid, max_characters);
			}
			if ( status >= 0 ) {
				switch(terminator)
				{
					case MYHDF5_NULLTERM: { status = H5Tset_strpad(adtypid, H5T_STR_NULLTERM); break; }
					case MYHDF5_NULLPAD: { status = H5Tset_strpad(adtypid, H5T_STR_NULLPAD); break; }
					case MYHDF5_SPACEPAD: { status = H5Tset_strpad(adtypid, H5T_STR_SPACEPAD); break; }
					default: { status = -1; break; }
				}
				if ( status >= 0 ) {
					switch(encoding)
					{
						case MYHDF5_ASCII: { status = H5Tset_cset(adtypid, H5T_CSET_ASCII); break; }
						case MYHDF5_UTF8: { status = H5Tset_cset(adtypid, H5T_CSET_UTF8); break; }
						default: { status = -1; break; }
					}
					if ( status >= 0 ) {
						//dimensionality == MYHDF5_LIST_SINGLE or LIST_LONG
						int rank = 1; hsize_t dims[1] = { val.size() };
						aspcid = H5Screate_simple( rank, dims, NULL );
						if ( aspcid != H5I_INVALID_HID ) {
							attrid = H5Acreate2(dsetid, keyword.c_str(), adtypid, aspcid, H5P_DEFAULT, H5P_DEFAULT);
							if ( attrid >= 0 ) {
								if ( dimensionality == MYHDF5_LIST_SINGLE ) {
									if ( memory == MYHDF5_IS_VARIABLE ) {
										if ( H5Awrite(attrid, adtypid, val.data() ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
											cout << "Writing attribute " << keyword << " success" << "\n";
#endif
											return nexus_attribute_close( MYHDF5_SUCCESS );
										}
										cerr << "Writing attribute " << keyword << " failed!" << "\n";
										return nexus_attribute_close( MYHDF5_ATTRWRITE_FAILED );
									}
									else { //MYHDF5_IS_FIXED
										char* adata = new char[max_characters];
										if ( adata != NULL ) {
											strcpy(adata, val[0].c_str());
											if ( H5Awrite(attrid, adtypid, adata) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
												cout << "Writing attribute " << keyword << " success" << "\n";
#endif
												delete[] adata; adata = NULL;
												return nexus_attribute_close( MYHDF5_SUCCESS );
											}
											cerr << "Writing attribute " << keyword << " failed!" << "\n";
											delete[] adata; adata = NULL;
											return nexus_attribute_close( MYHDF5_ATTRWRITE_FAILED );
										}
										cerr << "Writing attribute " << keyword << " failed!" << "\n";
										delete[] adata; adata = NULL;
										return nexus_attribute_close( MYHDF5_ATTRWRITE_FAILED );
									}
								} //list_single
								else { //MYHDF5_LIST_LONG
									if ( memory == MYHDF5_IS_VARIABLE ) {
										status = 0;
										char** adata = new char*[val.size()];
										if ( adata != NULL ) {
											for ( size_t i=0; i < val.size(); i++ ) {
												adata[i] = new char[val[i].length()+1];
												if ( adata[i] != NULL ) {
													strcpy(adata[i], val[i].c_str());
												}
												else {
													status = -1;
												}
											}
											if ( status >= 0 ) {
												if ( H5Awrite(attrid, adtypid, adata ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
													cout << "Writing attribute " << keyword << " success" << "\n";
#endif
													for ( size_t i=0; i < val.size(); i++) {
														delete[] adata[i]; adata[i] = NULL;
													}
													delete[] adata; adata = NULL;
													return nexus_attribute_close( MYHDF5_SUCCESS );
												}
												cerr << "Writing attribute " << keyword << " failed!" << "\n";
												for ( size_t i=0; i < val.size(); i++) {
													delete[] adata[i]; adata[i] = NULL;
												}
												delete[] adata; adata = NULL;
												return nexus_attribute_close( MYHDF5_ATTRWRITE_FAILED );
											}
											cerr << "Writing attribute " << keyword << " failed!" << "\n";
											for ( size_t i=0; i < val.size(); i++) {
												delete[] adata[i]; adata[i] = NULL;
											}
											delete[] adata; adata = NULL;
											return nexus_attribute_close( MYHDF5_ATTRWRITE_FAILED );
										}
										delete[] adata; adata = NULL;
										return nexus_attribute_close( MYHDF5_ATTRWRITE_FAILED );
									}
									else { //MYHDF_IS_FIXED
										char* adata = new char[val.size() * max_characters];
										if ( adata != NULL ) {
											for( size_t i = 0; i < val.size(); i++ ) {
												for( size_t j = 0; j < val.at(i).length(); j++ ) {
													adata[i*max_characters+j] = val[i][j];
												}
												for( size_t j = val.at(i).length(); j < max_characters; j++ ) {
													switch(terminator)
													{
														case MYHDF5_NULLTERM: { adata[i*max_characters+j] = '\0'; break; }
														case MYHDF5_NULLPAD: { adata[i*max_characters+j] = 0; break; }
														case MYHDF5_SPACEPAD: { adata[i*max_characters+j] = ' '; break; }
														default: { break; }
													}
												}
											}
											if ( H5Awrite(attrid, adtypid, adata) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
												cout << "Writing attribute " << keyword << " success" << "\n";
#endif
												delete[] adata; adata = NULL;
												return nexus_attribute_close( MYHDF5_SUCCESS );
											}
											cerr << "Writing attribute " << keyword << " failed!" << "\n";
											delete[] adata; adata = NULL;
											return nexus_attribute_close( MYHDF5_ATTRWRITE_FAILED );
										}
										cerr << "Writing attribute " << keyword << " failed!" << "\n";
										delete[] adata; adata = NULL;
										return nexus_attribute_close( MYHDF5_ATTRWRITE_FAILED );
									}
								} //list_long
								cerr << "Writing attribute " << keyword << " failed!" << "\n";
								return nexus_attribute_close( MYHDF5_ATTRWRITE_FAILED );
							}
							cerr << "Creating attribute " << keyword << " failed!" << "\n";
							return nexus_attribute_close( MYHDF5_ATTRCREATE_FAILED );
						}
						cerr << "Dspacing attribute " << keyword << " failed!" << "\n";
						return nexus_attribute_close( MYHDF5_ATTRCREATE_FAILED );
					}
					cerr << "Formatting cset attribute " << keyword << " failed!" << "\n";
					return nexus_attribute_close( MYHDF5_FAILED );
				}
				cerr << "Formatting str_strpad attribute " << keyword << " failed!" << "\n";
				return nexus_attribute_close( MYHDF5_FAILED );
			}
			cerr << "Formatting size attribute " << keyword << " failed!" << "\n";
			return nexus_attribute_close( MYHDF5_FAILED );
		}
		cerr << "Creating datatype attribute " << keyword << " failed!" << "\n";
		return nexus_attribute_close( MYHDF5_FAILED );
	}
	//empty lists are not written out and dont create problems
	return MYHDF5_SUCCESS;
}


//high-level
int HdfFiveSeqHdl::nexus_write_attributes( const string loc_name, ioAttributes const & attrs )
{
	if ( nexus_write_basic_scalar_attributes( loc_name, attrs ) != MYHDF5_SUCCESS )
		return MYHDF5_FAILED;
	if ( nexus_write_string_scalar_attributes( loc_name, attrs ) != MYHDF5_SUCCESS )
		return MYHDF5_FAILED;
	if ( nexus_write_basic_array_attributes( loc_name, attrs ) != MYHDF5_SUCCESS )
		return MYHDF5_FAILED;
	if ( nexus_write_string_array_attributes( loc_name, attrs ) != MYHDF5_SUCCESS )
		return MYHDF5_FAILED;

	return MYHDF5_SUCCESS;
}


int HdfFiveSeqHdl::nexus_write_basic_scalar_attributes( const string loc_name, ioAttributes const & attrs )
{
	//MK::this function must only be called within one of the nexus_write functions as it
	//assumes/reuses open(ed) resources on an existent file and existent dataset to reduce costs
	//for checks and unnecessary interaction with the HDF5 library !
	//assume fileid, dset are set to valid resources on an existent dataset entry loc_name in an open file

	htri_t status = 0;  //<0 (failure), == 0 (false), == 1 (true)
	//https://docs.hdfgroup.org/hdf5/develop/_h5public_8h.html#aa8f6c28736dbd0f18388c67911d38aca
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
	cout << "Writing attributes to loc_name " << loc_name << "\n";
#endif

//scalar attributes
	for( auto it = attrs.u8.begin(); it != attrs.u8.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status == 0 ) {
			if ( nexus_write_attribute_value( it->first, it->second ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding u8 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			} //do not return something because we want to add further attributes
			else {
				cerr << "Adding u8 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else if ( status > 0 ) {
			cerr << "Adding u8 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else { //status < 0
			cerr << "Adding u8 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.i8.begin(); it != attrs.i8.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status == 0 ) {
			if ( nexus_write_attribute_value( it->first, it->second ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding i8 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding i8 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else if ( status > 0 ) {
			cerr << "Adding i8 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else {
			cerr << "Adding i8 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.u16.begin(); it != attrs.u16.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status == 0 ) {
			if ( nexus_write_attribute_value( it->first, it->second ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding u16 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding u16 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else if ( status > 0 ) {
			cerr << "Adding u16 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else {
			cerr << "Adding u16 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.i16.begin(); it != attrs.i16.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status == 0 ) {
			if ( nexus_write_attribute_value( it->first, it->second ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding i16 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding i16 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else if ( status > 0 ) {
			cerr << "Adding i16 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else {
			cerr << "Adding i16 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.u32.begin(); it != attrs.u32.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status == 0 ) {
			if ( nexus_write_attribute_value( it->first, it->second ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding u32 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding u32 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else if ( status > 0 ) {
			cerr << "Adding u32 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else {
			cerr << "Adding u32 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.i32.begin(); it != attrs.i32.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status == 0 ) {
			if ( nexus_write_attribute_value( it->first, it->second ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding i32 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding i32 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else if ( status > 0 ) {
			cerr << "Adding i32 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else {
			cerr << "Adding i32 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.u64.begin(); it != attrs.u64.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status == 0 ) {
			if ( nexus_write_attribute_value( it->first, it->second ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding u64 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding u64 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else if ( status > 0 ) {
			cerr << "Adding u64 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else {
			cerr << "Adding u64 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.i64.begin(); it != attrs.i64.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status == 0 ) {
			vector<long> value_input = { it->second };
			if ( nexus_write_attribute_value( it->first, it->second ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding i64 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding i64 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else if ( status > 0 ) {
			cerr << "Adding i64 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else {
			cerr << "Adding i64 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.f32.begin(); it != attrs.f32.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status == 0 ) {
			if ( nexus_write_attribute_value( it->first, it->second ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding f32 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding f32 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else if ( status > 0 ) {
			cerr << "Adding f32 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else {
			cerr << "Adding f32 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.f64.begin(); it != attrs.f64.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status == 0 ) {
			if ( nexus_write_attribute_value( it->first, it->second ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding f64 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding f64 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else if ( status > 0 ) {
			cerr << "Adding f64 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else {
			cerr << "Adding f64 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	return MYHDF5_SUCCESS;
}


int HdfFiveSeqHdl::nexus_write_string_scalar_attributes( const string loc_name, ioAttributes const & attrs )
{
	//MK::this function must only be called within one of the nexus_write functions as it
	//assumes/reuses open(ed) resources on an existent file and existent dataset to reduce costs
	//for checks and unnecessary interaction with the HDF5 library !
	//assume fileid, dset are set to valid resources on an existent dataset entry loc_name in an open file

	htri_t status = 0;  //<0 (failure), == 0 (false), == 1 (true)
	//scalar string as an attribute
	for( auto it = attrs.str.begin(); it != attrs.str.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status > 0 ) {
			cerr << "Adding str attribute " << it->first << ", " << it->second << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else if ( status == 0 ) {
			if ( nexus_write_attribute_value( it->first, it->second ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding str attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding str attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else {
			cerr << "Adding str attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	return MYHDF5_SUCCESS;
}


int HdfFiveSeqHdl::nexus_write_basic_array_attributes( const string loc_name, ioAttributes const & attrs )
{
	//MK::this function must only be called within one of the nexus_write functions as it
	//assumes/reuses open(ed) resources on an existent file and existent dataset to reduce costs
	//for checks and unnecessary interaction with the HDF5 library !
	//assume fileid, dset are set to valid resources on an existent dataset entry loc_name in an open file

	htri_t status = 0;  //<0 (failure), == 0 (false), == 1 (true)
	//basic type array attributes
	for( auto it = attrs.u8_arr.begin(); it != attrs.u8_arr.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status > 0 ) {
			cerr << "Adding u8 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else if ( status == 0 ) {
			vector<unsigned char> value_input = it->second;
			if ( nexus_write_attribute_value( it->first, value_input ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding u8 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding u8 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else {
			cerr << "Adding u8 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.i8_arr.begin(); it != attrs.i8_arr.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status > 0 ) {
			cerr << "Adding i8 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else if ( status == 0 ) {
			vector<char> value_input = it->second;
			if ( nexus_write_attribute_value( it->first, value_input ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding i8 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding i8 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else {
			cerr << "Adding i8 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.u16_arr.begin(); it != attrs.u16_arr.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status > 0 ) {
			cerr << "Adding u16 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else if ( status == 0 ) {
			vector<unsigned short> value_input = it->second;
			if ( nexus_write_attribute_value( it->first, value_input ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding u16 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding u16 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else {
			cerr << "Adding u16 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.i16_arr.begin(); it != attrs.i16_arr.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status > 0 ) {
			cerr << "Adding i16 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else if ( status == 0 ) {
			vector<short> value_input = it->second;
			if ( nexus_write_attribute_value( it->first, value_input ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding i16 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding i16 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else {
			cerr << "Adding i16 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.u32_arr.begin(); it != attrs.u32_arr.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status > 0 ) {
			cerr << "Adding u32 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else if ( status == 0 ) {
			vector<unsigned int> value_input = it->second;
			if ( nexus_write_attribute_value( it->first, value_input ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding u32 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding u32 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else {
			cerr << "Adding u32 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.i32_arr.begin(); it != attrs.i32_arr.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status > 0 ) {
			cerr << "Adding i32 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else if ( status == 0 ) {
			vector<int> value_input = it->second;
			if ( nexus_write_attribute_value( it->first, value_input ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding i32 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding i32 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else {
			cerr << "Adding i32 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.u64_arr.begin(); it != attrs.u64_arr.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status > 0 ) {
			cerr << "Adding u64 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else if ( status == 0 ) {
			vector<unsigned long> value_input = it->second;
			if ( nexus_write_attribute_value( it->first, value_input ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding u64 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding u64 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else {
			cerr << "Adding u64 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.i64_arr.begin(); it != attrs.i64_arr.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status > 0 ) {
			cerr << "Adding i64 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else if ( status == 0 ) {
			vector<long> value_input = it->second;
			if ( nexus_write_attribute_value( it->first, value_input ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding i64 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding i64 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else {
			cerr << "Adding i64 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.f32_arr.begin(); it != attrs.f32_arr.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status > 0 ) {
			cerr << "Adding f32 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else if ( status == 0 ) {
			vector<float> value_input = it->second;
			if ( nexus_write_attribute_value( it->first, value_input ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding f32 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding f32 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else {
			cerr << "Adding f32 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}

	for( auto it = attrs.f64_arr.begin(); it != attrs.f64_arr.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status > 0 ) {
			cerr << "Adding f64 attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else if ( status == 0 ) {
			vector<double> value_input = it->second;
			if ( nexus_write_attribute_value( it->first, value_input ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding f64 attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding f64 attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else {
			cerr << "Adding f64 attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}
	return MYHDF5_SUCCESS;
}


int HdfFiveSeqHdl::nexus_write_string_array_attributes( const string loc_name, ioAttributes const & attrs )
{
	//MK::this function must only be called within one of the nexus_write functions as it
	//assumes/reuses open(ed) resources on an existent file and existent dataset to reduce costs
	//for checks and unnecessary interaction with the HDF5 library !
	//assume fileid, dset are set to valid resources on an existent dataset entry loc_name in an open file
	htri_t status = 0;  //<0 (failure), == 0 (false), == 1 (true)
	//string array attributes
	for( auto it = attrs.str_arr.begin(); it != attrs.str_arr.end(); it++ ) {
		status = H5Aexists_by_name( fileid, loc_name.c_str(), it->first.c_str(), 0 );
		if ( status > 0 ) {
			cerr << "Adding str attribute " << it->first << " to loc_name " << loc_name << " rejected because an attribute with such a name exists !" << "\n";
		}
		else if ( status == 0 ) {
			vector<string> value_input = it->second;
			if ( nexus_write_attribute_value( it->first, value_input ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
				cout << "Adding str attribute " << it->first << " to loc_name " << loc_name << " success" << "\n";
#endif
			}
			else {
				cerr << "Adding str attribute " << it->first << " to loc_name " << loc_name << " failed!" << "\n";
				return MYHDF5_ATTRWRITE_FAILED;
			}
		}
		else {
			cerr << "Adding str attribute " << it->first << " to loc_name " << loc_name << " failed, upon existence check !" << "\n";
			return MYHDF5_ATTRACCESS_FAILED;
		}
	}
	return MYHDF5_SUCCESS;
}


vector<string> HdfFiveSeqHdl::split_absolute_path( const string h5_absolute_path )
{
	vector<string> retval = vector<string>();
	char delimiter = '/';
	string fwslash = "/";

#ifdef MYHDF5_VERBOSE_LEVEL_ONE
	cout << "__" << h5_absolute_path << "__" << "\n";
#endif

#ifdef MYHDF5_VERBOSE_LEVEL_TWO
	cout << "h5_absolute_path.length() " << h5_absolute_path.length() << "\n";
#endif
	if ( h5_absolute_path.length() > 0 ) {
		size_t idx_s = 0;
		size_t idx_e = h5_absolute_path.length();
		//remove eventually preceeding and/or trailing fwslash
		if (h5_absolute_path.rfind(delimiter, 0) != string::npos) {
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
			cout << "/ at start" << "\n";
#endif
			idx_s++;
		}
		if (h5_absolute_path.compare(h5_absolute_path.length()-1, 1, fwslash) == 0) {
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		  cout << "/ at end" << "\n";
#endif
		  idx_e--;
		}
		stringstream sstream("");
		for ( size_t i = idx_s; i < idx_e; ++i) {
			sstream << h5_absolute_path.at(i);
		}
		string tmp = "";
		while( getline(sstream, tmp, delimiter) ) {
			retval.push_back(tmp);
		}
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "retval.size() " << retval.size() << "\n";
#endif
		/*
		cout << "debugging parts" << "\n";
		for( size_t i = 0; i < parts.size(); i++ ) {
			cout << parts[i] << "\n";
		}*/
	}
	return retval;
}


int HdfFiveSeqHdl::nexus_write_group( const string grpnm, ioAttributes const & attrs )
{
	if ( nexus_open( H5F_ACC_RDWR ) == H5I_INVALID_HID ) {
		cerr << "Opening file during write group attempt failed!" << "\n";
		return MYHDF5_FOPEN_FAILED;
	}
	//now the file is open otherwise we wouldnt have gotten here

	//check if the path in the open file pointed to by fileid exists
	vector<string> level_by_level = split_absolute_path( grpnm );

	//start at the file, root
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
	cout << "Checking a chain of hopefully existing links..." << "\n";
#endif
	hid_t curr_loc_id = fileid;
	string curr_group = "";

	//##MK::https://gist.github.com/jzrake/3025642
	for( size_t i = 0; i < level_by_level.size(); i++ ) {
		curr_group = level_by_level[i]; //nodes_visited += "/" + graph_nodes[i];
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "Currently visiting " << curr_group << "\n";
#endif
		htri_t status = H5Lexists( curr_loc_id, curr_group.c_str(), H5P_DEFAULT );
		if( status > 0 ) { //exists
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
			cout << curr_group << " --> exists" << "\n";
#endif
			hid_t next_loc_id = H5Gopen( curr_loc_id, curr_group.c_str(), H5P_DEFAULT );
			if ( next_loc_id == H5I_INVALID_HID ) {
				if ( curr_loc_id != fileid ) {
					if ( H5Gclose( curr_loc_id ) >= 0 ) {
						//nop
					}
					else {
						if ( H5Fclose( fileid ) >= 0 ) {
							cerr << "Successful H5Fclose after a failed H5Gclose while preparing walk after " << curr_group << "\n";
							return MYHDF5_GCLOSE_FAILED;
						}
						cerr << "Failed H5Fclose after a failed H5Gclose while preparing walk after " << curr_group << "\n";
						return MYHDF5_FCLOSE_FAILED;
					}
				}
				if ( H5Fclose( fileid ) >= 0 ) {
					cerr << "Successful H5Fclose after failed H5Gopen on " << curr_group << "\n";
					return MYHDF5_GOPEN_FAILED;
				}
				cerr << "Failed H5Fclose after failed H5Gopen on " << curr_group << "\n";
				return MYHDF5_FCLOSE_FAILED;
			}
			else {
				if ( i == (level_by_level.size() - 1) ) {
					//open
					dsetid = next_loc_id; //this line is a hack that should be implemented cleaner
					//nexus_write_attributes do not open resources but assume that there is an opened
					//loc_id with a valid location id handler dsetid against which attributes are written
					//usually the locations are datasets but here we wish these to be groups and thus
					//we need to pass the class this information, the attribute write function do
					//not currently and must not close the resources pointed to my dsetid!only the
					//they only close the locally created handles for the attributes
					int attr_status = nexus_write_attributes( grpnm, attrs );
					if ( attr_status == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
						cout << "Adding attributes to " << grpnm << " was successful" << "\n";
#endif
					}
					else {
						cerr << "Adding attributes to " << grpnm << " failed!" << "\n";
					}
				}
				if ( curr_loc_id != fileid ) {
					if ( H5Gclose( curr_loc_id ) >= 0 ) {
						//nop
					}
					else {
						if ( H5Fclose( fileid ) >= 0 ) {
							cerr << "Successful H5Fclose after a failed H5Gclose while preparing walk after >0 " << curr_group << "\n";
							return MYHDF5_GCLOSE_FAILED;
						}
						cerr << "Failed H5Fclose after a failed H5Gclose while preparing walk after >0 " << curr_group << "\n";
						return MYHDF5_FCLOSE_FAILED;
					}
				}
				curr_loc_id = next_loc_id;
				continue;
			}
		}
		else if ( status == 0 ) {
			//no error, library states that group under nodes_visited does not exist, so create it
			hid_t next_loc_id = H5Gcreate2(curr_loc_id, curr_group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			if ( next_loc_id == H5I_INVALID_HID ) {
				if ( curr_loc_id != fileid ) {
					if ( H5Gclose( curr_loc_id ) >= 0 ) {
						//nop
					}
					else {
						if ( H5Fclose( fileid ) >= 0 ) {
							cerr << "Successful H5Fclose after a failed H5Gclose while preparing walk after ==0 " << curr_group << "\n";
							return MYHDF5_GCLOSE_FAILED;
						}
						cerr << "Failed H5Fclose after a failed H5Gclose while preparing walk after ==0 " << curr_group << "\n";
						return MYHDF5_FCLOSE_FAILED;
					}
				}
				if ( H5Fclose( fileid ) >= 0 ) {
					cerr << "Successful H5Fclose after failed creating group ==0 " << curr_group << "\n";
					return MYHDF5_GCREATE_FAILED;
				}
				cerr << "Failed H5Fclose after failed create group ==0 " << curr_group << "\n";
				return MYHDF5_FCLOSE_FAILED;
			}
			else {
				if ( i == (level_by_level.size() - 1) ) {
					dsetid = next_loc_id;
					int attr_status = nexus_write_attributes( grpnm, attrs );
					if ( attr_status == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
						cout << "Adding attributes to " << grpnm << " was successful" << "\n";
#endif
					}
					else {
						cerr << "Adding attributes to " << grpnm << " failed!" << "\n";
					}
				}

				if ( curr_loc_id != fileid ) {
					if ( H5Gclose( curr_loc_id ) >= 0 ) {
						//nop
					}
					else {
						if ( H5Fclose( fileid ) >= 0 ) {
							cerr << "Successful H5Fclose after a failed H5Gclose after creating sub-group in " << curr_group << "\n";
							return MYHDF5_GCLOSE_FAILED;
						}
						cerr << "Failed H5Fclose after a failed H5Gclose after creating sub-group in " << curr_group << "\n";
						return MYHDF5_FCLOSE_FAILED;
					}
				}
				curr_loc_id = next_loc_id;
				continue;
			}
		}
		else { //status < 0
			cerr << "htri_t library error" << "\n";
			if ( curr_loc_id != fileid ) {
				if ( H5Gclose( curr_loc_id ) >= 0 ) {
					//nop
				}
				else {
					if ( H5Fclose( fileid ) >= 0 ) {
						cerr << "Successful H5Fclose after a failed H5Gclose when link did not exist " << curr_group << "\n";
						return MYHDF5_GCLOSE_FAILED;
					}
					cerr << "Failed H5Fclose after a failed H5Gclose when link did not exist " << curr_group << "\n";
					return MYHDF5_FCLOSE_FAILED;
				}
			}
			break;
		}
	}

	//##MK::add the annotations in the last group

	//success is achieved only when you have not previous run into a failure
	//such as having to stop the walk because of a missing path
	if ( H5Gclose( curr_loc_id ) < 0 ) {
		cerr << "Failed H5Gclose during cleanup" << "\n";
		return MYHDF5_GCLOSE_FAILED;
	}
	if ( H5Fclose( fileid ) < 0 ) {
		cerr << "Failed H5Fclose during cleanup" << "\n";
		return MYHDF5_FCLOSE_FAILED;
	}

	return MYHDF5_SUCCESS;
}


/*
//deprecated::prior 0.4.1 tools wrote even scalars always into a list as seen below
template<typename T> int HdfFiveSeqHdl::nexus_write( const string dsnm, T & val, ioAttributes const & attrs )
{
	vector<T> tmp = vector<T>( 1, val );
	return nexus_write( dsnm, io_info(), tmp, attrs );
}
*/

template<typename T> int HdfFiveSeqHdl::nexus_write( const string dsnm, const T val, ioAttributes const & attrs )
{
	//write scalar basic type, never compression for scalars
	dtypid = H5I_INVALID_HID;
	string msg = "";
	if constexpr(std::is_same<T, unsigned char>::value) { dtypid = H5T_STD_U8LE; msg = "u8"; }
	else if constexpr(std::is_same<T, char>::value) { dtypid = H5T_STD_I8LE; msg = "i8"; }
	else if constexpr(std::is_same<T, unsigned short>::value) { dtypid = H5T_STD_U16LE; msg = "u16"; }
	else if constexpr(std::is_same<T, short>::value) { dtypid = H5T_STD_I16LE; msg = "i16"; }
	else if constexpr(std::is_same<T, unsigned int>::value) { dtypid = H5T_STD_U32LE; msg = "u32"; }
	else if constexpr(std::is_same<T, int>::value) { dtypid = H5T_STD_I32LE; msg = "i32"; }
	else if constexpr(std::is_same<T, unsigned long>::value) { dtypid = H5T_STD_U64LE; msg = "u64"; }
	else if constexpr(std::is_same<T, long>::value) { dtypid = H5T_STD_I64LE; msg = "i64"; }
	else if constexpr(std::is_same<T, float>::value) { dtypid = H5T_IEEE_F32LE; msg = "f32"; }
	else if constexpr(std::is_same<T, double>::value) { dtypid = H5T_IEEE_F64LE; msg = "f64"; }
	else { /*nop*/ }
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "dtypid " << msg << " " << dtypid << "\n";
#endif

	if ( nexus_open( H5F_ACC_RDWR ) != H5I_INVALID_HID ) {
		//if ( nexus_path_exists( dsnm, true ) == false ) {
		dspcid = H5Screate( H5S_SCALAR );
		if ( dspcid != H5I_INVALID_HID ) {
			dsetid = H5Dcreate2( fileid, dsnm.c_str(), dtypid, dspcid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
			if ( dsetid != H5I_INVALID_HID ) {
				if ( H5Dwrite( dsetid, dtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &val ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
					cout << "Writing " << dsnm << " basic, scalar, success" << "\n";
#endif
					if ( nexus_write_attributes( dsnm, attrs ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
						cout << "Writing " << dsnm << " basic, scalar, attributes, success" << "\n";
#endif
						return nexus_close( MYHDF5_SUCCESS );
					}
					cerr << "Writing " << dsnm << " basic, scalar, attributes, scalar failed !" << "\n";
					return nexus_close( MYHDF5_ATTRWCOLL_FAILED );
				}
				cerr << "Writing " << dsnm << " basic, scalar, failed write !" << "\n";
				return nexus_close( MYHDF5_DSETWRITE_FAILED );
			}
			cerr << "Creating " << dsnm << " basic, scalar, failed dset !" << "\n";
			return nexus_close( MYHDF5_DSETOPEN_FAILED );
		}
		cerr << "Creating " << dsnm << " basic, scalar, failed dspc !" << "\n";
		return nexus_close( MYHDF5_DSPCOPEN_FAILED );
	}
	cerr << "Opening file " << h5resultsfn << " failed !" << "\n";
	return nexus_close( MYHDF5_FOPEN_FAILED );
}

template int HdfFiveSeqHdl::nexus_write( const string dsnm, const unsigned char val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, const char val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, const unsigned short val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, const short val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, const unsigned int val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, const int val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, const unsigned long val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, const long val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, const float val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, const double val, ioAttributes const & attrs );


template<typename T> int HdfFiveSeqHdl::nexus_write( const string dsnm, io_info const & ifo, vector<T> & val, ioAttributes const & attrs )
{
	//open file in read/write mode
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
	cout << "nexus_write array version val.size() " << val.size() << "\n";
#endif
	if ( val.size() > 0 ) {
		size_t n_values = 1;
		for( auto it = ifo.shape.begin(); it != ifo.shape.end(); it++ ) {
			n_values *= *it;
		}
		if ( ifo.is_valid == true && val.size() != 0 && val.size() >= n_values ) {
			dtypid = H5I_INVALID_HID;
			string msg = "";
			if constexpr(std::is_same<T, unsigned char>::value) { dtypid = H5T_STD_U8LE; msg = "u8"; }
			else if constexpr(std::is_same<T, char>::value) { dtypid = H5T_STD_I8LE; msg = "i8"; }
			else if constexpr(std::is_same<T, unsigned short>::value) { dtypid = H5T_STD_U16LE; msg = "u16"; }
			else if constexpr(std::is_same<T, short>::value) { dtypid = H5T_STD_I16LE; msg = "i16"; }
			else if constexpr(std::is_same<T, unsigned int>::value) { dtypid = H5T_STD_U32LE; msg = "u32"; }
			else if constexpr(std::is_same<T, int>::value) { dtypid = H5T_STD_I32LE; msg = "i32"; }
			else if constexpr(std::is_same<T, unsigned long>::value) { dtypid = H5T_STD_U64LE; msg = "u64"; }
			else if constexpr(std::is_same<T, long>::value) { dtypid = H5T_STD_I64LE; msg = "i64"; }
			else if constexpr(std::is_same<T, float>::value) { dtypid = H5T_IEEE_F32LE; msg = "f32"; }
			else if constexpr(std::is_same<T, double>::value) { dtypid = H5T_IEEE_F64LE; msg = "f64"; }
			else { /*nop*/ }
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "n_values " << n_values << "\n";
		cout << "dtypid " << dtypid << " " << msg << "\n";
#endif

			if ( nexus_open( H5F_ACC_RDWR ) != H5I_INVALID_HID ) {
				//if ( nexus_path_exists( dsnm, true ) == false ) {
					if ( ifo.chunk.size() == 0 ) { //no chunking but contiguous layout
						if ( ifo.shape.size() == 1 ) { //1d, vector
							int rank = 1;
							hsize_t dims[1] = { ifo.shape[0] };
							hsize_t maxdims[1] = { ifo.shape[0] };
							dspcid = H5Screate_simple( rank, dims, maxdims );
							if ( dspcid != H5I_INVALID_HID ) {
								dsetid = H5Dcreate2( fileid, dsnm.c_str(), dtypid, dspcid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
								if ( dsetid != H5I_INVALID_HID ) {
									if ( H5Dwrite( dsetid, dtypid, dspcid, dspcid, H5P_DEFAULT, val.data() ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
										cout << "Writing " << dsnm << " contiguous, 1d success" << "\n";
#endif

										int attr_status = nexus_write_attributes( dsnm, attrs );
										if ( attr_status == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
											cout << "Writing " << dsnm << " contiguous, 1d attributes success !" << "\n";
#endif
										}
										else {
											cerr << "Writing " << dsnm << " contiguous, 1d attributes failed " << attr_status << "\n";
										}
										return nexus_close( MYHDF5_SUCCESS );
									}
									cerr << "Writing " << dsnm << " contiguous, 1d failed !" << "\n";
									return nexus_close( MYHDF5_DSETWRITE_FAILED );
								}
								cerr << "Creating " << dsnm << " contiguous, 1d dset failed !" << "\n";
								return nexus_close( MYHDF5_DSETOPEN_FAILED );
							}
							cerr << "Creating " << dsnm << " contiguous, 1d dspc failed !" << "\n";
							return nexus_close( MYHDF5_DSPCOPEN_FAILED );
						}
						else if ( ifo.shape.size() == 2 ) { //2d, matrix
							int rank = 2;
							hsize_t dims[2] = { ifo.shape[0], ifo.shape[1] };
							hsize_t maxdims[2] = { ifo.shape[0], ifo.shape[1] };

							hsize_t offset[2] = { 0, 0 }; //hyperslab design
							hsize_t count[2] = { ifo.shape[0], ifo.shape[1] };
							hsize_t stride[2] = { 1, 1 };
							hsize_t block[2] = { 1, 1 };

							dspcid = H5Screate_simple( rank, dims, maxdims );
							if ( dspcid != H5I_INVALID_HID ) {
								dsetid = H5Dcreate2( fileid, dsnm.c_str(), dtypid, dspcid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
								if ( dsetid != H5I_INVALID_HID ) {
									if ( H5Sselect_hyperslab(dspcid, H5S_SELECT_SET, offset, stride, count, block) >= 0 ) {
										if ( H5Dwrite( dsetid, dtypid, dspcid, dspcid, H5P_DEFAULT, val.data() ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
											cout << "Writing " << dsnm << " contiguous, 2d success" << "\n";
#endif

											int attr_status = nexus_write_attributes( dsnm, attrs );
											if ( attr_status == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
												cout << "Writing " << dsnm << " contiguous, 2d attributes success !" << "\n";
#endif
											}
											else {
												cerr << "Writing " << dsnm << " contiguous, 2d attributes failed " << attr_status << "\n";
											}
											return nexus_close( MYHDF5_SUCCESS );
										}
										cerr << "Writing " << dsnm << " contiguous, 2d failed !" << "\n";
										return nexus_close( MYHDF5_DSETWRITE_FAILED );
									}
									cerr << "Creating " << dsnm << " contiguous, hslab 2d failed !" << "\n";
									return nexus_close( MYHDF5_DSETWRITE_FAILED );
								}
								cerr << "Creating " << dsnm << " contiguous, 2d dset failed !" << "\n";
								return nexus_close( MYHDF5_DSETOPEN_FAILED );
							}
							cerr << "Creating " << dsnm << " contiguous, 2d dspc failed !" << "\n";
							return nexus_close( MYHDF5_DSPCOPEN_FAILED );
						}
						cerr << "Creating " << dsnm << " failed, because d > 2 arrays not supported !" << "\n";
						return nexus_close( MYHDF5_NOT_IMPLEMENTED );
					} //done with contiguous
					else { // chunked layout, eventually using compression
						if ( ifo.shape.size() == 1 ) { //1d, vector
							//uncomment these variables to use SZIP compression
							//unsigned szip_options_mask;
							//unsigned szip_pixels_per_block;
							int rank = 1;
							hsize_t dims[1] = { ifo.shape[0] };
							hsize_t maxdims[1] = { ifo.shape[0] };
							dspcid = H5Screate_simple( rank, dims, maxdims );
							if ( dspcid != H5I_INVALID_HID ) {
								plistid = H5Pcreate( H5P_DATASET_CREATE );
								if ( plistid != H5I_INVALID_HID ) {
									int chunk_rank = 1;
									hsize_t chunk_dims[1] = { ifo.chunk[0] };
									if ( H5Pset_chunk( plistid, chunk_rank, chunk_dims ) >= 0 ) {
										if ( ifo.compression == MYHDF5_COMPRESSION_GZIP ) {
											if ( H5Pset_deflate( plistid, (unsigned) ifo.compression_opts ) >= 0 ) {
												//set ZLIB / DEFLATE Compression using compression level 6.
												//szip_options_mask = H5_SZIP_NN_OPTION_MASK;
												//szip_pixels_per_block = 16;
												//status = H5Pset_szip (plistid, szip_options_mask, szip_pixels_per_block);
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
												cout << "Setting deflate " << dsnm << " chunked, 1d" << "\n";
#endif
											}
											else {
												cerr << "Creating " << dsnm << " chunked, 1d set deflate failed !" << "\n";
												return nexus_close( MYHDF5_PLSTACCESS_FAILED );
											}
										}
										dsetid = H5Dcreate2( fileid, dsnm.c_str(), dtypid, dspcid, H5P_DEFAULT, plistid, H5P_DEFAULT );
										if ( dsetid != H5I_INVALID_HID ) {
											if ( H5Dwrite( dsetid, dtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, val.data() ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
												cout << "Writing " << dsnm << " chunked, 1d success" << "\n";
#endif

												int attr_status = nexus_write_attributes( dsnm, attrs );
												if ( attr_status == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
													cout << "Writing " << dsnm << " contiguous, 2d attributes success !" << "\n";
#endif
												}
												else {
													cerr << "Writing " << dsnm << " contiguous, 2d attributes failed " << attr_status << "\n";
												}
												return nexus_close( MYHDF5_SUCCESS );
											}
											cerr << "Writing " << dsnm << " chunked, 1d create failed !" << "\n";
											return nexus_close( MYHDF5_DSETWRITE_FAILED );
										}
										cerr << "Writing " << dsnm << " chunked, 1d create failed !" << "\n";
										return nexus_close( MYHDF5_DSETWRITE_FAILED );
									}
									cerr << "Creating " << dsnm << " chunked, 1d set chunk failed !" << "\n";
									return nexus_close( MYHDF5_DSETWRITE_FAILED );
								}
								cerr << "Creating " << dsnm << " chunked, 1d filter create failed !" << "\n";
								return nexus_close( MYHDF5_PLSTOPEN_FAILED );
							}
							cerr << "Creating " << dsnm << " chunked, 1d set dspc failed !" << "\n";
							return nexus_close( MYHDF5_DSPCOPEN_FAILED );
						}
						else if ( ifo.shape.size() == 2 ) { //2d, matrix
							//uncomment these variables to use SZIP compression
							//unsigned szip_options_mask;
							//unsigned szip_pixels_per_block;
							int rank = 2;
							hsize_t dims[2] = { ifo.shape[0], ifo.shape[1] };
							hsize_t maxdims[2] = { ifo.shape[0], ifo.shape[1] };
							dspcid = H5Screate_simple( rank, dims, maxdims );
							if ( dspcid != H5I_INVALID_HID ) {
								plistid = H5Pcreate( H5P_DATASET_CREATE );
								if ( plistid != H5I_INVALID_HID ) {
									int chunk_rank = 2;
									hsize_t chunk_dims[2] = { ifo.chunk[0], ifo.chunk[1] };
									if ( H5Pset_chunk( plistid, chunk_rank, chunk_dims ) >= 0 ) {
										if ( ifo.compression == MYHDF5_COMPRESSION_GZIP ) {
											if ( H5Pset_deflate( plistid, (unsigned) ifo.compression_opts ) >= 0 ) {
												//set ZLIB / DEFLATE Compression using compression level 6.
												//szip_options_mask = H5_SZIP_NN_OPTION_MASK;
												//szip_pixels_per_block = 16;
												//status = H5Pset_szip (plistid, szip_options_mask, szip_pixels_per_block);
#ifdef MYHDF5_VERBOSE_LEVEL_TWO
												cout << "Setting deflate " << dsnm << " chunked, 2d" << "\n";
#endif
											}
											else {
												cerr << "Creating " << dsnm << " chunked, 2d set deflate failed !" << "\n";
												return nexus_close( MYHDF5_PLSTACCESS_FAILED );
											}
										}
										dsetid = H5Dcreate2( fileid, dsnm.c_str(), dtypid, dspcid, H5P_DEFAULT, plistid, H5P_DEFAULT );
										if ( dsetid != H5I_INVALID_HID ) {
											if ( H5Dwrite( dsetid, dtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, val.data() ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
												cout << "Writing " << dsnm << " chunked, 2d success" << "\n";
#endif

												//##MK::write attributes
												int attr_status = nexus_write_attributes( dsnm, attrs );
												if ( attr_status == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
													cout << "Writing " << dsnm << " contiguous, 2d attributes success !" << "\n";
#endif
												}
												else {
													cerr << "Writing " << dsnm << " contiguous, 2d attributes failed " << attr_status << "\n";
												}
												return nexus_close( MYHDF5_SUCCESS );
											}
											cerr << "Writing " << dsnm << " chunked, 2d failed !" << "\n";
											return nexus_close( MYHDF5_DSETWRITE_FAILED );
										}
										cerr << "Writing " << dsnm << " chunked, 2d failed  !" << "\n";
										return nexus_close( MYHDF5_DSETWRITE_FAILED );
									}
									cerr << "Creating " << dsnm << " chunked, 2d chunking failed !" << "\n";
									return nexus_close( MYHDF5_DSETWRITE_FAILED );
								}
								cerr << "Creating " << dsnm << " chunked, 2d plist failed !" << "\n";
								return nexus_close( MYHDF5_PLSTOPEN_FAILED );
							}
							cerr << "Creating " << dsnm << " chunked, 2d dspc failed !" << "\n";
							return nexus_close( MYHDF5_DSPCOPEN_FAILED );
						}
						//ifo.shape.size() > 2 not supported !
						cerr << "Creating " << dsnm << " failed because ifo.shape.size() > 2 not implemented !" << "\n";
						return nexus_close( MYHDF5_NOT_IMPLEMENTED );
					} //done with chunked
				//} //done newly created dataset
				//cerr << "Creating " << dsnm << " failed because dsnm already exists !" << "\n";
				//return nexus_close( MYHDF5_ERROR_OBJ_EXISTS );
			}
			cerr << "Opening file " << h5resultsfn << " failed !" << "\n";
			return nexus_close( MYHDF5_FOPEN_FAILED );
		}
		cerr << "H5Core error on input size " << h5resultsfn << " failed !" << "\n";
		return MYHDF5_ERROR_INPUT_SIZE;
	}

	cerr << "WARNING: value array val for dataset " << dsnm << " is empty!" << "\n";
	return MYHDF5_SUCCESS;
}

template int HdfFiveSeqHdl::nexus_write( const string dsnm, io_info const & ifo, vector<unsigned char> & val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, io_info const & ifo, vector<char> & val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, io_info const & ifo, vector<unsigned short> & val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, io_info const & ifo, vector<short> & val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, io_info const & ifo, vector<unsigned int> & val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, io_info const & ifo, vector<int> & val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, io_info const & ifo, vector<unsigned long> & val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, io_info const & ifo, vector<long> & val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, io_info const & ifo, vector<float> & val, ioAttributes const & attrs );
template int HdfFiveSeqHdl::nexus_write( const string dsnm, io_info const & ifo, vector<double> & val, ioAttributes const & attrs );

//tool-level
int HdfFiveSeqHdl::nexus_write(	const string dsnm, io_info const & ifo, vector<string> & val, ioAttributes const & attrs )
{
	if ( val.size() > 0 ) {
		unsigned char dimensionality = MYHDF5_LIST_SINGLE;
		if ( val.size() > 1 ) {
			dimensionality = MYHDF5_LIST_LONG;
		}
		unsigned char memory = MYHDF5_IS_VARIABLE;
		unsigned char terminator = MYHDF5_NULLTERM;
		unsigned char encoding = MYHDF5_UTF8;

		size_t max_characters = 0;
		for( auto it = val.begin(); it != val.end(); it++ ) {
			if ( it->length() < max_characters )
				continue;
			else
				max_characters = it->length();
		}
		max_characters++;
		//we use Cstyle arrays which is why we need space for the terminator

#ifdef MYHDF5_VERBOSE_LEVEL_TWO
		cout << "dimensionality " << (int) dimensionality << "\n";
		cout << "memory " << (int) memory << "\n";
		cout << "terminator " << (int) terminator << "\n";
		cout << "encoding " << (int) encoding << "\n";
		cout << "ndims " << val.size() << "\n";
		cout << "sdims " << max_characters << "\n";
#endif
		return nexus_write_string_dataset( dsnm, val, max_characters, dimensionality, memory, terminator, encoding, attrs );
	}
	//cerr << "WARNING: value array val for dataset " << dsnm << " is empty!" << "\n";
	return MYHDF5_SUCCESS;
}


int HdfFiveSeqHdl::nexus_write(	const string dsnm, const string val, ioAttributes const & attrs )
{
	unsigned char dimensionality = MYHDF5_SCALAR;
	unsigned char memory = MYHDF5_IS_VARIABLE;
	unsigned char terminator = MYHDF5_NULLTERM;
	unsigned char encoding = MYHDF5_UTF8;
	size_t max_characters = val.length();
	max_characters++;
	//we use Cstyle arrays which is why we need space for the terminator

#ifdef MYHDF5_VERBOSE_LEVEL_TWO
	cout << "dimensionality " << (int) dimensionality << "\n";
	cout << "memory " << (int) memory << "\n";
	cout << "terminator " << (int) terminator << "\n";
	cout << "encoding " << (int) encoding << "\n";
	cout << "ndims " << val.size() << "\n";
	cout << "sdims " << max_characters << "\n";
#endif
	vector<string> values;
	values.push_back(val);
	return nexus_write_string_dataset( dsnm, values, max_characters, dimensionality, memory, terminator, encoding, attrs );
}


int HdfFiveSeqHdl::nexus_write(	const string dsnm, vector<string> & val,
		const unsigned char dimensionality, const unsigned char memory,
		const unsigned char terminator, const unsigned char encoding, ioAttributes const & attrs )
{
	size_t max_characters = 0;
	for( auto it = val.begin(); it != val.end(); it++ ) {
		if ( it->length() < max_characters )
			continue;
		else
			max_characters = it->length();
	}
	max_characters++;
	//we use Cstyle arrays which is why we need space for the terminator

#ifdef MYHDF5_VERBOSE_LEVEL_TWO
	cout << "dimensionality " << (int) dimensionality << "\n";
	cout << "memory " << (int) memory << "\n";
	cout << "terminator " << (int) terminator << "\n";
	cout << "encoding " << (int) encoding << "\n";
	cout << "ndims " << val.size() << "\n";
	cout << "sdims " << max_characters << "\n";
#endif
	return nexus_write_string_dataset( dsnm, val, max_characters, dimensionality, memory, terminator, encoding, attrs );
}


//low-level
int HdfFiveSeqHdl::nexus_write_string_dataset( const string dsnm, vector<string> & val, const size_t max_characters,
		const unsigned char dimensionality, const unsigned char memory, const unsigned char terminator, const unsigned char encoding,
		ioAttributes const & attrs )
{
	herr_t status = -1;
	dtypid = H5Tcopy(H5T_C_S1);
	if ( dtypid != H5I_INVALID_HID ) {
		if ( memory == MYHDF5_IS_VARIABLE ) {
			status = H5Tset_size(dtypid, H5T_VARIABLE);
		}
		else { //MYHDF5_IS_FIXED
			status = H5Tset_size(dtypid, max_characters);
		}
		if ( status >= 0 ) {
			switch(terminator)
			{
				case MYHDF5_NULLTERM: { status = H5Tset_strpad(dtypid, H5T_STR_NULLTERM); break; }
				case MYHDF5_NULLPAD: { status = H5Tset_strpad(dtypid, H5T_STR_NULLPAD); break; }
				case MYHDF5_SPACEPAD: { status = H5Tset_strpad(dtypid, H5T_STR_SPACEPAD); break; }
				default: { status = -1; break; }
			}
			if ( status >= 0 ) {
				switch(encoding)
				{
					case MYHDF5_ASCII: { status = H5Tset_cset(dtypid, H5T_CSET_ASCII); break; }
					case MYHDF5_UTF8: { status = H5Tset_cset(dtypid, H5T_CSET_UTF8); break; }
					default: { status = -1; break; }
				}
				if ( status >= 0 ) {
					if ( nexus_open( H5F_ACC_RDWR ) != H5I_INVALID_HID ) {
						if ( dimensionality == MYHDF5_SCALAR ) {
							dspcid = H5Screate( H5S_SCALAR );
						}
						else { //dimensionality == MYHDF5_LIST_SINGLE or LIST_LONG
							int rank = 1; hsize_t dims[1] = { val.size() };
							dspcid = H5Screate_simple( rank, dims, NULL );
						}
						if ( dspcid != H5I_INVALID_HID ) {
							dsetid = H5Dcreate2( fileid, dsnm.c_str(), dtypid, dspcid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
							if ( dsetid != H5I_INVALID_HID ) {
								if ( dimensionality == MYHDF5_SCALAR ) {
									if ( memory == MYHDF5_IS_VARIABLE ) {
										if ( H5Dwrite( dsetid, dtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, val.data() ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
											cout << "Writing " << dsnm << " dataset, string, scalar, variable success" << "\n";
#endif
											if ( nexus_write_attributes( dsnm, attrs ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
												cout << "Writing " << dsnm << " attributes, string, scalar, variable success" << "\n";
#endif
												return nexus_close( MYHDF5_SUCCESS );
											}
											cerr << "Writing " << dsnm << " attributes, string, scalar, variable failed!" << "\n";
											return nexus_close ( MYHDF5_ATTRWCOLL_FAILED );
										}
										cerr << "Writing " << dsnm << " dataset, string, scalar, variable failed!" << "\n";
										return nexus_close( MYHDF5_DSETWRITE_FAILED );
									}
									else { //MYHDF5_IS_FIXED
										status = 0;
										char* wdata = new char[max_characters];
										if ( wdata != NULL ) {
											strcpy(wdata, val[0].c_str());
											if ( H5Dwrite( dsetid, dtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
												cout << "Writing " << dsnm << " dataset, string, scalar, fixed success" << "\n";
#endif
												if ( nexus_write_attributes( dsnm, attrs ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
													cout << "Writing " << dsnm << " dataset, string, scalar, fixed success" << "\n";
#endif
													delete[] wdata; wdata = NULL;
													return nexus_close( MYHDF5_SUCCESS );
												}
												cerr << "Writing " << dsnm << " dataset, string, scalar, fixed failed!" << "\n";
												delete[] wdata; wdata = NULL;
												return nexus_close( MYHDF5_ATTRWCOLL_FAILED );
											}
											cerr << "Writing " << dsnm << " dataset, string, scalar, fixed failed!" << "\n";
											delete[] wdata; wdata = NULL;
											return nexus_close( MYHDF5_DSETWRITE_FAILED );
										}
										cerr << "Writing " << dsnm << " dataset, string, scalar, fixed failed!" << "\n";
										delete[] wdata; wdata = NULL;
										return nexus_close( MYHDF5_DSETWRITE_FAILED );
									}
								}
								else if ( dimensionality == MYHDF5_LIST_SINGLE ) {
									if ( memory == MYHDF5_IS_VARIABLE ) {
										if ( H5Dwrite( dsetid, dtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, val.data() ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
											cout << "Writing " << dsnm << " dataset, string, list_single, variable success" << "\n";
#endif
											if ( nexus_write_attributes( dsnm, attrs ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
												cout << "Writing " << dsnm << " attributes, string, list_single, variable success" << "\n";
#endif
												return nexus_close( MYHDF5_SUCCESS );
											}
											cerr << "Writing " << dsnm << " attributes, string, list_single, variable failed!" << "\n";
											return nexus_close( MYHDF5_ATTRWCOLL_FAILED );
										}
										cerr << "Writing " << dsnm << " dataset, string, list_single, variable failed!" << "\n";
										return nexus_close( MYHDF5_DSETWRITE_FAILED );
									}
									else { //MYHDF5_IS_FIXED
										status = 0;
										char* wdata = new char[max_characters];
										if ( wdata != NULL ) {
											for( size_t i = 0; i < val.at(0).length(); i++ ) {
												wdata[i] = val[0][i];
											}
											for( size_t j = val.at(0).length(); j < max_characters; j++ ) {
												switch(terminator)
												{
													case MYHDF5_NULLTERM: { wdata[j] = '\0'; break; }
													case MYHDF5_NULLPAD: { wdata[j] = 0; break; }
													case MYHDF5_SPACEPAD: { wdata[j] = ' '; break; }
													default: { break; }
												}
											}
										}
										else {
											status = -1;
										}
										if ( status >= 0 ) {
											if ( H5Dwrite( dsetid, dtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
												cout << "Writing " << dsnm << " dataset, string, list_single, fixed success" << "\n";
#endif
												if ( nexus_write_attributes( dsnm, attrs ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
													cout << "Writing " << dsnm << " attibutes, string, list_single, fixed success" << "\n";
#endif
													delete[] wdata; wdata = NULL;
													return nexus_close( MYHDF5_SUCCESS );
												}
												cerr << "Writing " << dsnm << " attributes, string, list_single, fixed failed!" << "\n";
												delete[] wdata; wdata = NULL;
												return nexus_close( MYHDF5_ATTRWCOLL_FAILED );
											}
											cerr << "Writing " << dsnm << " dataset, string, list_single, fixed failed!" << "\n";
											delete[] wdata; wdata = NULL;
											return nexus_close( MYHDF5_DSETWRITE_FAILED );
										}
										cerr << "Writing " << dsnm << " dataset, string, list_single, fixed failed!" << "\n";
										delete[] wdata; wdata = NULL;
										return nexus_close( MYHDF5_DSETWRITE_FAILED );
									}
								}
								else { //MYHDF5_LIST_LONG
									if ( memory == MYHDF5_IS_VARIABLE ) {
										status = 0;
										char** wdata = new char*[val.size()];
										if ( wdata != NULL ) {
											for ( size_t i=0; i < val.size(); i++ ) {
												wdata[i] = new char[val[i].size()+1];
												if ( wdata[i] != NULL ) {
													strcpy(wdata[i], val[i].c_str());
												}
												else {
													status = -1;
												}
											}
										}
										else {
											status = -1;
										}
										if ( status >= 0 ) {
											if ( H5Dwrite( dsetid, dtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
												cout << "Writing " << dsnm << " dataset, string, list_long, variable success" << "\n";
#endif
												if ( nexus_write_attributes( dsnm, attrs ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
													cout << "Writing " << dsnm << " attributes, string, list_long, variable success" << "\n";
#endif
												    for ( size_t i=0; i < val.size(); i++) {
												    	delete[] wdata[i]; wdata[i] = NULL;
												    }
												    delete[] wdata; wdata = NULL;
												    return nexus_close( MYHDF5_SUCCESS );
												}
												cerr << "Writing " << dsnm << " attributes, string, list_long, variable failed!" << "\n";
												for ( size_t i=0; i < val.size(); i++) {
											    	delete[] wdata[i]; wdata[i] = NULL;
											    }
											    delete[] wdata; wdata = NULL;
												return nexus_close( MYHDF5_ATTRWCOLL_FAILED );
											}
											cerr << "Writing " << dsnm << " dataset, string, list_long, variable failed!" << "\n";
										    for ( size_t i=0; i < val.size(); i++) {
										    	delete[] wdata[i]; wdata[i] = NULL;
										    }
										    delete[] wdata; wdata = NULL;
											return nexus_close( MYHDF5_DSETWRITE_FAILED );
										}
										cerr << "Writing " << dsnm << " dataset, string, list_long, variable failed!" << "\n";
										for ( size_t i=0; i < val.size(); i++) {
										    delete[] wdata[i]; wdata[i] = NULL;
										}
										delete[] wdata; wdata = NULL;
										return nexus_close( MYHDF5_DSETWRITE_FAILED );
									}
									else { //MYHDF5_IS_FIXED
										status = 0;
										char* wdata = new char[val.size() * max_characters];
										if ( wdata != NULL ) {
											for( size_t i = 0; i < val.size(); i++ ) {
												for( size_t j = 0; j < val.at(i).length(); j++ ) {
													wdata[i*max_characters+j] = val[i][j];
												}
												for( size_t j = val.at(i).length(); j < max_characters; j++ ) {
													switch(terminator)
													{
														case MYHDF5_NULLTERM: { wdata[i*max_characters+j] = '\0'; break; }
														case MYHDF5_NULLPAD: { wdata[i*max_characters+j] = 0; break; }
														case MYHDF5_SPACEPAD: { wdata[i*max_characters+j] = ' '; break; }
														default: { break; }
													}
												}
											}
										}
										else { status = -1; }
										if ( status >= 0 ) {
											if ( H5Dwrite( dsetid, dtypid, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata ) >= 0 ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
												cout << "Writing " << dsnm << " dataset, string, list_long, fixed success" << "\n";
#endif
												if ( nexus_write_attributes( dsnm, attrs ) == MYHDF5_SUCCESS ) {
#ifdef MYHDF5_VERBOSE_LEVEL_ONE
													cout << "Writing " << dsnm << " attributes, string, list_long, fixed success" << "\n";
#endif
													delete[] wdata; wdata = NULL;
													return nexus_close( MYHDF5_SUCCESS );
												}
												cerr << "Writing " << dsnm << " attributes, string, list_long, fixed failed!" << "\n";
												delete[] wdata; wdata = NULL;
												return nexus_close( MYHDF5_ATTRWCOLL_FAILED );
											}
											cerr << "Writing " << dsnm << " dataset, string, list_long, fixed failed!" << "\n";
											delete[] wdata; wdata = NULL;
											return nexus_close( MYHDF5_DSETWRITE_FAILED );
										}
										cerr << "Writing " << dsnm << " dataset, string, list_long, fixed failed!" << "\n";
										delete[] wdata; wdata = NULL;
										return nexus_close( MYHDF5_DSETWRITE_FAILED );
									}
								}
								//should not be encountered
								cerr << "Writing " << dsnm << " failed unexpectedly !" << "\n";
								return nexus_close( MYHDF5_FAILED );
							}
							cerr << "Creating dset2 " << dsnm << " failed !" << "\n";
							return nexus_close( MYHDF5_DSETOPEN_FAILED );
						}
						cerr << "Creating dset1 " << dsnm << " failed !" << "\n";
						return nexus_close( MYHDF5_DSETCREATE_FAILED );
					}
					cerr << "Opening " << fileid << " failed !" << "\n";
					return nexus_close( MYHDF5_FOPEN_FAILED );
				} //encoding
			} //terminator
		} //setsize
		return nexus_close( MYHDF5_FAILED );
	}
	return nexus_close( MYHDF5_SUCCESS );
}
