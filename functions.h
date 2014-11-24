#ifndef _FUNCTIONS_H
#define _FUNCTIONS_H

#include <exception>
#include <sstream>
using namespace std;
#include <vector>
#include <stdio.h>
#include <inttypes.h> 

template < typename T >
unsigned long totalRAM( const vector<T>& vec ) {
  return sizeof( vec ) + sizeof( T ) * vec.capacity();
}
//unsigned long totalRAM( const string& vec ) {
//  return sizeof( vec ) + sizeof( vec ) * vec.capacity();
//}

/* random number generator */
double unifrand();
double ran_normal();

/* get system variable */
string get_env_var( std::string const & key );

/* check if a file exists */
bool file_exist(string filename) ;

/* check if a file i binary */
bool is_binary(string filename);

/* check information from a pid file, or any file 
 * equivalent to grep $fields file */
string procpidstatus(string file, string fields);
/* check information from a pid file, /proc/$pid/status 
 * equivalent to grep $fields  /proc/$pid/status */
string procpidstatus(int pid, string fields);

/* check kmer length 
 */
int kmerlength(string& SEQ) ;


/* get a line from FILE, normally a pipe */
bool fgetline(FILE *fp, string& read);

void longestCommonSubstring(const string& str1, const string& str2, 
			    int& maxSubstr, int& p1, int& p2);

int edit_distance( const std::string& s1, const std::string& s2 );

string to_lower(string word);
string to_upper(string word);

/* case insensitive equal find */
size_t ci_find(const string& str1, const string& str2);
bool ci_equal(const string& str1, const string& str2);

/* split string into array; append at the end */
void stringsplit(string& line, vector<int>& arr, char delimiter);
void stringsplit(string& line, vector<float>& arr, char delimiter);
void stringsplit(string& line, vector<double>& arr, char delimiter);
void stringsplit(string& line, vector<string>& arr, char delimiter);
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/


/* get last/first line of a file */
string lastline(string file);
string firstline(string file);

string commify(long n1);

void arrayindex(vector<double>& val, vector<int>& rk, int order);
void arrayindex(vector<float>& val, vector<int>& rk, int order);
void arrayindex(vector<int>& val, vector<int>& rk, int order);
void arrayindex(vector<uint64_t>& val, vector<int>& rk, int order);
void arrayindex(vector<int64_t>& val, vector<int>& rk, int order);
void arrayrank(vector<double>& val, vector<int>& rk, int order);
void arrayrank(vector<float>& val, vector<int>& rk, int order);
void arrayrank(vector<int>& val, vector<int>& rk, int order);
void arrayrank(vector<uint64_t>& val, vector<int>& rk, int order);
void arrayrank(vector<int64_t>& val, vector<int>& rk, int order);

void runmed(double* y, double* smo, int n, int band);
void runmed(float* y, float* smo, int n, int band);
void runmed(int* y, int* smo, int n, int band);
void runmed(unsigned short int* y, unsigned short int* smo, int n, int band);

// convert number or anything to string
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
template <class T>
inline std::string to_string (const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

bool check_github_update(string compiledtime, string updatefile);

#endif
