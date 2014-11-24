/*************************************************************************
BEGIN OF LICENSE
Copyright (c) Yinghua Wu and Hongzhe Li (rsicnv project).

This program is free software; you can redistribute and/or modify
the codes written by the authors under the terms of the GNU General 
Public License as published by the Free Software Foundation 
(www.fsf.org); either version 2 of the License, or (at your option) 
any later version. 

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
END OF LICENSE

CONTACT: wu_yinghua@hotmail.com; hongzhe@upenn.edu
*************************************************************************/
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <complex>
#include <string>
#include <cstring>
#include <limits>
#include <algorithm>
using namespace std;

#include "functions.h"  // functions defined here


/* Simple normal random number generator, copied from genran.c;
   Simple normal random number generator, copied from wgsim.c */
double unifrand() { return rand() / double(RAND_MAX); }
double ran_normal()
{ 
  static int iset = 0; 
  static double gset; 
  double fac, rsq, v1, v2; 
  if (iset == 0) {
    do { 
      v1 = 2.0 * unifrand() - 1.0;
      v2 = 2.0 * unifrand() - 1.0; 
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0 * log(rsq) / rsq); 
    gset = v1 * fac; 
    iset = 1;
    return v2 * fac;
  } else {
    iset = 0;
    return gset;
  }
}

/* Get system variable, copied from
   http://stackoverflow.com/questions/631664/accessing-environment-variables-in-c */
string get_env_var( std::string const & key ) {                                 
  char * val;                                                                        
  val = getenv( key.c_str() );                                                       
  std::string retval = "";                                                           
  if (val != NULL) retval = val;
  return retval;                                                                        
}

bool file_exist(string filename) 
{
  FILE* fp = fopen(filename.c_str(), "r");
  if (fp) { fclose(fp); return(true); }
  else  return(false); 
}

bool is_binary(string filename)
{
  if ( filename=="-" || filename=="STDIN" || filename=="STDOUT" ) return false;
  
  FILE * file;
  file = fopen(filename.c_str(), "rb");
  unsigned char magicNumber[4]={0};
  size_t i= fread(magicNumber, 2, 1, file);
  fclose(file);                    
  //cerr << (int)magicNumber[0] << "\t" << (int)magicNumber[1] << endl;  
  if ( magicNumber[0]==31 && magicNumber[1]==139 ) return true;
  if ( ( magicNumber[0]<=31 || magicNumber[0]>127 ) && 
       ( magicNumber[1]<=31 || magicNumber[1]>127 ) ) return true;
  return false;
}

string procpidstatus(string file, string fields)
{
  string value="";
  ifstream inp(file.c_str());
  if ( !inp ) {
    cerr << "file " << file << " does not exist" << endl;
    return value;
  }
  
  while ( !inp.eof() ) {
    string tmps;
    getline(inp,tmps);
    if ( ci_find(tmps, fields ) != string::npos ) value+=tmps+"\n";
  }
  inp.close();
  return(value);
}
string procpidstatus(int pid, string fields)
{
  string value="";
  
  string file="/proc/"+to_string(pid)+"/status";
  ifstream inp(file.c_str());
  if ( !inp ) {
    cerr << "file " << file << " does not exist" << endl;
    return value;
  }
  
  while ( !inp.eof() ) {
    string tmps;
    getline(inp,tmps);
    if ( ci_find(tmps, fields ) != string::npos ) value+="#"+tmps+"\n";
  }
  inp.close();
  return(value);
}


int kmerlength(string& SEQ) 
{
  size_t i,k;
  size_t maxRepeat=5;
  int maxLength=0;
  string rep="";
  string maxrep="";
  for(i=0;i<SEQ.size();++i) {
    for(k=1;k<=maxRepeat;++k) {
      if ( i+k>=SEQ.size() ) continue;
      rep=SEQ.substr(i,k);
      
      int ndiff=0;
      int lRepeat=0;
      int idx=i+k;
      while(ndiff<=1 && idx+k<SEQ.size() ) {
	for(size_t i1=idx, i2=0; i1<idx+k; ++i1,++i2)
	  if ( rep[i2]!=SEQ[i1] ) ++ndiff;
	if ( ndiff<=1 ) lRepeat+=k;
	idx+=k;
      }
      
      if ( lRepeat > maxLength ) {
	maxLength=lRepeat;
	maxrep=rep;
      }
    }
  }
  return maxLength;

}


bool fgetline(FILE *fp, string& read)
{
  const int SIZEBUF=1024;
  char str[SIZEBUF];
  read="";
  while( fgets(str, sizeof (str), fp) ) {
    int len=strlen(str);
    read += str;
    if ( len <= (SIZEBUF-1) && str[len-1] == '\n') break;
    //cerr << len << "  keep reading " << endl;
  }
  if (!read.empty() && 
      read[read.length()-1] == '\n') read.erase(read.length()-1);  
  
  if ( !feof(fp) ) return true;
  else return false;
}

int edit_distance( const std::string& s1, const std::string& s2 )
{
  const int cost_del = 1;
  const int cost_ins = 1;
  const int cost_sub = 1;
  
  int c1[4096];
  int c2[4096];
  int* p = c1;
  int* q = c2;
  int* r;
  
  int i,j;
  int n1 = s1.length();
  int n2 = s2.length();
  
  if ( n1 > 4096 || n2 > 4096 ) {
    j=std::max(n1,n2)+1;
    int c11[j];
    int c22[j];
    p=c11;
    q=c22;
    cerr << "check edit_distance() in function.cpp, change size of new and recomiple\n";
    exit(0);
  }
  
  p[0] = 0;
  for( j = 1; j <= n2; ++j ) p[j] = p[j-1] + cost_ins;
  
  for( i = 1; i <= n1; ++i ){
    q[0] = p[0] + cost_del;
    for( j = 1; j <= n2; ++j ) {
      int d_del = p[j] + cost_del;
      int d_ins = q[j-1] + cost_ins;
      int d_sub = p[j-1] + ( s1[i-1] == s2[j-1] ? 0 :cost_sub );
      q[j] = std::min( std::min( d_del, d_ins ), d_sub );
    }
    r = p;
    p = q;
    q = r;
  }

  return p[n2];
  
}



void longestCommonSubstring(const string& str1, const string& str2, 
				int& maxSubstr, int& p1, int& p2)
{
  maxSubstr=p1=p2=-1;
  
  int c1[1024];
  int c2[1024];
  int *curr=c1;
  int *prev=c2;
  int *swap = 0;
  
  if ( str1.empty() || str2.empty() ) return;
  if ( str1.size() > 1024 || str2.size() > 1024 ) {
    int c11[str1.size()+str2.size()];
    int c22[str1.size()+str2.size()];
    curr=c11;
    prev=c22;
    cerr << "check file matchreads, change size of new and recomiple\n";
    exit(0);
  }
  
  for(int i = 0; i<(int)str1.size(); ++i)  {
    for(int j = 0; j<(int)str2.size(); ++j)  {
      
      if ( str1[i] != str2[j] ) curr[j] = 0;
      else {
	if ( i == 0 || j == 0) curr[j] = 1;
	else { curr[j] = 1 + prev[j-1]; }
	if ( maxSubstr < curr[j] ) {
	  maxSubstr = curr[j]; 
	  p1=i-maxSubstr+1;
	  p2=j-maxSubstr+1;
	}
      }
    }
    
    swap=curr;
    curr=prev;
    prev=swap;
  }
  
  return ;
}

/* split string into array; append at the end */
void stringsplit(string& line, vector<string>& arr, char delimiter)
{
  arr.clear();
  istringstream iss(line);
  string tmps;
  while( getline(iss, tmps, delimiter ) ) arr.push_back(tmps);
  return;
}
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
template<class Etype>
void stringsplit_tp(string& line, vector<Etype>& arr, char delimiter)
{
  arr.clear();
  istringstream iss(line);
  string tmps;
  while( getline(iss, tmps, delimiter ) ) {
    //    cerr << tmps << "|\t";
    istringstream iss2(tmps);
    Etype tmp;
    iss2 >> tmp;
    arr.push_back(tmp);
  }
  return;
}
void stringsplit(string& line, vector<int>& arr, char delimiter)
{ stringsplit_tp(line, arr, delimiter);  }
void stringsplit(string& line, vector<float>& arr, char delimiter)
{ stringsplit_tp(line, arr, delimiter);  }
void stringsplit(string& line, vector<double>& arr, char delimiter)
{ stringsplit_tp(line, arr, delimiter);  }
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/


string to_lower(string word) {
  std::transform(word.begin(), word.end(), word.begin(), ::tolower);
  return word;
}
string to_upper(string word) {
  std::transform(word.begin(), word.end(), word.begin(), ::toupper);
  return word;
}

/* case insensitive find 
 * find the position in first string that match the second
 */
bool ci_equal_char(char ch1, char ch2)
{ return toupper((unsigned char)ch1) == toupper((unsigned char)ch2); }
size_t ci_find(const string& str1, const string& str2)
{
  string::const_iterator pos = 
    std::search(str1.begin(), str1.end(), 
		str2.begin(), str2.end(), 
		ci_equal_char);
  if (pos == str1.end() )  return string::npos;
  else return pos-str1.begin();
}
/* case insensitive find */

/* case insensitive equal
 */
bool ci_equal(const string& str1, const string& str2)
{
  if ( str1.length() != str2.length() ) return(false); 
  for(size_t i=0;i<str1.length();++i)
    if ( ! ci_equal_char(str1[i], str2[i]) ) return(false); 
  return(true);
}
/* case insensitive equal */



/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
string lastline(string file) {
  ifstream inp(file.c_str());
  if ( !inp ) {
    cout << "File " << file << " not exist" << endl; 
    exit(0);
  }
  //cout << file << endl;

  inp.seekg(0, ios::end);
  long int current;
  current=inp.tellg();
  if ( current < 0 ) {
    cout << "file end position " << current << endl;
    cout << "change type of current to longer bytes to accomodate size of the file" << endl;
    cout << "change type of indices in other places too!" << endl;
    exit(0);
  }

  /*
    cout << current << endl;
    cout << std::numeric_limits<int>::min() << endl;
    cout << std::numeric_limits<int>::max() << endl;
    cout << std::numeric_limits<long int>::min() << endl;
    cout << std::numeric_limits<long int>::max() << endl;
    exit(0);
  */

  string lastLine;
  long int i=100;
  i= current>i? i:current-1;
  inp.seekg(-i, ios::end);
  char b;
  inp.read(&b,1);
  int foundnr=0;
  while ( b != '\n' || foundnr<2 ){ 
    // inp.seekg(length-i,ios::beg);
    i++;
    if (b != '\n') foundnr++;
    inp.seekg(-i,ios::end);
    inp.read(&b,1);
    //cout << b ;
  }
  while( !inp.eof() ) {
    string tmps1;
    getline(inp,tmps1);
    // cout << "---[" << tmps1 << "]---" << endl;
    if ( tmps1.length()>2 ) lastLine=tmps1;
  }
  inp.close();

  return(lastLine);
}  

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
string firstline(string file) {
  string firstLine;
  
  ifstream inp(file.c_str());
  if ( !inp ) {
    cout << "File " << file << " not exist" << endl; 
    exit(0);
  }
  
  while ( !inp.eof() ) {
    string tmps;
    getline(inp,tmps);
    if (tmps.length() < 2) continue;
    if (tmps[0] == '#' ) continue;
    firstLine=tmps;
    break; 
  };
  inp.close();
  
  return(firstLine);
}  

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
string commify(long n1)
{
  string s="";
  unsigned long n;
  n = n1>0 ? n1 : -n1;
  
  int cnt = 0;
  do
  {
    s.insert(0, 1, char('0' + n % 10));
    n /= 10;
    if (++cnt == 3 && n)
    {
      s.insert(0, 1, ',');
      cnt = 0;
    }
  } while (n);
  
  if ( n1<0 ) s="-"+s;
  
  return s;
}

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
template<class Etype>
void arrayindex_tp(vector<Etype>& val, vector<int>& rk, int order)
{
  
  if ( val.size() != rk.size() ) {
    cout << "array sizes not same" << endl
	 << "val.size()=" << val.size() << endl
	 << "rk.size()=" << rk.size() << endl
	 << "arrayindex_tp()" << endl;
    exit(0);
  }
  if ( order==0 ) {
    cout << "int order>0 : increasing order from idx 0~n" << endl
	 << "int order<0 : decreasing order from idx 0~n" << endl
	 << "int order=0 : no such option" << endl
	 << "order=" << order << endl;
    exit(0);
  }

  order = order>0 ? 1:-1;  

  size_t i;
  size_t itmp;
  for(i=0;i<rk.size();++i) rk[i]=i;
  
  bool swapped = false;
  long int k=0;
  do {
    swapped = false;
    if ( k<0 ) k=0;
    for(i=k; i<val.size()-1; ++i) {
      if ( order>0 && val[ rk[i] ] > val[ rk[i+1] ] ) swapped = true;
      if ( order<0 && val[ rk[i] ] < val[ rk[i+1] ] ) swapped = true;
      if ( swapped ) { itmp=rk[i]; rk[i]=rk[i+1]; rk[i+1]=itmp; k=i-order; break;}
    }
  } while ( swapped );
  
} //void arrayindex_tp()
void arrayindex(vector<double>& val, vector<int>& rk, int order)
{  arrayindex_tp(val,rk,order); }
void arrayindex(vector<float>& val, vector<int>& rk, int order)
{  arrayindex_tp(val,rk,order); }
void arrayindex(vector<int>& val, vector<int>& rk, int order)
{  arrayindex_tp(val,rk,order); }
void arrayindex(vector<uint64_t>& val, vector<int>& rk, int order)
{  arrayindex_tp(val,rk,order); }
void arrayindex(vector<int64_t>& val, vector<int>& rk, int order)
{  arrayindex_tp(val,rk,order); }

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
template<class Etype>
void arrayrank_tp(vector<Etype>& val, vector<int>& rk, int order)
{
  cout << "code not checked yet!!!!!!!!!!!!!!!" << endl;  
  cout << "code not checked yet!!!!!!!!!!!!!!!" << endl;
  cout << "code not checked yet!!!!!!!!!!!!!!!" << endl; exit(0);  
  if ( val.size() != rk.size() ) {
    cout << "array sizes not same" << endl
	 << "val.size()=" << val.size() << endl
	 << "rk.size()=" << rk.size() << endl
	 << "arrayrank_tp()" << endl;
    exit(0);
  }
  if ( order==0 ) {
    cout << "int order>0 : increasing order from idx 0~n" << endl
	 << "int order<0 : decreasing order from idx 0~n" << endl
	 << "int order=0 : no such option" << endl
	 << "order=" << order << endl;
    exit(0);
  }
  
  size_t i;
  int itmp;
  for(i=0;i<rk.size();++i) rk[i]=i;
  
  bool swapped = false;
  do {
    swapped = false;
    for(i=0;i<val.size()-1;++i) {
      if ( order>0 && val[ rk[i] ] > val[ rk[i+1] ] ) swapped = true;
      if ( order<0 && val[ rk[i] ] < val[ rk[i+1] ] ) swapped = true;
      if ( swapped ) {itmp=rk[i]; rk[i]=rk[i+1]; rk[i+1]=itmp; }
    }
  } while ( swapped );
  
  vector<int> idx(rk.size());
  for(i=0;i<rk.size();++i) idx[ rk[i] ] = i ;
  rk=idx;  
} //void arrayrank_tp()
void arrayrank(vector<double>& val, vector<int>& rk, int order)
{  arrayrank_tp(val,rk,order); }
void arrayrank(vector<float>& val, vector<int>& rk, int order)
{  arrayrank_tp(val,rk,order); }
void arrayrank(vector<int>& val, vector<int>& rk, int order)
{  arrayrank_tp(val,rk,order); }
void arrayrank(vector<uint64_t>& val, vector<int>& rk, int order)
{  arrayrank_tp(val,rk,order); }
void arrayrank(vector<int64_t>& val, vector<int>& rk, int order)
{  arrayrank_tp(val,rk,order); }
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/


/* begin runmed 
 *
 **/
template<class Etype>
void Srunmed(Etype* y, Etype* smo, int n, int band,
	     int end_rule, int debug)
{
  Etype rmed, rmin, temp, rnew, yout, yi;
  Etype rbe, rtb, rse, yin, rts;
  int imin, ismo, i, j, first, last, band2, kminus, kplus;
  int bw;

  if ( band % 2 != 1 ) {
    cout << "band must be odd" << endl;
    cout << "band =" << band << endl;
    cout << "runmed" << endl;
    exit(0);
  }
  if( band > n) {
    cout << "band = " << band << endl
	 << "n = " << n 
	 << "bandwidth/span of running medians is larger than n" << endl
	 << endl;
    exit(0);
  }
  if ( n<= 0 ) {
    cerr << "Srunmed(): n<=0" << endl;
    exit(0);
  }

  bw=band;

  //Array<Etype> scrat(bw);
  vector<Etype> scrat(bw);
  

  /* 1. Compute  'rmed' := Median of the first 'band' values
     ======================================================== */
  
  for (i = 0; i < bw; ++i ) scrat[i] = y[i];
  
  /* find minimal value  rmin = scrat[imin] <= scrat[j] */
  rmin = scrat[0]; imin = 0;
  for (i = 1; i < bw; ++i )
    if (scrat[i] < rmin) { rmin = scrat[i]; imin = i; }

  /* swap scrat[0] <-> scrat[imin] */
  temp = scrat[0]; scrat[0] = rmin; scrat[imin] = temp;
  
  /* sort the rest of of scrat[] by bubble (?) sort -- */
  for (i = 2; i < bw; ++i) {
    if (scrat[i] < scrat[i - 1]) {/* find the proper place for scrat[i] */
      temp = scrat[i];
      j = i;
      do {
	scrat[j] = scrat[j - 1];
	--j;
      } while (scrat[j - 1] > temp); /* now:  scrat[j-1] <= temp */
      scrat[j] = temp;
    }
  }
  
  band2 = bw / 2;
  rmed = scrat[band2]; /* == Median( y[(1:band2)-1] ) */
  /* "malloc" had  free( (char*) scrat);*/ /*-- release scratch memory --*/
  
  if( end_rule == 0) { /*-- keep DATA at end values */
    for (i = 0; i < band2; ++i) smo[i] = y[i];
  }
  else { /* if(*end_rule == 1)  copy median to CONSTANT end values */
    for (i = 0; i < band2; ++i) smo[i] = rmed;
  }
  smo[band2] = rmed;
  band2++; /* = bw / 2 + 1*/;
  
  if( debug ) cout << bw << "\t" << band2 << endl;
  
  /* Big	FOR Loop: RUNNING median, update the median 'rmed'
     ======================================================= */
  for (first = 1, last = bw, ismo = band2;
       last < n;
       ++first, ++last, ++ismo) {
    
    yin = y[last];
    yout = y[first - 1];
    
    if( debug ) 
      cout << "is=" << ismo << " y(in/out)= " << yin << "\t" << yout << endl;
    
    rnew = rmed; /* New median = old one   in all the simple cases --*/
	
    if (yin < rmed) {
      if (yout >= rmed) {
	kminus = 0;
	if (yout > rmed) {/*	--- yin < rmed < yout --- */
	  if( debug) cout << ": yin < rmed < yout " << endl;
	  rnew = yin;/* was -rinf */
	  for (i = first; i <= last; ++i) {
	    yi = y[i];
	    if (yi < rmed) {
	      ++kminus;
	      if (yi > rnew) rnew = yi;
	    }
	  }
	  if (kminus < band2)		rnew = rmed;
	}
	else {/*		--- yin < rmed = yout --- */
	  if( debug) cout << ": yin < rmed == yout " << endl;
	  rse = rts = yin;   /* was -rinf */
	  for (i = first; i <= last; ++i) {
	    yi = y[i];
	    if (yi <= rmed) {
	      if (yi < rmed) {
		++kminus;
		if (yi > rts)	rts = yi;
		if (yi > rse)	rse = yi;
	      } else		rse = yi;
	      
	    }
	  }
	  rnew = (kminus == band2) ? rts : rse ;
	  if(debug) cout << "k- : " << kminus << endl;
	}
      } /* else: both  yin, yout < rmed -- nothing to do .... */
    }
    else if (yin != rmed) { /* yin > rmed */
      if (yout <= rmed) {
	kplus = 0;
	if (yout < rmed) {/* -- yout < rmed < yin --- */
	  if(debug) cout << ": yout < rmed < yin " << endl;
	  rnew = yin; /* was rinf */
	  for (i = first; i <= last; ++i) {
	    yi = y[i];
	    if (yi > rmed) {
	      ++kplus;
	      if (yi < rnew)	rnew = yi;
	    }
	  }
	  if (kplus < band2)	rnew = rmed;
	  
	} else { /* -- yout = rmed < yin --- */
	  if(debug) cout << ": yout == rmed < yin " << endl;
	  rbe = rtb = yin; /* was rinf */
	  for (i = first; i <= last; ++i) {
	    yi = y[i];
	    if (yi >= rmed) {
	      if (yi > rmed) {
		++kplus;
		if (yi < rtb)	rtb = yi;
		if (yi < rbe)	rbe = yi;
	      } else		rbe = yi;
	    }
	  }
	  rnew = (kplus == band2) ? rtb : rbe;
	  if(debug) cout << "k+ : " << kplus << endl;
	}
      } /* else: both  yin, yout > rmed --> nothing to do */
    } /* else: yin == rmed -- nothing to do .... */
    if(debug) cout << "=> " << rmed << "\t" <<  rnew << endl;
    rmed = rnew;
    smo[ismo] = rmed;
  } /*     END FOR ------------ big Loop -------------------- */
  
  if(end_rule == 0) { /*-- keep DATA at end values */
    for (i = ismo; i < n; ++i)  smo[i] = y[i];
  }
  else { /* if(*end_rule == 1)  copy median to CONSTANT end values */
    for (i = ismo; i < n; ++i) smo[i] = rmed;
  }
} /* Srunmed */

void runmed(double* y, double* smo, int n, int band)
{ Srunmed(y, smo, n, band, 0, 0); }
void runmed(float* y, float* smo, int n, int band)
{ Srunmed(y, smo, n, band, 0, 0); }
void runmed(int* y, int* smo, int n, int band)
{ Srunmed(y, smo, n, band, 0, 0); }
void runmed(unsigned short int* y, unsigned short int* smo, int n, int band)
{ Srunmed(y, smo, n, band, 0, 0); }


bool check_github_update(string compiledTime, string gitupdate)
{
  string updatedTime="NULL";
  time_t compiled=0, gittime=0;
  struct tm tm1;
  
  if ( strptime(compiledTime.c_str(), "%a %b %d %H:%M:%S %Y", &tm1)!=NULL )
    compiled=mktime(&tm1); 
  
  /*
  cerr << tm1.tm_year << " "
       << tm1.tm_mon << " "
       << tm1.tm_mday << " "
       << tm1.tm_hour << " "
       << tm1.tm_min << " "
       << tm1.tm_sec << endl;
  */
  
  string cmd="curl --max-time 3 " + gitupdate + " 2>/dev/null";
  FILE* pipe = popen(cmd.c_str(), "r");
  if (pipe) {
    fgetline(pipe, updatedTime);
    if ( strptime(updatedTime.c_str(), "%a %b %d %H:%M:%S %Y", &tm1)!=NULL )
      gittime=mktime(&tm1);
  }
  pclose(pipe);
  
  if ( gittime==0 ) {
    cmd="wget --timeout=3 -O - " + gitupdate + " 2>/dev/null";
    pipe = popen(cmd.c_str(), "r");
    if (pipe) {
      fgetline(pipe, updatedTime);
      if ( strptime(updatedTime.c_str(), "%a %b %d %H:%M:%S %Y", &tm1)!=NULL )
	gittime=mktime(&tm1); 
    }
    pclose(pipe);
  }
  
  /*
  cerr << tm1.tm_year << " "
       << tm1.tm_mon << " "
       << tm1.tm_mday << " "
       << tm1.tm_hour << " "
       << tm1.tm_min << " "
       << tm1.tm_sec << endl;
  */
  
  
  //cerr << compiledTime << "\t" << compiled << "\n"
  //    << updatedTime << "\t" << gittime << endl
  //    << difftime(gittime, compiled) << endl;
  
  if ( gittime==0 ) gittime=compiled;
  
  if ( compiled != gittime ) {
    if ( compiled < gittime ) 
      cerr << "Update:\n"
	   << "  matchclips is updated at github on " << updatedTime << " UTC\n"
	   << "  your version was stamped " << compiledTime << " UTC\n"
	   << "  please download from https://github.com/yhwu/matchclips2\n"
	   << endl; 
    return true;
  }
  
  return false;
}

