#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <bitset>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>
#include <sys/time.h>
#include <cstdint>
#include <stdint.h>
#include <unistd.h>


// Includes for seeding
#include <boost/unordered_map.hpp>

using namespace std;

// Set the word/hash sizes here at compile time, must be able to cast from a ulong
float PAR_MIN_MATCH = 1;
bool PAR_DIAGNOSTICS = 0;
int PAR_GAP = 1;
int MAX_seeds = 0;
int GLOBAL_READ_WORDS = 0;
int GLOBAL_CURRENT_WORD = 0;
int GLOBAL_SKIPPED_WORDS = 0;
bool PAR_BIN_OUT = 0;
bool PAR_HAPLOID = 0;

ofstream FOUT;

const int CONST_READ_AHEAD = 10;
const int WORD_SIZE = 64;
typedef uint64_t hash_size;

struct Marker {
	string id;
	long pos;
	double cm;
	
	string print() {
		stringstream ss;
		ss << id << '\t' << pos << '\t' << cm << endl;
		return ss.str();
	}
};

class Individual {
	string id[2];
	int idnum;
public:
	bitset<WORD_SIZE> hap[ CONST_READ_AHEAD ];

	void clear( int w ) {
		hap[ w % CONST_READ_AHEAD ].reset();
	}
	void setMarker( int w , int bit ) {
		hap[ w % CONST_READ_AHEAD ].set( bit );
	}
	hash_size getWordHash( int w ) {
		return hap[ w % CONST_READ_AHEAD ].to_ulong();
	}
	string getWordString( int w ) {
		return hap[ w % CONST_READ_AHEAD ].to_string();
	}
	string print() { 
		stringstream ss;
		ss << id[0] << '\t' << id[1] << '\t' << endl;
		return ss.str();	
	}
	string getID() { return id[1]; }
	int getNum() { return idnum; }
	Individual(string,string,int);
};
Individual::Individual( string aid , string bid , int iid ) {
	id[0] = aid;
	id[1] = bid;
	idnum = iid;
	for ( int w = 0; w < CONST_READ_AHEAD ; w++ ) clear( w );
}

vector< Individual > all_ind;
vector< Marker > all_markers;
unsigned int num_ind;

// Convenience function to compute genetic distance between two words (start of w1 and end of w2)
double cmBetween( int w1 , int w2 ) {
	int end = WORD_SIZE * w2 + WORD_SIZE - 1;
	if ( end >= all_markers.size() ) end = all_markers.size() - 1;
	return all_markers[end].cm - all_markers[WORD_SIZE * w1].cm;
}
long getPos( int w ) {
	if ( w >= all_markers.size() ) w = all_markers.size() - 1;
	return all_markers[w].pos;
}

struct Match {
	int interval[2] = {0,0};
	// Print this match
	// pair : identifiers for the corresponding individuals in all_ind
	void print( pair<unsigned int,unsigned int> p ) {
		double mlen = cmBetween(interval[0],interval[1]);
		if ( mlen >= PAR_MIN_MATCH ) {
			if ( PAR_BIN_OUT ) {
                unsigned int pid[2];
                pid[0] = all_ind[p.first].getNum();
                pid[1] = all_ind[p.second].getNum();
                unsigned int sid[2];
                sid[0] = interval[0] * WORD_SIZE;
                sid[1] = interval[1] * WORD_SIZE + WORD_SIZE - 1;
                FOUT.write( (char*) &pid[0] , sizeof( unsigned int ) );
                FOUT.write( (char*) &pid[1] , sizeof( unsigned int ) );
                FOUT.write( (char*) &sid[0] , sizeof( unsigned int ) );
                FOUT.write( (char*) &sid[1] , sizeof( unsigned int ) );			
			} else {
				FOUT 	<< all_ind[p.first].getID() << "\t" 
						<< all_ind[p.second].getID() << "\t" 
						<< getPos( interval[0] * WORD_SIZE ) << "\t" 
						<< getPos( interval[1] * WORD_SIZE + WORD_SIZE - 1 ) << "\t" 
						<< mlen << "\t" 
						<< interval[1] - interval[0] + 1 << "\t";
				FOUT << endl;
			}
		}
	}

	void extend( int w ) {
		if ( interval[1] < w ) interval[1] = w;
	}
	Match(int);
	Match(void);
};
Match::Match( int i ) { interval[0] = interval[1] = i; }
Match::Match() { interval[0] = interval[1] = 0; }

/* Object for storing extension between pairs of individuals */
class ExtendHash {
	boost::unordered_map<unsigned int, Match > extend_hash;
	unsigned int num = 0;
	// Empty Match to insert into hash
	Match m;
	// Iterator for testing insertion
	std::pair<boost::unordered::iterator_detail::iterator<boost::unordered::detail::ptr_node<std::pair<const unsigned int, Match > > >, bool> extend_ret;
public:
	ExtendHash(unsigned int);
	// Compute pair of individuals from location indicator
	pair<unsigned int,unsigned int> locationToPair( unsigned int loc ) {
		pair<unsigned int,unsigned int> p;
		// round everyone down to the nearest haplotype
		if ( !PAR_HAPLOID ) {
			p.second = 2 * (loc % num);
			p.first = 2 * ((loc - p.second/2) / num);
		} else {
			p.second = loc % num;
			p.first = (loc - p.second) / num;		
		}
		return p;
	}
	// Compute location from pair of individuals
	unsigned int pairToLocation( unsigned int i , unsigned int j ) {
	
		if ( !PAR_HAPLOID ) {
			// round everyone down to the nearest haplotype
			i = (i - (i % 2)) / 2;
			j = (j - (j % 2)) / 2;	
		}
		unsigned int loc = (i > j) ? j * num + i : i * num + j;
		return loc;
	}
	// Extend or add a given pair in the current hash
	// unsigned int i,j : identifiers for the two individuals
	// int w : current word # to extend or add
	void extendPair( unsigned int i , unsigned int j , int w ) {
		m.interval[0] = GLOBAL_CURRENT_WORD;
		// Find/extend this location in the hash
		extend_ret = extend_hash.insert( pair< unsigned int , Match >( pairToLocation(i,j) , m) );
		(extend_ret.first->second).extend(w);
	}
	
	// Remove all pairs that were not extended beyond w
	// int w : word # to remove prior to
	void clearPairsPriorTo(int w) {
		for ( auto it = extend_hash.begin() ; it != extend_hash.end() ; ) {
			if ( it->second.interval[1] < w ) {
				it->second.print( locationToPair(it->first) );
				it = extend_hash.erase( it );
			} else it++;
		}
	}

	// Remove all pairs that were not extended beyond w
	// int w : word # to remove prior to
	void extendAllPairsTo(int w) {
		for ( auto it = extend_hash.begin() ; it != extend_hash.end() ; it++ ) it->second.interval[1] = w;
	}
	
	// Remove all pairs
	// int w : word # to remove prior to
	void clearAllPairs() {
		for ( auto it = extend_hash.begin() ; it != extend_hash.end() ; ) {
			it->second.print( locationToPair(it->first) );
			it = extend_hash.erase( it );
		}
	}
	int size() {
		return extend_hash.size();
	}
};
ExtendHash::ExtendHash(unsigned int n) { num = n; }

/* Object for storing initial word seeds */
class SeedHash {
	boost::unordered_map<hash_size, vector<unsigned int> > seed_hash;
	// Empty vector to insert into the seed hash
	vector<unsigned int> vec;
	// Iterator for testing insertion of elements
	// std::pair<boost::unordered::iterator_detail::iterator<boost::unordered::detail::ptr_node<std::pair<const unsigned int, vector<unsigned int> > > >, bool> seed_ret;
public:
	void insertIndividual( unsigned int i , hash_size word ) {
		auto seed_ret = seed_hash.insert( pair<hash_size, vector<unsigned int>> ( word , vec ) );
		(seed_ret.first->second).push_back( i );
	}
	void clear() {
		seed_hash.clear();
	}
	int size() { 
		return seed_hash.size();
	}
	
	// Generate a new hash for this vector of individuals
	unsigned long subHash( ExtendHash * e , vector<unsigned int> vec , int w ) {
		SeedHash cur_sh;
		// seed the next word from this subset of individuals
		for ( int i = 0 ; i < vec.size() ; i++ ) cur_sh.insertIndividual( vec[i] , all_ind[ vec[i] ].getWordHash( w ) );
		// recursion:
		// cerr << "\tsubhash seeds: " << w << " " << cur_sh.size() << endl;
		return cur_sh.extendAllPairs( e , w );
	}
	// Extend/save all pairs in the current hash
	// ExtendHash * e : Pointer to ExtendHash which will be called for each pair
	// returns : number of pairs evaluated
	unsigned long extendAllPairs( ExtendHash * e , int w ) {
		unsigned long tot_pairs = 0;
		for ( auto it = seed_hash.begin() ; it != seed_hash.end() ; ++it ) {
		
			// *** As long as the # of pairs is high, generate a sub-hash for the next word
			// *** Only store pairs of individuals that have collision in a small hash
			// *** Extend only to the haplotypes that seeded here
			if ( MAX_seeds != 0 && it->second.size() > MAX_seeds && w + 1 < GLOBAL_READ_WORDS ) {
				// recursively generate a sub-hash
				// IMPORTANT: if we run out of buffered words then this seed does not get analyzed
				if ( w + 1 < GLOBAL_READ_WORDS ) tot_pairs += subHash( e , it->second , w + 1 );
				else GLOBAL_SKIPPED_WORDS++;
			} else {
				tot_pairs += it->second.size() * (it->second.size() - 1) / 2;
				for ( int i = 0 ; i < it->second.size() ; i++ ) {
				for ( int ii = i+1 ; ii < it->second.size() ; ii++ ) {
					e->extendPair( it->second[i] , it->second[ii] , w );
				}
				}	
			}
		}
		return tot_pairs;	
	}
	// Debug print all pairs in the current hash
	// string word : formatted output for end of line
	// returns : number of pairs evaluated
	int printAllPairs(string word) {
		int tot_pairs = 0;
		for ( auto it = seed_hash.begin() ; it != seed_hash.end() ; ++it ) {
			tot_pairs += it->second.size() * (it->second.size() - 1) / 2;		
			for ( int i = 0 ; i < it->second.size() ; i++ ) {
			for ( int ii = i+1 ; ii < it->second.size() ; ii++ ) {
				cout << all_ind[it->second[i]].getID() << "\t" << all_ind[it->second[ii]].getID() << "\t" << word << endl;
			}			
			}
		}
		return tot_pairs;
	}	
};

double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

int main (int argc, char* argv[]) 
{	

    cout << R"(
	                      ___          ___ 
	  ___ ____ ______ _  / (_)__  ___ |_  |
	 / _ `/ -_) __/  ' \/ / / _ \/ -_) __/ 
	 \_, /\__/_/ /_/_/_/_/_/_//_/\__/____/ 
	/___/                                  
)" << endl;


	double TIME_start = get_cpu_time();
	double TIME_prev = TIME_start;
	float PAR_MIN_MAF = 0;
	float PAR_skip = 0;
	
	string line, discard;
	bool opt_error = 0;
	int c;
	// load switches
	while (!opt_error && (c = getopt (argc, argv, "hbd:s:g:f:m:")) != -1)
		switch (c)
		{
		  case 'b':
			PAR_BIN_OUT = 1;
			break;
		  case 'h':
			PAR_HAPLOID = 1;
			break;
		  case 'm':
			PAR_MIN_MATCH = atof( optarg );
			break;
		  case 'f':
			PAR_MIN_MAF = atof( optarg );
			break;			
		  case 'g':
			PAR_GAP = atoi( optarg );
			break;
		  case 's':
			PAR_skip = atof( optarg );
			break;			
		  case 'd':
			MAX_seeds = atoi( optarg );
			break;	
		  default:
		  	opt_error = 1;
		  }

	cout << endl << "Options:" << endl;
	cout << "---" << endl;
	cout << "-b\tBinary output [default off]                                                " << PAR_BIN_OUT << endl;
	cout << "-d\tDynamic hash seed cutoff [default 0/off]                                   " << MAX_seeds << endl;
	cout << "-f\tMinimum minor allele frequency [default 0.0]                               " << PAR_MIN_MAF << endl;
	cout << "-g\tAllowed gaps [default 1]                                                   " << PAR_GAP << endl;
	cout << "-h\tHaploid mode, do not allow switches between haplotypes [default off]       " << PAR_HAPLOID << endl;
	cout << "-m\tMinimum match length [default 1.0]                                         " << PAR_MIN_MATCH << endl;
	cout << "-s\tSkip words with (seeds/samples) less than than this value [default 0.0]    " << PAR_skip << endl;
	cout << endl;
	if ( opt_error == 1 ) abort ();
	
	// load parameters
	if(opt_error || argc - optind != 4){
		cerr << "ERROR: Incorrect number of parameters" << endl;
		cerr << "Usage: g2 [options] <haps file> <sample file> <genetic map file> <output file>" << endl;			
		return 0;
	}
	
	cout << "\thaps: " << argv[optind + 0] << endl;
	cout << "\tsample: " << argv[optind + 1] << endl;
	cout << "\tgmap: " << argv[optind + 2] << endl;
	cout << "\toutput: " << argv[optind + 3] << endl << endl;
	cout << "---" << endl << endl;
		
	ifstream file_haps(argv[optind + 0]);
	ifstream file_samp(argv[optind + 1]);
	ifstream file_genm(argv[optind + 2]);
	
	string out = string( argv[optind + 3] );
	if ( PAR_BIN_OUT ) {
		FOUT.open( ( out + ".bmatch" ).c_str() , ios::binary );
	} else {
		FOUT.open( argv[optind + 3] );
	}
	
	if(!FOUT ) { cerr << argv[optind + 3] << " could not be opened" << endl; return -1; }
	if(!file_haps ) { cerr << argv[optind + 0] << " could not be opened" << endl; return -1; }
	if(!file_samp ) { cerr << argv[optind + 1] << " could not be opened" << endl; return -1; }
	if(!file_genm ) { cerr << argv[optind + 2] << " could not be opened" << endl; return -1; }
	
	string map_field[3];
	stringstream ss;
	
	// *** read genetic map
	vector< pair<long,double> > genetic_map;
	int cur_g = 0;
	while(getline(file_genm,line)) {
		ss.clear(); ss.str(line);
		ss >> map_field[0] >> map_field[1] >> map_field[2];
		if ( map_field[0] == "position" || map_field[0] == "" ) continue;
		genetic_map.push_back( pair<long,double>( stol(map_field[0]) , stod(map_field[2]) ) );
		if ( cur_g > 0 && (genetic_map[ cur_g ].first < genetic_map[ cur_g - 1 ].first || genetic_map[ cur_g ].second < genetic_map[ cur_g - 1 ].second) ) {
			cerr << "ERROR: genetic map not in sorted order at line\n" << line << endl;
			return -1;
		}
		cur_g++;
	}
	file_genm.close();
	if ( genetic_map.size() < 2 ) { cerr << "ERROR: genetic map must have at least two valid entries" << endl; return -1; }

	cerr << "*** runtime : " << get_cpu_time() - TIME_start << "\t";
	cerr << genetic_map.size() << " genetic map entries read" << endl;
	
	// *** read individual information
	// skip first two lines
	getline(file_samp,line);
	getline(file_samp,line);

	int idctr = 0;
	while(getline(file_samp,line)) {
		ss.clear(); ss.str( line );
		ss >> map_field[0] >> map_field[1];
		if ( PAR_HAPLOID ) {
			all_ind.push_back( Individual(map_field[0],(map_field[1]+".0").c_str(),idctr) );
			idctr++;		
			all_ind.push_back( Individual(map_field[0],(map_field[1]+".1").c_str(),idctr) );
		} else {
			all_ind.push_back( Individual(map_field[0],map_field[1],idctr) );
			all_ind.push_back( Individual(map_field[0],map_field[1],idctr) );
		}
		idctr++;		
	}
	file_samp.close();
	num_ind = all_ind.size();
	
	cerr << "*** runtime : " << get_cpu_time() - TIME_start << "\t";
	cerr << num_ind / 2 << " sample identifiers read" << endl;
	
	Marker cur_marker;
	// track position through genetic map
	cur_g = 0;
	int snp_ctr;
	char al[2] , inp;
	string cur_al;
	
	// Storage for seeds
	SeedHash seeds;
	
	hash_size word[2];
	
	// Storage for extensions
	ExtendHash extend( PAR_HAPLOID ? num_ind : num_ind / 2 );

	// Hash individual words
	GLOBAL_READ_WORDS = 0;
	GLOBAL_CURRENT_WORD = 0;
	
	while( 1 ) {
		snp_ctr = 0;
		while( getline(file_haps,line) )
		{
			// read the meta data
			ss.clear(); ss.str( line );
			ss >> map_field[0] >> cur_marker.id >> cur_marker.pos >> al[0] >> al[1];
			if(map_field[0] == "") continue;
			// loop until we reach this marker
			while ( cur_marker.pos > genetic_map[cur_g].first && cur_g < genetic_map.size() - 1 ) cur_g++;
			if ( cur_marker.pos >= genetic_map[cur_g].first ) {
				// we found this exact marker, or we reached the end of the map
				cur_marker.cm = genetic_map[cur_g].second;
			} else if ( cur_g == 0 ) {
				// if we haven't hit the map yet, store first map entry
				cur_marker.cm = genetic_map[cur_g].second;
			} else {
				// interpolate from previous marker
				cur_marker.cm = genetic_map[cur_g-1].second + (cur_marker.pos - genetic_map[cur_g-1].first) * ( genetic_map[cur_g].second - genetic_map[cur_g-1].second ) / ( genetic_map[cur_g].first - genetic_map[cur_g-1].first );
				// cerr << "interpolating " << cur_marker.id << " : " << cur_marker.pos << " between " << genetic_map[cur_g-1].first << "," << genetic_map[cur_g-1].second << " and " << genetic_map[cur_g].first << "," << genetic_map[cur_g].second << " to " << cur_marker.cm << endl;
			}
		
			// restrict on MAF
			if ( PAR_MIN_MAF > 0 ) {
				int maf_ctr = 0;
				for( int i=0; i< num_ind ; i++) {
					ss >> inp;
					if ( inp == '1' ) maf_ctr++;
				}
				float maf = (float) maf_ctr / num_ind;
				if ( maf < PAR_MIN_MAF || maf > 1 - PAR_MIN_MAF ) continue;

				// re-load the data
				ss.clear(); ss.str( line );
				ss >> map_field[0] >> map_field[0] >> map_field[0] >> al[0] >> al[1];
			}
			
			all_markers.push_back( cur_marker );
		
			// read haplotype
			for( unsigned int i=0; i<num_ind ; i++ )
			{
				ss >> inp;
				if ( inp == '1' ) all_ind[i].setMarker( GLOBAL_READ_WORDS , snp_ctr );
			}
			snp_ctr++;
			
			if ( snp_ctr % WORD_SIZE == 0 ) {
				// cerr << "*** word " << GLOBAL_CURRENT_WORD << " " << cur_marker.pos << " " << cur_marker.cm << endl;
				if ( ++GLOBAL_READ_WORDS >= CONST_READ_AHEAD ) break;
				else cerr << "*** loading word buffer " << GLOBAL_READ_WORDS << " / " << CONST_READ_AHEAD << endl;
				snp_ctr = 0 ;
			}
		}
		
		// end if read all data
		if( GLOBAL_CURRENT_WORD >= GLOBAL_READ_WORDS ) {
			cerr << "processed " << GLOBAL_CURRENT_WORD * WORD_SIZE << " / " << all_markers.size() << " SNPs" << endl;
			break;
		}
	
		for ( unsigned int i = 0 ; i < num_ind ; i++ ) {
			// cerr << i << "\t" << all_ind[i].getWordString( w ) << endl;
            seeds.insertIndividual( i , all_ind[i].getWordHash( GLOBAL_CURRENT_WORD ) );
		}

		GLOBAL_SKIPPED_WORDS = 0;
		int cur_seeds = seeds.size();
		unsigned long cur_pairs = 0;

		// skip low-complexity words
		if ( (float) cur_seeds / num_ind > PAR_skip ) {
			cur_pairs = seeds.extendAllPairs( &extend , GLOBAL_CURRENT_WORD );
			extend.clearPairsPriorTo( GLOBAL_CURRENT_WORD - PAR_GAP );
		} else {
			cerr << "low complexity word - " << cur_seeds << " - skipping" << endl;
			extend.extendAllPairsTo( GLOBAL_CURRENT_WORD );
		}


		cerr << "*** runtime : " << get_cpu_time() - TIME_start << "\t"
			<< GLOBAL_CURRENT_WORD << "\t"
			<< cur_seeds << " main seeds, "
			<< cur_pairs << " pairs, "
			<< extend.size() << " active seeds, "
			<< GLOBAL_SKIPPED_WORDS << " skipped words" << endl;
		seeds.clear();
		
		for( unsigned int i=0; i< num_ind ; i++) all_ind[i].clear( GLOBAL_CURRENT_WORD );
		GLOBAL_CURRENT_WORD++;
	}
	extend.clearAllPairs();
	file_haps.close();	
	FOUT.close();
	
	if ( PAR_BIN_OUT ) {
		ofstream bmid_out( ( out + ".bmid" ).c_str() );
		for ( int i = 0 ; i < all_markers.size() ; i++ ) bmid_out << all_markers[i].print();
		bmid_out.close();
		
		ofstream bsid_out( ( out + ".bsid" ).c_str() );
		for ( int i = 0 ; i < all_ind.size() ; i++ ) {
			bsid_out << all_ind[i].print();
			if ( !PAR_HAPLOID ) i++;
		}
		bsid_out.close();
		
	}
	cerr << "*** runtime : " << get_cpu_time() - TIME_start << endl;
		
	return 0;
}
