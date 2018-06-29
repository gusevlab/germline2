// 1. for every true segment, how much of it is detected?
// 2. for every detected segment, how much of it is true?

#include <iostream>
#include <sstream>
#include <string>
#include <list>

using namespace std;

struct Match {
	string line;
	long pos[2];
	long match;
	// most overlapping segment
	long best_dif = 0, best_len = 0;
	Match() { match = 0; }
	void addOverlap( long d , long len );
};

void Match::addOverlap( long d , long len ) {
	match += d;
	if ( d > best_dif ) {
		best_dif = d;
		best_len = len;
	}
}

int main (int argc, char* argv[]) 
{

	// cin should be <0/1> <pair> <start> <end>
	string line , cur , prev_pair , pair;
	int pre;

	Match cur_match;
	list< Match > segments[2];
	list< Match >::iterator segments_it[2];
	
	stringstream ss;
	prev_pair = "";
	
	long ss_total[2];
	bool last_pair = false;
	
	while( true ) {
		if ( last_pair || getline( cin , line ) ) {
		
			if ( ! last_pair ) {
				ss.clear(); ss.str( line );
				cur_match.line = line;
				ss >> pre >> pair >> cur_match.pos[0] >> cur_match.pos[1];
			}
			
			if ( last_pair || pair != prev_pair ) {					
				if ( prev_pair != "" ) {
					// compute overlap
					for ( segments_it[0] = segments[0].begin() ; segments_it[0] != segments[0].end(); segments_it[0]++ )
					{
						for ( segments_it[1] = segments[1].begin() ; segments_it[1] != segments[1].end(); segments_it[1]++ )
						{
							long dif = 0;
							if ( (*segments_it[1]).pos[1] > (*segments_it[0]).pos[0] && (*segments_it[1]).pos[0] < (*segments_it[0]).pos[1] )
							{
								if ( (*segments_it[1]).pos[1] < (*segments_it[0]).pos[1] )
								{
									if ( (*segments_it[1]).pos[0] > (*segments_it[0]).pos[0] )
										dif = (*segments_it[1]).pos[1] - (*segments_it[1]).pos[0];
									else
										dif = (*segments_it[1]).pos[1] - (*segments_it[0]).pos[0];
								} else
								if ( (*segments_it[1]).pos[0] > (*segments_it[0]).pos[0] )
								{
									if ( (*segments_it[1]).pos[1] < (*segments_it[0]).pos[1] )
										dif = (*segments_it[1]).pos[1] - (*segments_it[1]).pos[0];
									else
										dif = (*segments_it[0]).pos[1] - (*segments_it[1]).pos[0];
								} else dif = (*segments_it[0]).pos[1] - (*segments_it[0]).pos[0];
							}
							(*segments_it[0]).addOverlap( dif , (*segments_it[1]).pos[1] - (*segments_it[1]).pos[0] );
							(*segments_it[1]).addOverlap( dif , (*segments_it[0]).pos[1] - (*segments_it[0]).pos[0] );
						}
					}

					for ( segments_it[0] = segments[0].begin() ; segments_it[0] != segments[0].end(); segments_it[0]++ ) {
						cout << "0\t" << (float) (*segments_it[0]).match / ((*segments_it[0]).pos[1] - (*segments_it[0]).pos[0]) << '\t' << (*segments_it[0]).line << endl;
						cout << "0_cor\t" << (*segments_it[0]).best_len << '\t' << ((*segments_it[0]).pos[1] - (*segments_it[0]).pos[0]) << endl;
					}
					for ( segments_it[0] = segments[1].begin() ; segments_it[0] != segments[1].end(); segments_it[0]++ )
						cout << "1\t" << (float) (*segments_it[0]).match / ((*segments_it[0]).pos[1] - (*segments_it[0]).pos[0]) << '\t' << (*segments_it[0]).line << endl;
					cout << "total" << '\t' << prev_pair << '\t' << ss_total[0] << '\t' << ss_total[1] << endl;
				}
				// clear out all the matches
				ss_total[0] = ss_total[1] = 0;
				segments[0].clear();
				segments[1].clear();
			}
			
			if ( last_pair ) break;
			
			segments[ pre ].push_back( cur_match );
			ss_total[ pre ] += cur_match.pos[1] - cur_match.pos[0];
			prev_pair = pair;
		} else {
			last_pair = true;
		}
	}
	return 1;
}	
