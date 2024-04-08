#include "util.h"
#include "gsum.h"

#include <algorithm>
#include <iostream>
#include <map>

std::string get_report(std::string algo, std::string dataset, std::string data_tag,
	double read, double merge, double encoding)
{
	char buff[233];
	snprintf(buff, sizeof buff,
		"%s, %s, %s\n"\
		"read: %.3f(s), merge: %.3f(s), encoding: %.3f(s)\n",
		algo.c_str(), dataset.c_str(), data_tag.c_str(),
		read, merge, encoding);
	return buff;
}

int main(int argc, char **argv) {
	if (argc == 2) {
		int round = 50, hash_num = 40;

		std::string algo_name = "Mags";
		std::string data_tag = "hash_num: " + std::to_string(hash_num) + ", round: " + std::to_string(round);

		// create graph from input path
		hclock::time_point s1 = hclock::now();
		Graph *g = from_text(argv[1]);
		
		// get the super-graph id of each node
		// sid[v] := the id of the super-node that contains v
		hclock::time_point s2 = hclock::now();
		VI sid = mags(*g, hash_num, round);
		
		// build optimal super-graph encoding from sid
		hclock::time_point s3 = hclock::now();
		SuperGraph sg(sid, *g);
		hclock::time_point s4 = hclock::now();

		// report time and relative size of the algorithm
		double read = get_duration(s1, s2);
		double merge = get_duration(s2, s3);
		double encoding = get_duration(s3, s4);
		std::cout << get_report(algo_name, argv[1], data_tag, read, merge, encoding);

		LL cost = sg.get_cost();
		LL en = g->edges_num_ * 2;
		printf("relative size: %lld/%lld = %.9f\n\n", cost, en, double(cost) / en);
	}
	return 0;
}