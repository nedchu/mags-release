#include "util.h"
#include "pgsum.h"
#include <omp.h>

#include <algorithm>
#include <iostream>

std::string get_report(std::string algo, std::string dataset, std::string data_tag, int num_thread,
	double read, double merge, double encoding)
{
	char buff[233];
	snprintf(buff, sizeof buff,
		"%s, %s, %s, #thread: %d\n"\
		"read: %.3f(s), merge: %.3f(s), encoding: %.3f(s)\n",
		algo.c_str(), dataset.c_str(), data_tag.c_str(), num_thread,
		read, merge, encoding);
	return buff;
}

int main(int argc, char **argv) {
	if (argc == 2 || argc == 3) {
		int num_thread = 40;
		if (argc == 3) {
			num_thread = std::stoi(argv[2]);
		}
		omp_set_num_threads(num_thread);
		int round = 50, hash_num = 40;

		std::string algo = "Para-Mags-DM";
		std::string data_tag = "hash_num: " + std::to_string(hash_num) + ", round: " + std::to_string(round);

		hclock::time_point s1 = hclock::now();
		Graph *g = from_text(argv[1]);
		hclock::time_point s2 = hclock::now();
		VI sid = parallel_mags_dm(*g, hash_num, round);
		hclock::time_point s3 = hclock::now();
		ParaSuperGraph sg2(sid, *g);
		hclock::time_point s4 = hclock::now();

		LL cost = sg2.get_cost();
		LL en = g->edges_num_ * 2;

		double read = get_duration(s1, s2);
		double merge = get_duration(s2, s3);
		double encoding = get_duration(s3, s4);
		std::cout << get_report(algo, argv[1], data_tag, num_thread, read, merge, encoding);
		printf("relative size: %lld/%lld  %.9f\n\n", cost, en, double(cost) / en);
	}
	return 0;
}

