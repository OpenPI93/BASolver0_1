#include "BAProblem.h"
#include <fstream>
#include <string>
#include <sophus\se3.hpp>

#include <thread>

using std::ifstream;
using std::string;

string p3d_file = "./p3d.txt";
string p2d_file = "./p2d.txt";

int main(int argc, char** argv)
{
	Sophus::Vector6d T = Sophus::SE3d().log();

	double K[4] = { 520.9, 521.0, 325.1, 249.7 };
	BASolver solver(K);
	baseCamera camera = {T[0], T[1], T[2], T[3], T[4], T[5] };
	solver.addCamera(camera);

	ifstream inp2p(p2d_file);

	if (!inp2p.is_open()) {
		cout << "can not find inp2p\n";
		return 1;
	}
	int index = 0;
	while (!inp2p.eof()) {
		observer edge;
		inp2p >> edge.u >> edge.v;
		edge.cam = 0;
		edge.pt = index;
		solver.addObserver(edge);
		++index;
	}
	inp2p.close();

	ifstream inp3p(p3d_file);
	if (!inp3p.is_open()) {
		cout << "can not find inp3p\n";
		return 1;
	}
	double fx = 520.9 * 520.9;
	double fy = 521.0 * 521.0;
	while (!inp3p.eof()) {
		point pt;
		inp3p >> pt[0] >> pt[1] >> pt[2];
		
		solver.addPoint(pt);
	}
	inp3p.close();

	solver.solveProblem(100);

#ifdef _MSC_VER
	system("pause");
#endif
}

