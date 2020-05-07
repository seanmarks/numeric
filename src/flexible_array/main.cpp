
#include <cstddef>
#include <iostream>
#include <vector>

#include "CommonTypes.h"
#include "VectorOfArrays.h"

int main(int argc, char* argv[])
{
	unsigned size = 15;
	VectorOfArrays<double> vec(size);
	vec.printLayout();
	std::cout << std::endl;

	vec.resize(5);
	vec.printLayout();
	std::cout << std::endl;

	vec.resize(20);
	vec.printLayout();
	std::cout << std::endl;

	vec.resize(3);
	vec.printLayout();
	std::cout << std::endl;

	vec.resize(7);
	vec.printLayout();
	std::cout << std::endl;

	/*
	std::cout << "size = " << vec.size() << " (.data().size() = " << vec.data().size() << ")\n";
	for ( unsigned i=0; i<vec.size(); ++i ) {
		std::cout << "  capacity[" << i << "] = " << vec.getArrayCapacity(i) << "\n";
	}
	*/

	return 0;
}
