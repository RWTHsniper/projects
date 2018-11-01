
#include "helpers.h"

using namespace std;

namespace helpers{


string put_leading_zeros(unsigned int width,unsigned int number){

	std::string str_i = std::to_string(number);
	std::string leadingzero_i = std::string(width - str_i.length(), '0') + str_i;

	return leadingzero_i;
}


}
