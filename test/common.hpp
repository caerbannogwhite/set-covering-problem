
#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <time.h>
#include <set>
#include <unordered_set>
#include <vector>

#include <armadillo>
#include "catch.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

#define DEBUG_VERBOSITY 0

#define SC_EPSILON_MED 1e-5
#define SC_EPSILON_SMALL 1e-12
#define BIG_M 1e20

