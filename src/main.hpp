
#ifndef SC_MAIN_H
#define SC_MAIN_H

#include "common.hpp"

int main_initialization(SCinstance &inst);
int main_read_params(SCinstance &inst, int argc, char *argv[]);
int main_read_instance_dns(SCinstance &inst);
int main_read_instance_spr(SCinstance &inst);

#endif //SC_MAIN_H
