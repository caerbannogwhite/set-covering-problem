#include "cpx_main.hpp"
#include "cpx_solver.hpp"

int main(int argc, char *argv[])
{
	SCinstance inst;

	cpxcomm_initialization(inst);
	if (cpxcomm_read_params(inst, argc, argv))
	{
		std::printf("Type --help to display available options.\n");
		return 0;
	}

	cpxsol(inst);

	return 0;
}