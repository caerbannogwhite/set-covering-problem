
#include "main.hpp"
#include "sc.hpp"

int main(int argc, char *argv[])
{
	SCinstance inst;

	comm_initialization(inst);
	if (comm_read_params(inst, argc, argv))
	{
		std::printf("Type --help to display available options.\n");
		return 0;
	}

	sc_solver(inst);

	return 0;
}