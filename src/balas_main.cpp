#include "balas_main.hpp"
#include "balas_common.hpp"
#include "balas_solver.hpp"

int main(int argc, char *argv[])
{
    BALSOLEnv inst;

    balcomm_initialization(inst);
    if (balcomm_read_params(inst, argc, argv))
    {
        std::printf("Type --help to display available options.\n");
        return 0;
    }

    balsol(inst);

    return 0;
}