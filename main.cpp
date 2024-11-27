#include "WFC.hpp"

using namespace PSE::WFC;

int main(int argc, char *argv[])
{
    std::string image_path = "Tile2.png";
    std::string output_path = "generated_map.png";

    WaveFunctionCollapse::RunParameters Params{image_path, output_path, 4, 50};
    WaveFunctionCollapse::Run(Params);

    return 0;
}