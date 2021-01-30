
#pragma once

#include <vector>
#include "RabiMCWF.h"


int main(){

    RabiMCWF t;
    std::vector<Vec2R> res = t.RunExperiment(1000, 500, 1,0.001,1000);
    // std::cout << std::to_string(real(res[200].x)) << std::endl;
    return 0;
}
