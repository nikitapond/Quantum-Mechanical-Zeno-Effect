
#pragma ONCE

#include "RabiMCWF.h"
#include <vector>


int main(){

    RabiMCWF t;
    std::vector<float> res = t.RunExperiment(100, 0, 0,0,0);
    std::cout << std::to_string(res[200]) << std::endl;
    return 0;
}
