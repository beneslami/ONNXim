#pragma once
#include <string>
#include <vector>
#include <memory>

// Forward declaration — hides all 3D-ICE internals
struct ThreeDICEImpl;

class ThreeDICEWrapper {
public:
    explicit ThreeDICEWrapper(const std::string& stk_path);
    ~ThreeDICEWrapper();

    std::vector<double> computeTemperatures(const std::vector<double>& power_w);
    double   getStepTimeSec()     const;
    uint32_t getNumFlpElements()  const;

private:
    std::unique_ptr<ThreeDICEImpl> _impl;
};