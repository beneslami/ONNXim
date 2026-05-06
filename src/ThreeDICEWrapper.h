#pragma once
#include <string>
#include <vector>
#include <memory>

// Forward declaration — hides all 3D-ICE internals
struct ThreeDICEImpl;

struct ComponentThermalStats {
    std::string name;
    double temp_min   = std::numeric_limits<double>::max();
    double temp_max   = std::numeric_limits<double>::lowest();
    double temp_sum   = 0.0;
    uint64_t samples  = 0;

    void update(double temp_c) {
        temp_min  = std::min(temp_min, temp_c);
        temp_max  = std::max(temp_max, temp_c);
        temp_sum += temp_c;
        samples++;
    }

    double temp_avg() const {
        return samples > 0 ? temp_sum / samples : 0.0;
    }
};

class ThreeDICEWrapper {
public:
    explicit ThreeDICEWrapper(const std::string& stk_path);
    ~ThreeDICEWrapper();

    std::vector<double> computeTemperatures(const std::vector<double>& power_w);
    double   getStepTimeSec()     const;
    uint32_t getNumFlpElements()  const;
    void updateStats(const std::vector<double>& temps);
    void printFinalStats() const;
    const std::vector<ComponentThermalStats>& getStats() const;
private:
    std::unique_ptr<ThreeDICEImpl> _impl;
};