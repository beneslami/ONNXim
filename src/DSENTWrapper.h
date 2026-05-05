// Developed by Benyamin Eslami, 1 May, 2026 - Washington State University
#pragma once

#pragma once

#include <map>
#include <string>
#include <filesystem>
#include <unistd.h>
#include <fcntl.h>
#include "DSENT.h"
#include "libutil/String.h"

namespace fs = std::filesystem;

class DSENTWrapper {
public:
    DSENTWrapper(const std::string& config_path) {
        // Save current working directory
        std::string original_dir = fs::current_path().string();

        // Change to config file's directory so relative paths resolve correctly
        std::string config_dir = fs::path(config_path).parent_path().string();
        fs::current_path(config_dir);

        // Initialize using just the filename (relative paths now work)
        std::string config_filename = fs::path(config_path).filename().string();
        _ms_model = DSENT::initialize(config_filename.c_str(), _config);

        // Restore original working directory
        fs::current_path(original_dir);
    }

    ~DSENTWrapper() {
        DSENT::finalize(_config, _ms_model);
    }

    double computePower(double injection_rate) {
        _config["InjectionRate"] = LibUtil::String(injection_rate);
        std::map<std::string, double> outputs;
        // Suppress DSENT's print statement console output
        int saved_stdout = dup(STDOUT_FILENO);
        int devnull = open("/dev/null", O_WRONLY);
        dup2(devnull, STDOUT_FILENO);
        close(devnull);

        DSENT::run(_config, _ms_model, outputs);

        // Restore stdout
        dup2(saved_stdout, STDOUT_FILENO);
        close(saved_stdout);
        auto it = outputs.find("total");
        return (it != outputs.end()) ? it->second : 0.0;
    }

    void printSummary(double avg_injection_rate) {
        auto summary_it = _config.find("EvaluateStringSummary");
        if (summary_it == _config.end()) {
            spdlog::warn("[DSENT] EvaluateStringSummary not found in config");
            return;
        }

        // Swap EvaluateString with EvaluateStringSummary temporarily
        LibUtil::String original = _config["EvaluateString"];
        _config["EvaluateString"] = summary_it->second;
        _config["InjectionRate"] = LibUtil::String(avg_injection_rate);

        std::map<std::string, double> outputs;
        DSENT::run(_config, _ms_model, outputs);

        // Restore original
        _config["EvaluateString"] = original;
    }

private:
    std::map<LibUtil::String, LibUtil::String> _config;
    DSENT::Model* _ms_model = nullptr;
};