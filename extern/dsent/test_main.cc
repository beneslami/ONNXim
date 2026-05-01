#include <iostream>
#include <map>
#include <string>

#include "DSENT.h"
#include "libutil/String.h"

using namespace std;
using namespace LibUtil;

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <config_file>" << endl;
        cerr << "Example: " << argv[0] << " configs/router.cfg" << endl;
        return 1;
    }

    // Step 1: Initialize — loads config, builds model
    map<String, String> config;
    DSENT::Model* ms_model = DSENT::initialize(argv[1], config);

    cout << "=== DSENT Model Built Successfully ===" << endl;

    // Step 2: Run — processes EvaluateString queries from config
    map<string, double> outputs;
    DSENT::run(config, ms_model, outputs);

    // Step 3: Print outputs
    cout << "=== DSENT Outputs ===" << endl;
    for (const auto& it : outputs) {
        cout << it.first << " = " << it.second << endl;
    }

    // Step 4: Finalize — cleanup
    DSENT::finalize(config, ms_model);

    return 0;
}
