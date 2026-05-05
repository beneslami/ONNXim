#include "ThreeDICEWrapper.h"
#include <stdexcept>
#include <filesystem>
#include <spdlog/spdlog.h>

extern "C" {
#include "stack_description.h"
#include "stack_file_parser.h"
#include "analysis.h"
#include "output.h"
#include "thermal_data.h"
#include "powers_queue.h"
#include "floorplan.h"
}

namespace fs = std::filesystem;

struct ThreeDICEImpl {
    StackDescription_t stkd;
    Analysis_t         analysis;
    Output_t           output;
    ThermalData_t      tdata;
    uint32_t           n_flp_elements = 0;
    std::string        original_dir;
};

ThreeDICEWrapper::ThreeDICEWrapper(const std::string& stk_path) : _impl(std::make_unique<ThreeDICEImpl>())
{
    _impl->original_dir = fs::current_path().string();
    fs::current_path(fs::path(stk_path).parent_path());
    std::string stk_filename = fs::path(stk_path).filename().string();

    stack_description_init(&_impl->stkd);
    analysis_init         (&_impl->analysis);
    output_init           (&_impl->output);

    if (parse_stack_description_file((String_t)stk_filename.c_str(), &_impl->stkd, &_impl->analysis, &_impl->output) != TDICE_SUCCESS) {
        fs::current_path(_impl->original_dir);
        throw std::runtime_error("3D-ICE: failed to parse " + stk_filename);
    }

    generate_output_headers(&_impl->output, _impl->stkd.Dimensions, (String_t)"% ");

    thermal_data_init(&_impl->tdata);

    if (thermal_data_build(&_impl->tdata, &_impl->stkd.StackElements, _impl->stkd.Dimensions, &_impl->analysis, &_impl->stkd.Materials) != TDICE_SUCCESS) {
        fs::current_path(_impl->original_dir);
        throw std::runtime_error("3D-ICE: thermal_data_build failed");
    }

    fs::current_path(_impl->original_dir);
    _impl->n_flp_elements = get_total_number_of_floorplan_elements(&_impl->stkd);
    spdlog::info("[3DICE] Initialized: {} floorplan elements", _impl->n_flp_elements);
}

ThreeDICEWrapper::~ThreeDICEWrapper() {
    thermal_data_destroy      (&_impl->tdata);
    stack_description_destroy (&_impl->stkd);
    output_destroy            (&_impl->output);
}

double ThreeDICEWrapper::getStepTimeSec() const {
    return (double)_impl->analysis.StepTime;
}

uint32_t ThreeDICEWrapper::getNumFlpElements() const {
    return _impl->n_flp_elements;
}

std::vector<double> ThreeDICEWrapper::computeTemperatures(const std::vector<double>& power_w) {
    if ((uint32_t)power_w.size() != _impl->n_flp_elements)
        throw std::runtime_error("3D-ICE: wrong power vector size");

    PowersQueue_t pq;
    powers_queue_init(&pq);
    powers_queue_build(&pq, _impl->n_flp_elements);
    for (double p : power_w)
        put_into_powers_queue(&pq, (Power_t)p);

    insert_power_values(&_impl->tdata.PowerGrid, &pq);
    powers_queue_destroy(&pq);

    SimResult_t result = emulate_step(&_impl->tdata, _impl->stkd.Dimensions, &_impl->analysis);
    if (result == TDICE_SOLVER_ERROR)
        throw std::runtime_error("3D-ICE: solver error");

    std::vector<double> temps;
    std::vector<StackElement_t*> stkel_vec;

    for (StackElementListNode_t *node = stack_element_list_begin(&_impl->stkd.StackElements);
         node != NULL;
         node = stack_element_list_next(node))
        stkel_vec.push_back(stack_element_list_data(node));

    for (int i = (int)stkel_vec.size()-1; i >= 0; i--) {
        StackElement_t *stkel = stkel_vec[i];
        if (stkel->SEType != TDICE_STACK_ELEMENT_DIE) continue;

        Floorplan_t *flp = &stkel->Pointer.Die->Floorplan;
        CellIndex_t layer_area = get_layer_area(_impl->stkd.Dimensions);
        CellIndex_t offset = (stkel->Offset + stkel->Pointer.Die->SourceLayerOffset) * layer_area;
        Quantity_t n = 0;

        Temperature_t *tmax = get_all_max_temperatures_floorplan(flp, _impl->stkd.Dimensions, _impl->tdata.Temperatures + offset, &n, NULL);

        for (Quantity_t j = 0; j < n; j++)
            temps.push_back((double)tmax[j] - 273.15);
        free(tmax);
    }
    return temps;
}