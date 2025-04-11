// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once

#include "PushPropagationModel.hpp"

#include <stdint.h>

class PushPropagationModelBuilder
{
public:
    PushPropagationModelBuilder(double pushoutStrength, uint64_t rng_seed);
    PushPropagationModel Build();

private:
    double pushoutStrength{};
    uint64_t rng_seed{};
};
