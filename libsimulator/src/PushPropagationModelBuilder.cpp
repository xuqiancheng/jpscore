// SPDX-License-Identifier: LGPL-3.0-or-later
#include "PushPropagationModelBuilder.hpp"

PushPropagationModelBuilder::PushPropagationModelBuilder(
    double pushoutStrength,
    uint64_t rng_seed)
    : pushoutStrength(pushoutStrength), rng_seed(rng_seed)
{
}

PushPropagationModel PushPropagationModelBuilder::Build()
{
    return PushPropagationModel(pushoutStrength, rng_seed);
}
