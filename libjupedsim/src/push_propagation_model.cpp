// SPDX-License-Identifier: LGPL-3.0-or-later
#include "jupedsim/push_propagation_model.h"
#include "jupedsim/error.h"

#include "Conversion.hpp"
#include "ErrorMessage.hpp"

#include <PushPropagationModel.hpp>
#include <PushPropagationModelBuilder.hpp>
#include <PushPropagationModelData.hpp>

#include <stdint.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Push Propagation Model Builder
////////////////////////////////////////////////////////////////////////////////////////////////////
JUPEDSIM_API JPS_PushPropagationModelBuilder
JPS_PushPropagationModelBuilder_Create(double pushoutStrength, uint64_t rng_seed)
{
    return reinterpret_cast<JPS_PushPropagationModelBuilder>(
        new PushPropagationModelBuilder(pushoutStrength, rng_seed));
}

JUPEDSIM_API JPS_OperationalModel JPS_PushPropagationModelBuilder_Build(
    JPS_PushPropagationModelBuilder handle,
    JPS_ErrorMessage* errorMessage)
{
    assert(handle != nullptr);
    auto builder = reinterpret_cast<PushPropagationModelBuilder*>(handle);
    JPS_OperationalModel result{};
    try {
        result = reinterpret_cast<JPS_OperationalModel>(new PushPropagationModel(builder->Build()));
    } catch(const std::exception& ex) {
        if(errorMessage) {
            *errorMessage = reinterpret_cast<JPS_ErrorMessage>(new JPS_ErrorMessage_t{ex.what()});
        }
    } catch(...) {
        if(errorMessage) {
            *errorMessage = reinterpret_cast<JPS_ErrorMessage>(
                new JPS_ErrorMessage_t{"Unknown internal error."});
        }
    }
    return result;
}

JUPEDSIM_API void JPS_PushPropagationModelBuilder_Free(JPS_PushPropagationModelBuilder handle)
{
    delete reinterpret_cast<PushPropagationModelBuilder*>(handle);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// PushPropagationModelState
////////////////////////////////////////////////////////////////////////////////////////////////////
double JPS_PushPropagationModelState_GetStrengthNeighborRepulsion(JPS_PushPropagationModelState handle)
{
    assert(handle);
    const auto state = reinterpret_cast<const PushPropagationModelData*>(handle);
    return state->strengthNeighborRepulsion;
}

void JPS_PushPropagationModelState_SetStrengthNeighborRepulsion(
    JPS_PushPropagationModelState handle,
    double strengthNeighborRepulsion)
{
    assert(handle);
    auto state = reinterpret_cast<PushPropagationModelData*>(handle);
    state->strengthNeighborRepulsion = strengthNeighborRepulsion;
}

double JPS_PushPropagationModelState_GetRangeNeighborRepulsion(JPS_PushPropagationModelState handle)
{
    assert(handle);
    const auto state = reinterpret_cast<const PushPropagationModelData*>(handle);
    return state->rangeNeighborRepulsion;
}

void JPS_PushPropagationModelState_SetRangeNeighborRepulsion(
    JPS_PushPropagationModelState handle,
    double rangeNeighborRepulsion)
{
    assert(handle);
    auto state = reinterpret_cast<PushPropagationModelData*>(handle);
    state->rangeNeighborRepulsion = rangeNeighborRepulsion;
}

double JPS_PushPropagationModelState_GetAnticipationTime(JPS_PushPropagationModelState handle)
{
    assert(handle);
    const auto state = reinterpret_cast<const PushPropagationModelData*>(handle);
    return state->anticipationTime;
}

void JPS_PushPropagationModelState_SetAnticipationTime(
    JPS_PushPropagationModelState handle,
    double anticipationTime)
{
    assert(handle);
    auto state = reinterpret_cast<PushPropagationModelData*>(handle);
    state->anticipationTime = anticipationTime;
}

double JPS_PushPropagationModelState_GetReactionTime(JPS_PushPropagationModelState handle)
{
    assert(handle);
    const auto state = reinterpret_cast<const PushPropagationModelData*>(handle);
    return state->reactionTime;
}

void JPS_PushPropagationModelState_SetReactionTime(
    JPS_PushPropagationModelState handle,
    double reactionTime)
{
    assert(handle);
    auto state = reinterpret_cast<PushPropagationModelData*>(handle);
    state->reactionTime = reactionTime;
}

double JPS_PushPropagationModelState_GetWallBufferDistance(JPS_PushPropagationModelState handle)
{
    assert(handle);
    const auto state = reinterpret_cast<const PushPropagationModelData*>(handle);
    return state->wallBufferDistance;
}

void JPS_PushPropagationModelState_SetWallBufferDistance(
    JPS_PushPropagationModelState handle,
    double wallBufferDistance)
{
    assert(handle);
    auto state = reinterpret_cast<PushPropagationModelData*>(handle);
    state->wallBufferDistance = wallBufferDistance;
}

double JPS_PushPropagationModelState_GetTimeGap(JPS_PushPropagationModelState handle)
{
    assert(handle);
    const auto state = reinterpret_cast<const PushPropagationModelData*>(handle);
    return state->timeGap;
}

void JPS_PushPropagationModelState_SetTimeGap(JPS_PushPropagationModelState handle,
    double time_gap)
{
    assert(handle);
    auto state = reinterpret_cast<PushPropagationModelData*>(handle);
    state->timeGap = time_gap;
}

double JPS_PushPropagationModelState_GetV0(JPS_PushPropagationModelState handle)
{
    assert(handle);
    const auto state = reinterpret_cast<const PushPropagationModelData*>(handle);
    return state->v0;
}

void JPS_PushPropagationModelState_SetV0(JPS_PushPropagationModelState handle, double v0)
{
    assert(handle);
    auto state = reinterpret_cast<PushPropagationModelData*>(handle);
    state->v0 = v0;
}
double JPS_PushPropagationModelState_GetRadius(JPS_PushPropagationModelState handle)
{
    assert(handle);
    const auto state = reinterpret_cast<const PushPropagationModelData*>(handle);
    return state->radius;
}

void JPS_PushPropagationModelState_SetRadius(JPS_PushPropagationModelState handle,
    double radius)
{
    assert(handle);
    auto state = reinterpret_cast<PushPropagationModelData*>(handle);
    state->radius = radius;
}
