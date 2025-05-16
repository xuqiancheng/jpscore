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

double
JPS_PushPropagationModelState_GetExternalForce(JPS_PushPropagationModelState handle)
{
    assert(handle);
    const auto state = reinterpret_cast<const PushPropagationModelData*>(handle);
    return state->externalForce;
}

void JPS_PushPropagationModelState_SetExternalForce(
    JPS_PushPropagationModelState handle,
    double external_force)
{
    assert(handle);
    auto state = reinterpret_cast<PushPropagationModelData*>(handle);
    state->externalForce = external_force;
}


double JPS_PushPropagationModelState_GetWeight(JPS_PushPropagationModelState handle)
{
    assert(handle);
    const auto state = reinterpret_cast<const PushPropagationModelData*>(handle);
    return state->weight;
}

void JPS_PushPropagationModelState_SetWeight(
    JPS_PushPropagationModelState handle,
    double weight)
{
    assert(handle);
    auto state = reinterpret_cast<PushPropagationModelData*>(handle);
    state->weight = weight;
}

double JPS_PushPropagationModelState_GetHeight(JPS_PushPropagationModelState handle)
{
    assert(handle);
    const auto state = reinterpret_cast<const PushPropagationModelData*>(handle);
    return state->height;
}

void JPS_PushPropagationModelState_SetHeight(JPS_PushPropagationModelState handle, double height)
{
    assert(handle);
    auto state = reinterpret_cast<PushPropagationModelData*>(handle);
    state->height = height;
}

JUPEDSIM_API double
JPS_PushPropagationModelState_GetForceAmplify(JPS_PushPropagationModelState handle)
{
    assert(handle);
    const auto state = reinterpret_cast<const PushPropagationModelData*>(handle);
    return state->forceAmplify;
}

JUPEDSIM_API void JPS_PushPropagationModelState_SetForceAmplify(
    JPS_PushPropagationModelState handle,
    double forceAmplify)
{
    assert(handle);
    auto state = reinterpret_cast<PushPropagationModelData*>(handle);
    state->forceAmplify = forceAmplify;
}

JUPEDSIM_API double
JPS_PushPropagationModelState_GetPushReactionTime(JPS_PushPropagationModelState handle)
{
    assert(handle);
    const auto state = reinterpret_cast<const PushPropagationModelData*>(handle);
    return state->pushReactionTime;
}

JUPEDSIM_API void JPS_PushPropagationModelState_SetPushReactionTime(
    JPS_PushPropagationModelState handle,
    double pushReactionTime)
{
    assert(handle);
    auto state = reinterpret_cast<PushPropagationModelData*>(handle);
    state->pushReactionTime = pushReactionTime;
}

JUPEDSIM_API double JPS_PushPropagationModelState_GetDeceConst(JPS_PushPropagationModelState handle)
{
    assert(handle);
    const auto state = reinterpret_cast<const PushPropagationModelData*>(handle);
    return state->deceConst;
}

JUPEDSIM_API void
JPS_PushPropagationModelState_SetDeceConst(JPS_PushPropagationModelState handle, double deceConst)
{
    assert(handle);
    auto state = reinterpret_cast<PushPropagationModelData*>(handle);
    state->deceConst = deceConst;
}

JUPEDSIM_API double JPS_PushPropagationModelState_GetDeceVelo(JPS_PushPropagationModelState handle)
{
    assert(handle);
    const auto state = reinterpret_cast<const PushPropagationModelData*>(handle);
    return state->deceVelo;
}

JUPEDSIM_API void
JPS_PushPropagationModelState_SetDeceVelo(JPS_PushPropagationModelState handle, double deceVelo)
{
    assert(handle);
    auto state = reinterpret_cast<PushPropagationModelData*>(handle);
    state->deceVelo = deceVelo;
}

JUPEDSIM_API double JPS_PushPropagationModelState_GetKpush(JPS_PushPropagationModelState handle)
{
    assert(handle);
    const auto state = reinterpret_cast<const PushPropagationModelData*>(handle);
    return state->Kpush;
}

JUPEDSIM_API void
JPS_PushPropagationModelState_SetKpush(JPS_PushPropagationModelState handle, double kpush)
{
    assert(handle);
    auto state = reinterpret_cast<PushPropagationModelData*>(handle);
    state->Kpush = kpush;
}

JUPEDSIM_API double JPS_PushPropagationModelState_GetDpush(JPS_PushPropagationModelState handle)
{
    assert(handle);
    const auto state = reinterpret_cast<const PushPropagationModelData*>(handle);
    return state->Dpush;
}

JUPEDSIM_API void
JPS_PushPropagationModelState_SetDpush(JPS_PushPropagationModelState handle, double dpush)
{
    assert(handle);
    auto state = reinterpret_cast<PushPropagationModelData*>(handle);
    state->Dpush = dpush;
}
