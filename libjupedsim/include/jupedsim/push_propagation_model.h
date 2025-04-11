// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once

#include "error.h"
#include "export.h"
#include "operational_model.h"
#include "types.h"

#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

/**
 * Opaque type for a Push Propagation Model Builder.
 */
typedef struct JPS_PushPropagationModelBuilder_t* JPS_PushPropagationModelBuilder;

/**
 * Creates a Push Propagation Model builder.
 * @param pushoutStrength strength of repulsive force of walls.
 * @param rng_seed Seed value for random number generator.
 * @return the builder
 */
JUPEDSIM_API JPS_PushPropagationModelBuilder
JPS_PushPropagationModelBuilder_Create(double pushoutStrength, uint64_t rng_seed);

/**
 * Creates a JPS_OperationalModel of type Push Propagation Model from the
 * JPS_PushPropagationModelBuilder.
 * @param handle the builder to operate on
 * @param[out] errorMessage if not NULL: will be set to a JPS_ErrorMessage in case of an error
 * @return a PS_PushPropagationModel or NULL if an error occurred.
 */
JUPEDSIM_API JPS_OperationalModel JPS_PushPropagationModelBuilder_Build(
    JPS_PushPropagationModelBuilder handle,
    JPS_ErrorMessage* errorMessage);

/**
 * Frees a JPS_PushPropagationModelBuilder
 * @param handle to the JPS_PushPropagationModelBuilder to free.
 */
JUPEDSIM_API void JPS_PushPropagationModelBuilder_Free(JPS_PushPropagationModelBuilder handle);

/**
 * Opaque type of Push Propagation Model state
 */
typedef struct JPS_PushPropagationModelState_t* JPS_PushPropagationModelState;

/**
 * Read strength neighbor repulsion of this agent.
 * @param handle of the Agent to access.
 * @return strength neighbor repulsion of this agent
 */
JUPEDSIM_API double JPS_PushPropagationModelState_GetStrengthNeighborRepulsion(JPS_PushPropagationModelState handle);

/**
 * Write strength neighbor repulsion of this agent.
 * @param handle of the Agent to access.
 * @param strengthNeighborRepulsion of this agent.
 */
JUPEDSIM_API void JPS_PushPropagationModelState_SetStrengthNeighborRepulsion(
    JPS_PushPropagationModelState handle,
    double strengthNeighborRepulsion);

/**
 * Read range neighbor repulsion of this agent.
 * @param handle of the Agent to access.
 * @return range neighbor repulsion of this agent
 */
JUPEDSIM_API double
JPS_PushPropagationModelState_GetRangeNeighborRepulsion(JPS_PushPropagationModelState handle);

/**
 * Write range neighbor repulsion of this agent.
 * @param handle of the Agent to access.
 * @param rangeNeighborRepulsion of this agent.
 */
JUPEDSIM_API void JPS_PushPropagationModelState_SetRangeNeighborRepulsion(
    JPS_PushPropagationModelState handle,
    double rangeNeighborRepulsion);

/**
 * Read anticipation time of this agent.
 * @param handle of the Agent to access.
 * @return anticipation time of this agent
 */
JUPEDSIM_API double
JPS_PushPropagationModelState_GetAnticipationTime(JPS_PushPropagationModelState handle);

/**
 * Write anticipation time of this agent.
 * @param handle of the Agent to access.
 * @param anticipationTime of this agent.
 */
JUPEDSIM_API void JPS_PushPropagationModelState_SetAnticipationTime(
    JPS_PushPropagationModelState handle,
    double anticipationTime);

/**
 * Read reaction time of this agent.
 * @param handle of the Agent to access.
 * @return reaction time of this agent
 */
JUPEDSIM_API double
JPS_PushPropagationModelState_GetReactionTime(JPS_PushPropagationModelState handle);

/**
 * Write reaction time of this agent.
 * @param handle of the Agent to access.
 * @param reactionTime of this agent.
 */
JUPEDSIM_API void JPS_PushPropagationModelState_SetReactionTime(
    JPS_PushPropagationModelState handle,
    double reactionTime);

/**
 * Write wall buffer distance of this agent.
 * @param handle of the Agent to access.
 * @param wallBufferDistance of this agent.
 */
JUPEDSIM_API void JPS_PushPropagationModelState_SetWallBufferDistance(
    JPS_PushPropagationModelState handle,
    double wallBufferDistance);

/**
 * Read wall buffer distance of this agent.
 * @param handle of the Agent to access.
 * @return wall buffer distance of this agent
 */
JUPEDSIM_API double
JPS_PushPropagationModelState_GetWallBufferDistance(JPS_PushPropagationModelState handle);

/**
 * Read e0 of this agent.
 * @param handle of the Agent to access.
 * @return e0 of this agent
 */
JUPEDSIM_API JPS_Point JPS_PushPropagationModelState_GetE0(JPS_PushPropagationModelState handle);

/**
 * Write e0 of this agent.
 * @param handle of the Agent to access.
 * @param e0 of this agent.
 */
JUPEDSIM_API void
JPS_PushPropagationModelState_SetE0(JPS_PushPropagationModelState handle, JPS_Point e0);

/**
 * Read time gap of this agent.
 * @param handle of the Agent to access.
 * @return time gap of this agent
 */
JUPEDSIM_API double JPS_PushPropagationModelState_GetTimeGap(JPS_PushPropagationModelState handle);

/**
 * Write time gap of this agent.
 * @param handle of the Agent to access.
 * @param time_gap of this agent.
 */
JUPEDSIM_API void JPS_PushPropagationModelState_SetTimeGap(JPS_PushPropagationModelState handle,
    double time_gap);

/**
 * Read v0 of this agent.
 * @param handle of the Agent to access.
 * @return v0 of this agent
 */
JUPEDSIM_API double JPS_PushPropagationModelState_GetV0(JPS_PushPropagationModelState handle);

/**
 * Write v0 of this agent.
 * @param handle of the Agent to access.
 * @param v0 of this agent.
 */
JUPEDSIM_API void
JPS_PushPropagationModelState_SetV0(JPS_PushPropagationModelState handle, double v0);

/**
 * Read radius of this agent.
 * @param handle of the Agent to access.
 * @return radius of this agent
 */
JUPEDSIM_API double JPS_PushPropagationModelState_GetRadius(JPS_PushPropagationModelState handle);

/**
 * Write radius of this agent in meters.
 * @param handle of the Agent to access.
 * @param radius (m) of this agent.
 */
JUPEDSIM_API void
JPS_PushPropagationModelState_SetRadius(JPS_PushPropagationModelState handle,
    double radius);

/**
 * Describes parameters of an Agent in Push Propagation Model
 */
typedef struct JPS_PushPropagationModelAgentParameters {
    /**
     * Position of the agent.
     * The position needs to inside the accessible area.
     */
    JPS_Point position{0, 0};
    /**
     * Defines the journey this agent will take use
     */
    JPS_JourneyId journeyId = 0;
    /**
     * Defines the current stage of its journey
     */
    JPS_StageId stageId = 0;

    /**
     * @param time_gap of the agents using this profile (T in the OV-function)
     */
    double time_gap = 1.06;
    /**
     * @param v0 of the agents using this profile(desired speed) double radius;
     */
    double v0 = 1.2;
    /**
     * @param radius of the agent in 'meters'
     */
    double radius = 0.15;

    /**
     *  Strength of the repulsion from neighbors
     */
    double strengthNeighborRepulsion{8.0};

    /**
     * Range of the repulsion from neighbors
     */
    double rangeNeighborRepulsion{0.1};

    /**
     * Wall buffer distance to geometry boundaries
     */
    double wallBufferDistance{0.1};

    /**
     * Anticipation time in seconds
     */
    double anticipationTime{0.5};

    /**
     * Reaction time in seconds
     */
    double reactionTime{0.1};

} JPS_PushPropagationModelAgentParameters;

#ifdef __cplusplus
}
#endif
