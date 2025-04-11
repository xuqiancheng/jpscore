// SPDX-License-Identifier: LGPL-3.0-or-later
#include "AnticipationVelocityModelUpdate.hpp"
#include "PushPropagationModelUpdate.hpp"
#include "CollisionFreeSpeedModelUpdate.hpp"
#include "CollisionFreeSpeedModelV2Update.hpp"
#include "GeneralizedCentrifugalForceModelUpdate.hpp"
#include "SocialForceModelUpdate.hpp"

#include <variant>

using OperationalModelUpdate = std::variant<
    GeneralizedCentrifugalForceModelUpdate,
    CollisionFreeSpeedModelUpdate,
    CollisionFreeSpeedModelV2Update,
    AnticipationVelocityModelUpdate,
    PushPropagationModelUpdate,
    SocialForceModelUpdate>;