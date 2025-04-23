// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once
#include "Point.hpp"

struct PushPropagationModelUpdate {
    Point position{};
    Point velocity{};
    Point orientation{};
    double simTime{};
};
